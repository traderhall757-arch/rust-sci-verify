//! Running a target function over a dataset and matching expected failures.

use std::panic::{AssertUnwindSafe, catch_unwind};
use std::path::Path;
use std::sync::{Arc, Mutex, mpsc};
use std::thread;
use std::time::{Duration, Instant};

use crate::compare::{Verdict, classify};
use crate::data::{Dataset, Header, Point};
use crate::isolated;

/// What the target crate returned at one point.
#[derive(Debug, Clone, PartialEq)]
pub enum Outcome {
    /// Output values, one list per output name in header order.
    Value(Vec<Vec<f64>>),
    /// The crate reported an error (e.g. `Err(..)` for an invalid input).
    Error(String),
    /// The crate panicked; always counted as crate wrong.
    Panic(String),
    /// The crate returned nothing within the time limit of
    /// [`Check::run_with_time_limit`] (e.g. an endless loop); always counted
    /// as crate wrong.
    Timeout(Duration),
    /// The process evaluating this point died, e.g. aborted on a stack
    /// overflow, in [`Check::run_isolated`]; always counted as crate wrong.
    Crash(String),
}

impl Outcome {
    /// A single scalar output.
    #[must_use]
    pub fn scalar(y: f64) -> Self {
        Self::Value(vec![vec![y]])
    }
}

/// Points where the crate is known to be wrong, documented in
/// `findings/<finding>.md`. Every listed point must fail and every failing
/// point must be listed, so the list is exact.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Expected {
    /// Findings file stem, e.g. `001-puruspe-betai`.
    pub finding: &'static str,
    /// Each point's inputs, flattened in header order.
    pub points: &'static [&'static [f64]],
}

impl Expected {
    /// The finding that lists `inputs`, if any.
    pub(crate) fn finding_for(expected: &[Self], inputs: &[f64]) -> Option<&'static str> {
        expected
            .iter()
            .find(|e| e.points.iter().any(|p| same_inputs(p, inputs)))
            .map(|e| e.finding)
    }
}

/// One function of a target crate checked against one dataset.
#[derive(Debug, Clone, Copy)]
pub struct Check {
    /// The function as users call it, e.g. `puruspe::erfcx`.
    pub function: &'static str,
    /// Dataset name under `reference/data/`.
    pub dataset: &'static str,
    /// Known crate-wrong points.
    pub expected: &'static [Expected],
}

/// The judgement of one point.
#[derive(Debug, Clone, PartialEq)]
pub struct PointResult {
    /// The reference point.
    pub point: Point,
    /// What the crate returned.
    pub outcome: Outcome,
    /// How it was judged.
    pub verdict: Verdict,
    /// The crate's relative error (NaN when not judged numerically).
    pub error: f64,
}

/// Every point of one [`Check`].
#[derive(Debug, Clone)]
pub struct FunctionResult {
    /// The check that was run.
    pub check: Check,
    /// The dataset's header.
    pub header: Header,
    /// Per-point results, in dataset order.
    pub points: Vec<PointResult>,
}

impl FunctionResult {
    /// Number of points with verdict `v`.
    #[must_use]
    pub fn count(&self, v: Verdict) -> usize {
        self.points.iter().filter(|p| p.verdict == v).count()
    }
}

fn panic_message(payload: &(dyn std::any::Any + Send)) -> String {
    payload
        .downcast_ref::<&str>()
        .map(ToString::to_string)
        .or_else(|| payload.downcast_ref::<String>().cloned())
        .unwrap_or_else(|| "non-string panic".into())
}

impl Check {
    /// Loads the dataset from the workspace and evaluates `f` at every point.
    ///
    /// # Panics
    ///
    /// If the dataset cannot be loaded.
    #[must_use]
    pub fn run(self, f: impl Fn(&Point) -> Outcome) -> FunctionResult {
        let data = Dataset::load(&crate::workspace_root(), self.dataset)
            .unwrap_or_else(|e| panic!("{}: {e}", self.dataset));
        self.run_on(data, f)
    }

    /// Evaluates `f` at every point of an already loaded dataset.
    #[must_use]
    pub fn run_on(self, data: Dataset, f: impl Fn(&Point) -> Outcome) -> FunctionResult {
        if isolated::in_child() {
            return self.judge(data, Vec::new());
        }
        let outcomes = silenced_panics(|| {
            data.points
                .iter()
                .map(|point| evaluate(&f, point))
                .collect()
        });
        self.judge(data, outcomes)
    }

    /// Like [`Check::run`], but a point that has not returned within `limit`
    /// gets [`Outcome::Timeout`], so a crate that loops forever cannot hang
    /// the test.
    ///
    /// Every point runs on its own thread, all at once, and the limit counts
    /// from the start of the run. A thread that has not finished is
    /// abandoned: it keeps running (and using a CPU, if it loops) until the
    /// test process exits, so put such a check last in its test.
    ///
    /// # Panics
    ///
    /// If the dataset cannot be loaded or a thread cannot be started.
    #[must_use]
    pub fn run_with_time_limit(
        self,
        limit: Duration,
        f: impl Fn(&Point) -> Outcome + Send + Sync + 'static,
    ) -> FunctionResult {
        let data = Dataset::load(&crate::workspace_root(), self.dataset)
            .unwrap_or_else(|e| panic!("{}: {e}", self.dataset));
        self.run_on_with_time_limit(data, limit, f)
    }

    /// [`Check::run_with_time_limit`] on an already loaded dataset.
    ///
    /// # Panics
    ///
    /// If a thread cannot be started.
    #[must_use]
    pub fn run_on_with_time_limit(
        self,
        data: Dataset,
        limit: Duration,
        f: impl Fn(&Point) -> Outcome + Send + Sync + 'static,
    ) -> FunctionResult {
        if isolated::in_child() {
            return self.judge(data, Vec::new());
        }
        let f = Arc::new(f);
        let outcomes = silenced_panics(|| {
            let deadline = Instant::now() + limit;
            let receivers: Vec<_> = data
                .points
                .iter()
                .map(|point| {
                    let (tx, rx) = mpsc::channel();
                    let (f, point) = (Arc::clone(&f), point.clone());
                    thread::Builder::new()
                        .spawn(move || {
                            // The receiver is gone only after a timeout.
                            let _ = tx.send(evaluate(f.as_ref(), &point));
                        })
                        .expect("start a thread");
                    rx
                })
                .collect();
            receivers
                .into_iter()
                .map(|rx| {
                    let left = deadline.saturating_duration_since(Instant::now());
                    rx.recv_timeout(left).unwrap_or(Outcome::Timeout(limit))
                })
                .collect()
        });
        self.judge(data, outcomes)
    }

    /// Like [`Check::run`], but the points are evaluated in child processes,
    /// so a call that ends the process (a stack overflow aborts it) gets
    /// [`Outcome::Crash`] and the run continues at the next point; a point
    /// with no result within `limit` of the previous one gets
    /// [`Outcome::Timeout`].
    ///
    /// Each child is the current test executable running only `test`, the
    /// calling test's full name (e.g. `verify::special_matches_reference`),
    /// and evaluates this check alone: every other check in that test
    /// returns no points there, and the child exits once this check is done.
    /// Put the check first in its test, so that the checks before it, though
    /// skipped, do not even load their data again in every child.
    ///
    /// # Panics
    ///
    /// If the dataset cannot be loaded, a child cannot be started, or a child
    /// never reaches this check (e.g. a wrong `test`).
    #[must_use]
    pub fn run_isolated(
        self,
        test: &str,
        limit: Duration,
        f: impl Fn(&Point) -> Outcome,
    ) -> FunctionResult {
        let data = Dataset::load(&crate::workspace_root(), self.dataset)
            .unwrap_or_else(|e| panic!("{}: {e}", self.dataset));
        self.run_on_isolated(data, test, limit, f)
    }

    /// [`Check::run_isolated`] on an already loaded dataset; the child must
    /// build the same dataset.
    ///
    /// # Panics
    ///
    /// As [`Check::run_isolated`].
    #[must_use]
    pub fn run_on_isolated(
        self,
        data: Dataset,
        test: &str,
        limit: Duration,
        f: impl Fn(&Point) -> Outcome,
    ) -> FunctionResult {
        if isolated::is_child_for(self) {
            isolated::serve(&data, &f);
        }
        if isolated::in_child() {
            return self.judge(data, Vec::new());
        }
        let outcomes = isolated::run(self, &data, test, limit);
        self.judge(data, outcomes)
    }

    fn judge(self, data: Dataset, outcomes: Vec<Outcome>) -> FunctionResult {
        let points = data
            .points
            .into_iter()
            .zip(outcomes)
            .map(|(point, outcome)| {
                let (verdict, error) = classify(
                    &outcome,
                    &point,
                    &data.header.outputs,
                    data.header.tolerance.rel,
                );
                PointResult {
                    point,
                    outcome,
                    verdict,
                    error,
                }
            })
            .collect();
        FunctionResult {
            check: self,
            header: data.header,
            points,
        }
    }
}

/// `f(point)`, with a panic turned into [`Outcome::Panic`].
pub(crate) fn evaluate(f: &impl Fn(&Point) -> Outcome, point: &Point) -> Outcome {
    catch_unwind(AssertUnwindSafe(|| f(point)))
        .unwrap_or_else(|p| Outcome::Panic(panic_message(p.as_ref())))
}

/// Runs `body` with the panic hook silenced. Panics inside the target crate
/// are results, not test failures, and must not flood the output. The hook is
/// global, so concurrent runs take turns to restore it correctly.
fn silenced_panics(body: impl FnOnce() -> Vec<Outcome>) -> Vec<Outcome> {
    static HOOK: Mutex<()> = Mutex::new(());
    let _turn = HOOK
        .lock()
        .unwrap_or_else(std::sync::PoisonError::into_inner);
    let hook = std::panic::take_hook();
    std::panic::set_hook(Box::new(|_| {}));
    let outcomes = body();
    std::panic::set_hook(hook);
    outcomes
}

pub(crate) fn same_inputs(a: &[f64], b: &[f64]) -> bool {
    a.len() == b.len() && a.iter().zip(b).all(|(x, y)| x.to_bits() == y.to_bits())
}

/// Differences between the crate-wrong points and the declared expected
/// failures, plus expected entries whose findings file does not exist.
/// Empty means the run matches what is documented.
#[must_use]
pub fn expectation_problems(results: &[FunctionResult], findings_dir: &Path) -> Vec<String> {
    let mut problems = Vec::new();
    for r in results {
        let name = r.check.function;
        for e in r.check.expected {
            if !findings_dir.join(format!("{}.md", e.finding)).is_file() {
                problems.push(format!(
                    "{name}: expected failure cites missing findings/{}.md",
                    e.finding
                ));
            }
            for inputs in e.points {
                let failing = r.points.iter().any(|p| {
                    p.verdict.is_crate_failure()
                        && same_inputs(&p.point.flat_inputs(&r.header), inputs)
                });
                if !failing {
                    problems.push(format!(
                        "{name}{inputs:?}: listed in {} but no longer fails; update the finding",
                        e.finding
                    ));
                }
            }
        }
        for p in r.points.iter().filter(|p| p.verdict.is_crate_failure()) {
            let inputs = p.point.flat_inputs(&r.header);
            if Expected::finding_for(r.check.expected, &inputs).is_none() {
                problems.push(format!(
                    "{name}{inputs:?}: unexpected {:?} (crate {:?}, mpmath {:?})",
                    p.verdict, p.outcome, p.point.mpmath
                ));
            }
        }
    }
    problems
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data::tests::sample;

    fn data() -> Dataset {
        Dataset::parse(&sample()).unwrap()
    }

    const NONE: Check = Check {
        function: "t::f",
        dataset: "f",
        expected: &[],
    };

    #[test]
    fn runs_values_errors_and_panics() {
        let r = NONE.run_on(data(), |p| {
            if p.scalar("x") < 0.0 {
                panic!("negative")
            } else {
                Outcome::scalar(0.1)
            }
        });
        assert_eq!(r.points[0].verdict, Verdict::Pass);
        assert_eq!(r.points[1].outcome, Outcome::Panic("negative".into()));
        assert_eq!(r.count(Verdict::CrateWrong), 1);
        let r = NONE.run_on(data(), |_| std::panic::panic_any(format!("owned {}", 1)));
        assert_eq!(r.points[0].outcome, Outcome::Panic("owned 1".into()));
        let r = NONE.run_on(data(), |_| std::panic::panic_any(7_i32));
        assert_eq!(
            r.points[0].outcome,
            Outcome::Panic("non-string panic".into())
        );
    }

    #[test]
    fn time_limit_turns_endless_calls_into_timeouts() {
        let limit = Duration::from_millis(500);
        let r = NONE.run_on_with_time_limit(data(), limit, |p| {
            if p.scalar("x") < 0.0 {
                // Never returns, without using a CPU.
                loop {
                    thread::park();
                }
            }
            // Slow but within the limit: still a result.
            thread::sleep(Duration::from_millis(50));
            Outcome::scalar(0.1)
        });
        assert_eq!(r.points[0].verdict, Verdict::Pass);
        assert_eq!(r.points[1].outcome, Outcome::Timeout(limit));
        // Point 1 is an overflow point: a timeout is never an accepted value.
        assert_eq!(r.points[1].verdict, Verdict::CrateWrong);
        // Panics are still results, and points keep their dataset order.
        let r = NONE.run_on_with_time_limit(data(), limit, |p| {
            assert!(p.scalar("x") > 0.0, "negative");
            Outcome::scalar(0.1)
        });
        assert_eq!(r.points[0].outcome, Outcome::scalar(0.1));
        assert_eq!(r.points[1].outcome, Outcome::Panic("negative".into()));
    }

    #[test]
    fn run_with_time_limit_loads_committed_data() {
        let check = Check {
            function: "t::erfcx",
            dataset: "erfcx",
            expected: &[],
        };
        let r = check.run_with_time_limit(Duration::from_secs(60), |_| Outcome::scalar(0.0));
        assert_eq!(r.points.len(), r.header.points);
        assert!(r.points.iter().all(|p| p.outcome == Outcome::scalar(0.0)));
    }

    #[test]
    fn run_loads_committed_data() {
        let check = Check {
            function: "t::erfcx",
            dataset: "erfcx",
            expected: &[],
        };
        let r = check.run(|_| Outcome::scalar(0.0));
        assert_eq!(r.points.len(), r.header.points);
    }

    #[test]
    #[should_panic(expected = "no-such-dataset")]
    fn run_panics_on_missing_dataset() {
        let check = Check {
            function: "t::f",
            dataset: "no-such-dataset",
            expected: &[],
        };
        let _ = check.run(|_| Outcome::scalar(0.0));
    }

    fn findings_dir(with: &[&str]) -> std::path::PathBuf {
        let dir = std::env::temp_dir().join(format!(
            "rsv-findings-{}-{}",
            std::process::id(),
            with.len()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        for f in with {
            std::fs::write(dir.join(format!("{f}.md")), "x").unwrap();
        }
        dir
    }

    #[test]
    fn expectations_must_match_failures_exactly_and_cite_findings() {
        const EXP: &[Expected] = &[Expected {
            finding: "001-t-f",
            points: &[&[0.5]],
        }];
        let wrong = |_: &Point| Outcome::scalar(0.2);
        let dir = findings_dir(&["001-t-f"]);

        // Unexpected failures: both points fail, nothing expected.
        let r = NONE.run_on(data(), wrong);
        let p = expectation_problems(&[r], &dir);
        assert_eq!(p.len(), 2);
        assert!(p[0].contains("unexpected CrateWrong"));

        // Point 0 expected and failing; point 1 (overflow, finite) still unexpected.
        let check = Check {
            expected: EXP,
            ..NONE
        };
        let p = expectation_problems(&[check.run_on(data(), wrong)], &dir);
        assert_eq!(p.len(), 1);
        assert!(p[0].contains("[-30.0]"));

        // Expected failure that now passes.
        let fixed = |p: &Point| {
            Outcome::scalar(if p.scalar("x") > 0.0 {
                0.1
            } else {
                f64::INFINITY
            })
        };
        let p = expectation_problems(&[check.run_on(data(), fixed)], &dir);
        assert_eq!(p.len(), 1);
        assert!(p[0].contains("no longer fails"));

        // Expected failure citing a findings file that does not exist.
        let empty = findings_dir(&[]);
        let all_but_0 = |p: &Point| {
            Outcome::scalar(if p.scalar("x") > 0.0 {
                0.2
            } else {
                f64::INFINITY
            })
        };
        let p = expectation_problems(&[check.run_on(data(), all_but_0)], &empty);
        assert_eq!(
            p,
            vec!["t::f: expected failure cites missing findings/001-t-f.md".to_string()]
        );
        assert!(expectation_problems(&[check.run_on(data(), all_but_0)], &dir).is_empty());
        std::fs::remove_dir_all(dir).unwrap();
        std::fs::remove_dir_all(empty).unwrap();
    }

    #[test]
    fn one_finding_lists_many_points() {
        const EXP: &[Expected] = &[
            Expected {
                finding: "001-t-f",
                points: &[&[1.0], &[0.5], &[-30.0]],
            },
            Expected {
                finding: "002-t-f",
                points: &[&[2.0]],
            },
        ];
        assert_eq!(Expected::finding_for(EXP, &[0.5]), Some("001-t-f"));
        assert_eq!(Expected::finding_for(EXP, &[-30.0]), Some("001-t-f"));
        assert_eq!(Expected::finding_for(EXP, &[2.0]), Some("002-t-f"));
        assert_eq!(Expected::finding_for(EXP, &[3.0]), None);
        let dir = findings_dir(&["001-t-f", "002-t-f", "x"]);
        let check = Check {
            expected: EXP,
            ..NONE
        };
        // Both dataset points fail and are listed; [1.0] and [2.0] are listed
        // but are not in the dataset, so they "no longer fail".
        let p = expectation_problems(&[check.run_on(data(), |_| Outcome::scalar(0.2))], &dir);
        assert_eq!(
            p,
            vec![
                "t::f[1.0]: listed in 001-t-f but no longer fails; update the finding".to_string(),
                "t::f[2.0]: listed in 002-t-f but no longer fails; update the finding".to_string(),
            ]
        );
        std::fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn same_inputs_compares_bits_and_length() {
        assert!(same_inputs(&[1.0, 2.0], &[1.0, 2.0]));
        assert!(!same_inputs(&[1.0], &[1.0, 2.0]));
        assert!(!same_inputs(&[0.0], &[-0.0]));
    }
}
