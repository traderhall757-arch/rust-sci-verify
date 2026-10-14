//! Running a target function over a dataset and matching expected failures.

use std::panic::{AssertUnwindSafe, catch_unwind};
use std::path::Path;
use std::sync::Mutex;

use crate::compare::{Verdict, classify};
use crate::data::{Dataset, Header, Point};

/// What the target crate returned at one point.
#[derive(Debug, Clone, PartialEq)]
pub enum Outcome {
    /// Output values, one list per output name in header order.
    Value(Vec<Vec<f64>>),
    /// The crate reported an error (e.g. `Err(..)` for an invalid input).
    Error(String),
    /// The crate panicked; always counted as crate wrong.
    Panic(String),
}

impl Outcome {
    /// A single scalar output.
    #[must_use]
    pub fn scalar(y: f64) -> Self {
        Self::Value(vec![vec![y]])
    }
}

/// A point where the crate is known to be wrong, documented in
/// `findings/<finding>.md`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Expected {
    /// Findings file stem, e.g. `001-puruspe-betai`.
    pub finding: &'static str,
    /// The point's inputs, flattened in header order.
    pub inputs: &'static [f64],
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
        // Panics inside the target crate are results, not test failures;
        // silence the panic hook so they do not flood the output. The hook is
        // global, so concurrent runs take turns to restore it correctly.
        static HOOK: Mutex<()> = Mutex::new(());
        let _turn = HOOK
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        let hook = std::panic::take_hook();
        std::panic::set_hook(Box::new(|_| {}));
        let points = data
            .points
            .into_iter()
            .map(|point| {
                let outcome = catch_unwind(AssertUnwindSafe(|| f(&point)))
                    .unwrap_or_else(|p| Outcome::Panic(panic_message(p.as_ref())));
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
        std::panic::set_hook(hook);
        FunctionResult {
            check: self,
            header: data.header,
            points,
        }
    }
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
            let failing = r.points.iter().any(|p| {
                p.verdict.is_crate_failure()
                    && same_inputs(&p.point.flat_inputs(&r.header), e.inputs)
            });
            if !failing {
                problems.push(format!(
                    "{name}{:?}: listed in {} but no longer fails; update the finding",
                    e.inputs, e.finding
                ));
            }
        }
        for p in r.points.iter().filter(|p| p.verdict.is_crate_failure()) {
            let inputs = p.point.flat_inputs(&r.header);
            if !r
                .check
                .expected
                .iter()
                .any(|e| same_inputs(&inputs, e.inputs))
            {
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
            inputs: &[0.5],
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
    fn same_inputs_compares_bits_and_length() {
        assert!(same_inputs(&[1.0, 2.0], &[1.0, 2.0]));
        assert!(!same_inputs(&[1.0], &[1.0, 2.0]));
        assert!(!same_inputs(&[0.0], &[-0.0]));
    }
}
