//! Running a check's points in child processes, so that a crate call that
//! ends the process (a stack overflow aborts it; so does `abort`) or never
//! returns becomes a result instead of ending or hanging the test.
//!
//! The parent starts the current test executable again, running only the
//! named test, with two environment variables: which check to evaluate and
//! from which point. In that child every other check is skipped
//! ([`in_child`]); the isolated check evaluates its points in order, writes
//! one line per point to stdout and exits. When the child dies or a point
//! takes too long, the parent records [`Outcome::Crash`] or
//! [`Outcome::Timeout`] for that point and starts a new child at the next.

use std::io::{BufRead, BufReader, Read, Write};
use std::process::{Child, Command, ExitStatus, Stdio};
use std::sync::mpsc::{self, Receiver, RecvTimeoutError};
use std::thread;
use std::time::Duration;
use std::{env, io, process};

use serde::{Deserialize, Serialize};

use crate::check::{Check, Outcome, evaluate};
use crate::data::{Dataset, Point};

/// Which check a child evaluates: `function<TAB>dataset`.
const CHECK_VAR: &str = "RSV_ISOLATED_CHECK";
/// The index of the first point a child evaluates.
const START_VAR: &str = "RSV_ISOLATED_START";
/// Written by a child once it has found its check: `RSV-ISOLATED-READY <start>`.
const READY: &str = "RSV-ISOLATED-READY ";
/// Written by a child after its last point.
const DONE: &str = "RSV-ISOLATED-DONE";
/// Starts each result line: `RSV-ISOLATED-POINT <index> <json>`.
const POINT: &str = "RSV-ISOLATED-POINT ";

/// True in a child process started by [`run`]. There, every check other than
/// the isolated one is skipped.
pub(crate) fn in_child() -> bool {
    env::var_os(CHECK_VAR).is_some()
}

fn key(check: Check) -> String {
    format!("{}\t{}", check.function, check.dataset)
}

/// True in the child process that evaluates `check`.
pub(crate) fn is_child_for(check: Check) -> bool {
    env::var(CHECK_VAR).is_ok_and(|k| k == key(check))
}

/// An [`Outcome`] as a child sends it. Values travel as bit patterns, so
/// NaN payloads, signed zeros and infinities arrive exactly.
#[derive(Serialize, Deserialize)]
enum Wire {
    Value(Vec<Vec<u64>>),
    Error(String),
    Panic(String),
}

impl Wire {
    fn from_outcome(o: Outcome) -> Self {
        match o {
            Outcome::Value(v) => Self::Value(
                v.iter()
                    .map(|xs| xs.iter().map(|x| x.to_bits()).collect())
                    .collect(),
            ),
            Outcome::Error(e) => Self::Error(e),
            // A child evaluates directly: no time limit, and a crash ends it.
            Outcome::Panic(e) | Outcome::Crash(e) => Self::Panic(e),
            Outcome::Timeout(_) => Self::Panic(String::new()),
        }
    }

    fn into_outcome(self) -> Outcome {
        match self {
            Self::Value(v) => Outcome::Value(
                v.iter()
                    .map(|xs| xs.iter().map(|&x| f64::from_bits(x)).collect())
                    .collect(),
            ),
            Self::Error(e) => Outcome::Error(e),
            Self::Panic(e) => Outcome::Panic(e),
        }
    }
}

/// Child side: evaluates `f` at every point from `RSV_ISOLATED_START` on,
/// writes one line per point, and exits the process.
pub(crate) fn serve(data: &Dataset, f: &impl Fn(&Point) -> Outcome) -> ! {
    let start: usize = env::var(START_VAR)
        .ok()
        .and_then(|s| s.parse().ok())
        .expect("RSV_ISOLATED_START is an index");
    // Panics inside the crate are results; this process does nothing else.
    std::panic::set_hook(Box::new(|_| {}));
    let mut out = io::stdout().lock();
    let mut send = |line: String| {
        writeln!(out, "{line}")
            .and_then(|()| out.flush())
            .expect("write to the parent");
    };
    // libtest may have started a line ("test name ... ") without ending it.
    send(format!("\n{READY}{start}"));
    for (i, point) in data.points.iter().enumerate().skip(start) {
        let wire = Wire::from_outcome(evaluate(f, point));
        send(format!(
            "{POINT}{i} {}",
            serde_json::to_string(&wire).expect("encode an outcome")
        ));
    }
    send(DONE.to_string());
    process::exit(0)
}

/// Parent side: the outcome at every point of `data`, each evaluated in a
/// child process running test `test` of the current executable. A point
/// with no result within `limit` after the previous one gets
/// [`Outcome::Timeout`]; a point during which the child dies gets
/// [`Outcome::Crash`].
pub(crate) fn run(check: Check, data: &Dataset, test: &str, limit: Duration) -> Vec<Outcome> {
    // A child that started children of its own could do so without end.
    assert!(
        !in_child(),
        "an isolated child reached another isolated check's parent side"
    );
    let exe = env::current_exe().expect("the test executable's path");
    let mut outcomes = Vec::with_capacity(data.points.len());
    loop {
        let mut child = Command::new(&exe)
            .args([test, "--exact", "--nocapture", "--test-threads=1"])
            .env(CHECK_VAR, key(check))
            .env(START_VAR, outcomes.len().to_string())
            .stdin(Stdio::null())
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .expect("start a child process");
        let lines = lines_of(child.stdout.take().expect("piped stdout"));
        let mut stderr = child.stderr.take().expect("piped stderr");
        let stderr = thread::spawn(move || {
            let mut text = String::new();
            // Whatever could be read is enough to describe a crash.
            let _ = stderr.read_to_string(&mut text);
            text
        });
        let end = follow(&lines, &mut child, limit, &mut outcomes);
        let stderr = stderr.join().expect("read the child's stderr");
        match end {
            End::Done => break,
            End::TimedOut => {}
            End::Died(status) => outcomes.push(Outcome::Crash(crash_reason(status, &stderr))),
            End::NotReached(status) => panic!(
                "{test}: the child process did not reach check {} on {} \
                 (is the test name right?); {status}; stderr:\n{stderr}",
                check.function, check.dataset
            ),
        }
    }
    assert_eq!(outcomes.len(), data.points.len(), "one outcome per point");
    outcomes
}

/// How a child process ended.
enum End {
    /// It evaluated every point from its start on.
    Done,
    /// A point took longer than the limit; the child was killed.
    TimedOut,
    /// It died during a point.
    Died(ExitStatus),
    /// It ended before reaching the check: no such test, or the test failed
    /// before the check.
    NotReached(ExitStatus),
}

/// Collects a child's results into `outcomes` until it ends. Before READY the
/// child only starts and loads data, so there is no time limit; after it,
/// each point must arrive within `limit` of the previous one, or it gets
/// [`Outcome::Timeout`] and the child is killed.
fn follow(
    lines: &Receiver<String>,
    child: &mut Child,
    limit: Duration,
    outcomes: &mut Vec<Outcome>,
) -> End {
    let mut ready = false;
    loop {
        let line = if ready {
            lines.recv_timeout(limit)
        } else {
            lines.recv().map_err(|_| RecvTimeoutError::Disconnected)
        };
        match line {
            Ok(line) => {
                // Prefix matches throughout: libtest's own lines never start
                // with these markers.
                if line.starts_with(DONE) {
                    child.wait().expect("wait for the child");
                    return End::Done;
                } else if let Some(rest) = line.strip_prefix(POINT) {
                    outcomes.push(parse_point(rest, outcomes.len()));
                } else if let Some(start) = line.strip_prefix(READY) {
                    assert_eq!(
                        start.parse::<usize>().ok(),
                        Some(outcomes.len()),
                        "the child starts at the next point"
                    );
                    ready = true;
                }
            }
            Err(RecvTimeoutError::Timeout) => {
                // Killing fails only if the child has already exited.
                let _ = child.kill();
                child.wait().expect("wait for the child");
                outcomes.push(Outcome::Timeout(limit));
                return End::TimedOut;
            }
            Err(RecvTimeoutError::Disconnected) => {
                let status = child.wait().expect("wait for the child");
                return if ready {
                    End::Died(status)
                } else {
                    End::NotReached(status)
                };
            }
        }
    }
}

/// The lines of `out`, read on a thread; the channel closes at end of file.
fn lines_of(out: impl Read + Send + 'static) -> Receiver<String> {
    let (tx, rx) = mpsc::channel();
    thread::spawn(move || {
        for line in BufReader::new(out).lines().map_while(Result::ok) {
            if tx.send(line).is_err() {
                break;
            }
        }
    });
    rx
}

fn parse_point(rest: &str, expected_index: usize) -> Outcome {
    let (index, json) = rest.split_once(' ').expect("index and outcome");
    assert_eq!(
        index.parse::<usize>().expect("a point index"),
        expected_index,
        "points arrive in order"
    );
    serde_json::from_str::<Wire>(json)
        .expect("decode an outcome")
        .into_outcome()
}

/// How a child that died during a point is described in the report.
fn crash_reason(status: ExitStatus, stderr: &str) -> String {
    if stderr.contains("has overflowed its stack") {
        "stack overflow".into()
    } else {
        format!("process ended ({status})")
    }
}

#[cfg(test)]
mod tests {
    use std::hint::black_box;
    use std::time::Instant;

    use super::*;
    use crate::compare::Verdict;
    use crate::data::tests::sample;

    fn data() -> Dataset {
        Dataset::parse(&sample()).unwrap()
    }

    fn named(function: &'static str) -> Check {
        Check {
            function,
            dataset: "f",
            expected: &[],
        }
    }

    /// Recurses `n` times with a real frame each time (no tail calls).
    fn deep(n: u64) -> u64 {
        let frame = black_box([n; 64]);
        if n == 0 { 0 } else { deep(n - 1) + frame[0] }
    }

    const LIMIT: Duration = Duration::from_secs(60);

    /// False in a child. Checks are skipped there, so assertions about any
    /// check but the test's last one run only in the parent.
    fn parent() -> bool {
        env::var_os(CHECK_VAR).is_none()
    }

    #[test]
    fn a_child_skips_every_other_check() {
        // Fails the test if a child evaluates it.
        let plain = named("t::plain").run_on(data(), |_| {
            assert!(parent(), "a child evaluated a plain check");
            Outcome::scalar(0.1)
        });
        assert_eq!(plain.points.len(), if parent() { 2 } else { 0 });
        let test = "isolated::tests::a_child_skips_every_other_check";
        // The next check's children must skip this one, or they would
        // report its values.
        let first =
            named("t::first").run_on_isolated(data(), test, LIMIT, |_| Outcome::scalar(0.1));
        if parent() {
            assert_eq!(first.points[1].outcome, Outcome::scalar(0.1));
        }
        let r = named("t::overflow").run_on_isolated(data(), test, LIMIT, |p| {
            Outcome::scalar(deep(if p.scalar("x") < 0.0 { u64::MAX } else { 0 }) as f64)
        });
        assert_eq!(r.points[0].outcome, Outcome::scalar(0.0));
        assert_eq!(r.points[1].outcome, Outcome::Crash("stack overflow".into()));
        assert_eq!(r.points[1].verdict, Verdict::CrateWrong);
    }

    #[test]
    fn after_a_crash_the_next_child_starts_at_the_next_point() {
        let r = named("t::abort").run_on_isolated(
            data(),
            "isolated::tests::after_a_crash_the_next_child_starts_at_the_next_point",
            LIMIT,
            |p| {
                if p.scalar("x") > 0.0 {
                    process::abort()
                }
                Outcome::scalar(f64::INFINITY)
            },
        );
        let Outcome::Crash(reason) = &r.points[0].outcome else {
            panic!("{:?}", r.points[0].outcome)
        };
        assert!(reason.starts_with("process ended (signal: 6"), "{reason}");
        assert_eq!(r.points[1].outcome, Outcome::scalar(f64::INFINITY));
        assert_eq!(r.points[1].verdict, Verdict::Convention);
    }

    #[test]
    fn a_point_that_never_returns_times_out() {
        let limit = Duration::from_millis(500);
        let started = Instant::now();
        let r = named("t::endless").run_on_isolated(
            data(),
            "isolated::tests::a_point_that_never_returns_times_out",
            limit,
            |p| {
                while p.scalar("x") > 0.0 {
                    thread::park();
                }
                Outcome::Error("e".into())
            },
        );
        assert_eq!(r.points[0].outcome, Outcome::Timeout(limit));
        assert_eq!(r.points[1].outcome, Outcome::Error("e".into()));
        assert!(started.elapsed() < LIMIT, "the endless child was killed");
    }

    #[test]
    fn values_and_panics_arrive_exactly() {
        let nan = f64::from_bits(0x7ff8_0000_0000_0123);
        let r = named("t::exact").run_on_isolated(
            data(),
            "isolated::tests::values_and_panics_arrive_exactly",
            LIMIT,
            move |p| {
                assert!(p.scalar("x") > 0.0, "negative");
                Outcome::Value(vec![vec![-0.0, nan], vec![]])
            },
        );
        let Outcome::Value(v) = &r.points[0].outcome else {
            panic!("{:?}", r.points[0].outcome)
        };
        let bits: Vec<u64> = v[0].iter().map(|x| x.to_bits()).collect();
        assert_eq!(bits, [(-0.0f64).to_bits(), nan.to_bits()]);
        assert!(v[1].is_empty());
        assert_eq!(r.points[1].outcome, Outcome::Panic("negative".into()));
    }

    #[test]
    fn run_isolated_loads_committed_data() {
        let check = Check {
            function: "t::erfcx",
            dataset: "erfcx",
            expected: &[],
        };
        let r = check.run_isolated(
            "isolated::tests::run_isolated_loads_committed_data",
            LIMIT,
            |_| Outcome::scalar(0.0),
        );
        assert_eq!(r.points.len(), r.header.points);
        assert!(r.points.iter().all(|p| p.outcome == Outcome::scalar(0.0)));
    }

    #[test]
    #[should_panic(expected = "is the test name right?")]
    fn a_wrong_test_name_is_reported() {
        let _ = named("t::f")
            .run_on_isolated(data(), "no::such::test", LIMIT, |_| Outcome::scalar(0.0));
    }

    #[test]
    fn wire_keeps_errors_and_reports_crashes_as_panics() {
        let back = |o| Wire::from_outcome(o).into_outcome();
        assert_eq!(back(Outcome::Error("e".into())), Outcome::Error("e".into()));
        assert_eq!(back(Outcome::Crash("c".into())), Outcome::Panic("c".into()));
        assert_eq!(back(Outcome::Timeout(LIMIT)), Outcome::Panic(String::new()));
    }

    #[test]
    fn crash_reasons() {
        use std::os::unix::process::ExitStatusExt;
        let abort = ExitStatus::from_raw(6);
        assert_eq!(
            crash_reason(abort, "thread 'x' has overflowed its stack\n"),
            "stack overflow"
        );
        assert!(crash_reason(abort, "").starts_with("process ended (signal: 6"));
    }
}
