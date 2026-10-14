use std::time::Duration;

use harness::{Check, Expected, Outcome, finish};

use crate::expected;
use statrs::distribution::{ContinuousCDF, Normal, StudentsT};
use statrs::function::{erf, gamma};

/// The standard normal distribution (the reference data has mean 0, sd 1).
fn standard_normal() -> Normal {
    Normal::new(0.0, 1.0).expect("valid parameters")
}

/// How long one point of a time-limited check may take. A correct answer
/// takes microseconds; the limit only has to tell an endless loop apart on a
/// busy CI machine.
const TIME_LIMIT: Duration = Duration::from_secs(5);

/// Student's t with `df` degrees of freedom (the reference data has location
/// 0, scale 1).
fn students_t(df: f64) -> StudentsT {
    StudentsT::new(0.0, 1.0, df).expect("valid parameters")
}

#[test]
fn statrs_matches_reference() {
    let results = [
        Check {
            function: "statrs::function::gamma::gamma",
            dataset: "gamma",
            expected: &[Expected {
                finding: "003-statrs-gamma",
                points: expected::GAMMA,
            }],
        }
        .run(|p| Outcome::scalar(gamma::gamma(p.scalar("x")))),
        Check {
            function: "statrs::function::erf::erfc",
            dataset: "erfc",
            expected: &[Expected {
                finding: "004-statrs-erfc",
                points: expected::ERFC,
            }],
        }
        .run(|p| Outcome::scalar(erf::erfc(p.scalar("x")))),
        Check {
            function: "statrs::distribution::Normal::cdf",
            dataset: "norm_cdf",
            expected: &[Expected {
                finding: "009-statrs-normal",
                points: expected::NORMAL_CDF,
            }],
        }
        .run(|p| Outcome::scalar(standard_normal().cdf(p.scalar("x")))),
        Check {
            function: "statrs::distribution::Normal::sf",
            dataset: "norm_sf",
            expected: &[Expected {
                finding: "009-statrs-normal",
                points: expected::NORMAL_SF,
            }],
        }
        .run(|p| Outcome::scalar(standard_normal().sf(p.scalar("x")))),
        Check {
            function: "statrs::distribution::Normal::inverse_cdf",
            dataset: "norm_ppf",
            expected: &[],
        }
        .run(|p| Outcome::scalar(standard_normal().inverse_cdf(p.scalar("p")))),
        Check {
            function: "statrs::distribution::StudentsT::cdf",
            dataset: "t_cdf",
            expected: &[Expected {
                finding: "010-statrs-students-t",
                points: expected::T_CDF,
            }],
        }
        .run(|p| Outcome::scalar(students_t(p.scalar("df")).cdf(p.scalar("x")))),
        Check {
            function: "statrs::distribution::StudentsT::sf",
            dataset: "t_sf",
            expected: &[Expected {
                finding: "010-statrs-students-t",
                points: expected::T_SF,
            }],
        }
        .run(|p| Outcome::scalar(students_t(p.scalar("df")).sf(p.scalar("x")))),
        // Last: it loops forever at some points (finding 011), and the
        // abandoned threads keep a CPU busy until the test ends.
        Check {
            function: "statrs::distribution::StudentsT::inverse_cdf",
            dataset: "t_ppf",
            expected: &[Expected {
                finding: "011-statrs-students-t-inverse-cdf",
                points: expected::T_PPF,
            }],
        }
        .run_with_time_limit(TIME_LIMIT, |p| {
            Outcome::scalar(students_t(p.scalar("df")).inverse_cdf(p.scalar("p")))
        }),
    ];
    finish("statrs", &results);
}
