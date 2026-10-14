use harness::{Check, Expected, Outcome, finish};

use crate::expected;
use statrs::distribution::{ContinuousCDF, Normal};
use statrs::function::{erf, gamma};

/// The standard normal distribution (the reference data has mean 0, sd 1).
fn standard_normal() -> Normal {
    Normal::new(0.0, 1.0).expect("valid parameters")
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
    ];
    finish("statrs", &results);
}
