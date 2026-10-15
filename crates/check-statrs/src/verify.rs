use std::time::Duration;

use harness::{Check, Dataset, Expected, Outcome, finish, workspace_root};

use crate::expected;
use statrs::distribution::{
    Beta, ChiSquared, ContinuousCDF, Erlang, Exp, FisherSnedecor, Gamma, Normal, StudentsT,
};
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

/// The gamma distribution with shape `a` and rate 1 (the reference data's
/// parametrisation).
fn gamma_dist(a: f64) -> Gamma {
    Gamma::new(a, 1.0).expect("valid parameters")
}

/// Chi-squared with `df = 2a` degrees of freedom: X is chi-squared(2a) exactly
/// when X / 2 is gamma(a, rate 1), so the gamma datasets check it at 2x (cdf,
/// sf) and through half its quantile (ppf). Doubling is exact for every grid
/// value (none is above 1e300); halving is exact except for the lowest bit of
/// a subnormal quantile, where errors are judged relative to 2.2e-308 anyway.
fn chi_squared(a: f64) -> ChiSquared {
    ChiSquared::new(2.0 * a).expect("valid parameters")
}

/// The exponential distribution with rate 1 (the reference data's).
fn exponential() -> Exp {
    Exp::new(1.0).expect("valid parameters")
}

/// The points of a gamma dataset whose shape `a` is a whole number, the shapes
/// an Erlang distribution can have (1, 2, 3, 5, 10, 20, 50, 100, 1e3, 1e4, 1e5
/// and 1e6 on these grids). The report counts only these points, and its
/// grid line says so.
fn integer_shapes(name: &str) -> Dataset {
    let mut data = Dataset::load(&workspace_root(), name).unwrap_or_else(|e| panic!("{name}: {e}"));
    data.points.retain(|p| p.scalar("a").fract() == 0.0);
    data.header.points = data.points.len();
    data.header.grid = format!(
        "the whole-number shapes of the {name} grid (a = 1, 2, 3, 5, 10, 20, 50, 100, 1e3, \
         1e4, 1e5, 1e6), {} of its points. The full grid: {}",
        data.points.len(),
        data.header.grid
    );
    data
}

/// The Erlang distribution with shape `a` (a whole number) and rate 1; the
/// gamma distribution with the same parameters.
#[expect(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    reason = "`a` is a whole number from 1 to 1e6 (`integer_shapes`)"
)]
fn erlang(a: f64) -> Erlang {
    Erlang::new(a as u64, 1.0).expect("valid parameters")
}

/// The beta distribution with shapes `a` and `b` (on [0, 1], as in the
/// reference data).
fn beta_dist(a: f64, b: f64) -> Beta {
    Beta::new(a, b).expect("valid parameters")
}

/// The F distribution with `d1` and `d2` degrees of freedom.
fn fisher_snedecor(d1: f64, d2: f64) -> FisherSnedecor {
    FisherSnedecor::new(d1, d2).expect("valid parameters")
}

/// Evaluates `f` at the point's degrees of freedom and argument `x` (or `p`).
fn f_dist(p: &harness::Point, arg: &str, f: impl Fn(&FisherSnedecor, f64) -> f64) -> Outcome {
    Outcome::scalar(f(
        &fisher_snedecor(p.scalar("d1"), p.scalar("d2")),
        p.scalar(arg),
    ))
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
        Check {
            function: "statrs::distribution::Gamma::cdf",
            dataset: "gamma_cdf",
            expected: &[Expected {
                finding: "012-statrs-gamma-cdf",
                points: expected::GAMMA_CDF,
            }],
        }
        .run(|p| Outcome::scalar(gamma_dist(p.scalar("a")).cdf(p.scalar("x")))),
        Check {
            function: "statrs::distribution::Gamma::sf",
            dataset: "gamma_sf",
            expected: &[Expected {
                finding: "012-statrs-gamma-cdf",
                points: expected::GAMMA_SF,
            }],
        }
        .run(|p| Outcome::scalar(gamma_dist(p.scalar("a")).sf(p.scalar("x")))),
        Check {
            function: "statrs::distribution::ChiSquared::cdf",
            dataset: "gamma_cdf",
            expected: &[Expected {
                finding: "012-statrs-gamma-cdf",
                points: expected::GAMMA_CDF,
            }],
        }
        .run(|p| Outcome::scalar(chi_squared(p.scalar("a")).cdf(2.0 * p.scalar("x")))),
        Check {
            function: "statrs::distribution::ChiSquared::sf",
            dataset: "gamma_sf",
            expected: &[Expected {
                finding: "012-statrs-gamma-cdf",
                points: expected::GAMMA_SF,
            }],
        }
        .run(|p| Outcome::scalar(chi_squared(p.scalar("a")).sf(2.0 * p.scalar("x")))),
        // Time-limited although no point is known to hang: the series and
        // continued fraction inside `cdf` have no iteration limit. Before
        // `t_ppf`, whose abandoned threads would compete for the CPU.
        Check {
            function: "statrs::distribution::Gamma::inverse_cdf",
            dataset: "gamma_ppf",
            expected: &[
                Expected {
                    finding: "013-statrs-gamma-inverse-cdf",
                    points: expected::GAMMA_PPF_WHOLE,
                },
                Expected {
                    finding: "013-statrs-gamma-inverse-cdf",
                    points: expected::GAMMA_PPF_FRACTIONAL,
                },
            ],
        }
        .run_with_time_limit(TIME_LIMIT, |p| {
            Outcome::scalar(gamma_dist(p.scalar("a")).inverse_cdf(p.scalar("p")))
        }),
        Check {
            function: "statrs::distribution::ChiSquared::inverse_cdf",
            dataset: "gamma_ppf",
            expected: &[Expected {
                finding: "013-statrs-gamma-inverse-cdf",
                points: expected::CHI2_PPF,
            }],
        }
        .run_with_time_limit(TIME_LIMIT, |p| {
            Outcome::scalar(chi_squared(p.scalar("a")).inverse_cdf(p.scalar("p")) / 2.0)
        }),
        Check {
            function: "statrs::distribution::Exp::cdf",
            dataset: "expon_cdf",
            expected: &[Expected {
                finding: "023-statrs-exp-cdf",
                points: expected::EXP_CDF,
            }],
        }
        .run(|p| Outcome::scalar(exponential().cdf(p.scalar("x")))),
        Check {
            function: "statrs::distribution::Exp::sf",
            dataset: "expon_sf",
            expected: &[],
        }
        .run(|p| Outcome::scalar(exponential().sf(p.scalar("x")))),
        Check {
            function: "statrs::distribution::Exp::inverse_cdf",
            dataset: "expon_ppf",
            expected: &[],
        }
        .run(|p| Outcome::scalar(exponential().inverse_cdf(p.scalar("p")))),
        // `Erlang` calls `Gamma` for `cdf`, `sf` and `inverse_cdf`, so these
        // three run on the gamma datasets' whole-number shapes and fail
        // exactly where `Gamma` does.
        Check {
            function: "statrs::distribution::Erlang::cdf",
            dataset: "gamma_cdf",
            expected: &[Expected {
                finding: "012-statrs-gamma-cdf",
                points: expected::GAMMA_CDF,
            }],
        }
        .run_on(integer_shapes("gamma_cdf"), |p| {
            Outcome::scalar(erlang(p.scalar("a")).cdf(p.scalar("x")))
        }),
        Check {
            function: "statrs::distribution::Erlang::sf",
            dataset: "gamma_sf",
            expected: &[Expected {
                finding: "012-statrs-gamma-cdf",
                points: expected::GAMMA_SF,
            }],
        }
        .run_on(integer_shapes("gamma_sf"), |p| {
            Outcome::scalar(erlang(p.scalar("a")).sf(p.scalar("x")))
        }),
        // Time-limited like `Gamma::inverse_cdf`, which it calls.
        Check {
            function: "statrs::distribution::Erlang::inverse_cdf",
            dataset: "gamma_ppf",
            expected: &[Expected {
                finding: "013-statrs-gamma-inverse-cdf",
                points: expected::GAMMA_PPF_WHOLE,
            }],
        }
        .run_on_with_time_limit(integer_shapes("gamma_ppf"), TIME_LIMIT, |p| {
            Outcome::scalar(erlang(p.scalar("a")).inverse_cdf(p.scalar("p")))
        }),
        Check {
            function: "statrs::distribution::Beta::cdf",
            dataset: "beta_cdf",
            expected: &[Expected {
                finding: "019-statrs-beta-cdf",
                points: expected::BETA_CDF,
            }],
        }
        .run(|p| Outcome::scalar(beta_dist(p.scalar("a"), p.scalar("b")).cdf(p.scalar("x")))),
        Check {
            function: "statrs::distribution::Beta::sf",
            dataset: "beta_sf",
            expected: &[Expected {
                finding: "019-statrs-beta-cdf",
                points: expected::BETA_SF,
            }],
        }
        .run(|p| Outcome::scalar(beta_dist(p.scalar("a"), p.scalar("b")).sf(p.scalar("x")))),
        Check {
            function: "statrs::distribution::FisherSnedecor::cdf",
            dataset: "f_cdf",
            expected: &[Expected {
                finding: "021-statrs-fisher-snedecor-cdf",
                points: expected::F_CDF,
            }],
        }
        .run(|p| f_dist(p, "x", FisherSnedecor::cdf)),
        Check {
            function: "statrs::distribution::FisherSnedecor::sf",
            dataset: "f_sf",
            expected: &[Expected {
                finding: "021-statrs-fisher-snedecor-cdf",
                points: expected::F_SF,
            }],
        }
        .run(|p| f_dist(p, "x", FisherSnedecor::sf)),
        // Near the end: it loops forever at 7 points (finding 020); the
        // abandoned threads only share the CPU with `t_ppf`, whose correct
        // points need microseconds.
        Check {
            function: "statrs::distribution::Beta::inverse_cdf",
            dataset: "betaincinv",
            expected: &[Expected {
                finding: "020-statrs-beta-inverse-cdf",
                points: expected::BETA_PPF,
            }],
        }
        .run_with_time_limit(TIME_LIMIT, |p| {
            Outcome::scalar(beta_dist(p.scalar("a"), p.scalar("b")).inverse_cdf(p.scalar("p")))
        }),
        // It calls `inv_beta_reg` like `Beta::inverse_cdf`, so it can loop
        // forever too.
        Check {
            function: "statrs::distribution::FisherSnedecor::inverse_cdf",
            dataset: "f_ppf",
            expected: &[Expected {
                finding: "022-statrs-fisher-snedecor-inverse-cdf",
                points: expected::F_PPF,
            }],
        }
        .run_with_time_limit(TIME_LIMIT, |p| f_dist(p, "p", FisherSnedecor::inverse_cdf)),
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
