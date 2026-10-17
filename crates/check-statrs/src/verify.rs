use std::time::Duration;

use harness::{Check, Dataset, Expected, Outcome, finish, workspace_root};

use crate::expected;
use statrs::distribution::{
    Beta, Chi, ChiSquared, ContinuousCDF, Erlang, Exp, FisherSnedecor, Gamma, InverseGamma,
    LogNormal, Normal, Pareto, StudentsT, Weibull,
};
use statrs::function::{beta, erf, gamma};
use statrs::stats_tests::chisquare::chisquare;
use statrs::stats_tests::f_oneway::f_oneway;
use statrs::stats_tests::skewtest::skewtest;
use statrs::stats_tests::ttest_onesample::ttest_onesample;
use statrs::stats_tests::{Alternative, NaNPolicy};

/// The points of the gammaln dataset where Gamma(x) is positive or x is a pole.
/// statrs documents `ln_gamma` as the logarithm of Gamma, not of |Gamma|, so
/// where Gamma(x) < 0 (x < 0 with an odd floor) its NaN is its documented
/// answer and those points are not run. The report counts only these points,
/// and its grid line says so.
fn gamma_not_negative() -> Dataset {
    let mut data =
        Dataset::load(&workspace_root(), "gammaln").unwrap_or_else(|e| panic!("gammaln: {e}"));
    data.points.retain(|p| {
        let x = p.scalar("x");
        x > 0.0 || p.tags.contains(&harness::Tag::Pole) || x.floor().rem_euclid(2.0) == 0.0
    });
    data.header.points = data.points.len();
    data.header.grid = format!(
        "the points of the gammaln grid where Gamma(x) > 0 or x is a pole, {} of its points \
         (statrs's ln_gamma is the logarithm of Gamma, not of |Gamma|). The full grid: {}",
        data.points.len(),
        data.header.grid
    );
    data
}

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

/// The inverse gamma distribution with shape `a` and rate 1 (the reference
/// data's parametrisation).
fn inverse_gamma(a: f64) -> InverseGamma {
    InverseGamma::new(a, 1.0).expect("valid parameters")
}

/// The chi distribution with `k` degrees of freedom (an integer in the data).
fn chi(k: f64) -> Chi {
    assert!(
        k.fract() == 0.0 && k >= 1.0,
        "chi data has a non-integer k: {k}"
    );
    Chi::new(k as u64).expect("valid parameters")
}

/// The Weibull distribution with shape `k` and scale `lambda` (a power of two
/// in the data).
fn weibull(k: f64, lambda: f64) -> Weibull {
    Weibull::new(k, lambda).expect("valid parameters")
}

/// The log-normal distribution with location `mu` and scale `sigma`.
fn lognormal(mu: f64, sigma: f64) -> LogNormal {
    LogNormal::new(mu, sigma).expect("valid parameters")
}

/// The Pareto distribution with scale `xm` (a power of two in the data) and
/// shape `alpha`; statrs takes the scale first.
fn pareto(alpha: f64, xm: f64) -> Pareto {
    Pareto::new(xm, alpha).expect("valid parameters")
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

/// The alternative hypothesis coded as in the reference data: -1 less,
/// 0 two-sided, 1 greater.
fn alternative(code: f64) -> Alternative {
    match code {
        -1.0 => Alternative::Less,
        0.0 => Alternative::TwoSided,
        1.0 => Alternative::Greater,
        _ => panic!("alternative code {code}"),
    }
}

/// Turns a test's `(statistic, p-value)` result into the two outputs.
fn statistic_and_pvalue<E: std::fmt::Debug>(r: Result<(f64, f64), E>) -> Outcome {
    match r {
        Ok((statistic, pvalue)) => Outcome::Value(vec![vec![statistic], vec![pvalue]]),
        Err(e) => Outcome::Error(format!("{e:?}")),
    }
}

/// The counts of a chisquare point (whole numbers in the data).
#[expect(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    reason = "the counts are whole numbers from 0 to about 1e6"
)]
fn counts(p: &harness::Point) -> Vec<usize> {
    p.inputs["f_obs"].iter().map(|&c| c as usize).collect()
}

/// The `ddof` of a chisquare point (0 or 1 in the data).
#[expect(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    reason = "ddof is 0 or 1"
)]
fn ddof(p: &harness::Point) -> usize {
    p.scalar("ddof") as usize
}

/// The points of the chisquare dataset whose expected frequencies are all
/// equal, the ones statrs assumes when `f_exp` is `None`. The report counts
/// only these points, and its grid line says so.
fn equal_expected() -> Dataset {
    let mut data =
        Dataset::load(&workspace_root(), "chisquare").unwrap_or_else(|e| panic!("chisquare: {e}"));
    // Exactly equal: the generator gives every category the same double.
    data.points.retain(|p| {
        let first = p.inputs["f_exp"][0].to_bits();
        p.inputs["f_exp"].iter().all(|e| e.to_bits() == first)
    });
    data.header.points = data.points.len();
    data.header.grid = format!(
        "the points of the chisquare grid with equal expected frequencies, {} of its points, \
         run with f_exp = None. The full grid: {}",
        data.points.len(),
        data.header.grid
    );
    data
}

/// The groups of an f_oneway point: `x` split by `sizes`.
#[expect(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    reason = "the group sizes are whole numbers from 2 to 30"
)]
fn groups(p: &harness::Point) -> Vec<Vec<f64>> {
    let mut rest = p.inputs["x"].as_slice();
    p.inputs["sizes"]
        .iter()
        .map(|&size| {
            let (group, tail) = rest.split_at(size as usize);
            rest = tail;
            group.to_vec()
        })
        .collect()
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
            function: "statrs::function::gamma::ln_gamma",
            dataset: "gammaln",
            expected: &[Expected {
                finding: "030-statrs-ln-gamma",
                points: expected::LN_GAMMA,
            }],
        }
        .run_on(gamma_not_negative(), |p| {
            Outcome::scalar(gamma::ln_gamma(p.scalar("x")))
        }),
        Check {
            function: "statrs::function::gamma::digamma",
            dataset: "digamma",
            expected: &[Expected {
                finding: "031-statrs-digamma",
                points: expected::DIGAMMA,
            }],
        }
        .run(|p| Outcome::scalar(gamma::digamma(p.scalar("x")))),
        Check {
            function: "statrs::function::gamma::gamma_lr",
            dataset: "gammainc",
            expected: &[
                Expected {
                    finding: "032-statrs-gamma-lr-zero",
                    points: expected::GAMMA_LR_ZERO,
                },
                Expected {
                    finding: "012-statrs-gamma-cdf",
                    points: expected::GAMMA_LR,
                },
            ],
        }
        .run(|p| Outcome::scalar(gamma::gamma_lr(p.scalar("a"), p.scalar("x")))),
        Check {
            function: "statrs::function::gamma::gamma_lr",
            dataset: "gammainc_large",
            expected: &[Expected {
                finding: "012-statrs-gamma-cdf",
                points: expected::GAMMA_LR_LARGE,
            }],
        }
        .run(|p| Outcome::scalar(gamma::gamma_lr(p.scalar("a"), p.scalar("x")))),
        Check {
            function: "statrs::function::beta::beta_reg",
            dataset: "betainc",
            expected: &[Expected {
                finding: "019-statrs-beta-cdf",
                points: expected::BETA_REG,
            }],
        }
        .run(|p| Outcome::scalar(beta::beta_reg(p.scalar("a"), p.scalar("b"), p.scalar("x")))),
        Check {
            function: "statrs::function::beta::beta_reg",
            dataset: "betainc_large",
            expected: &[Expected {
                finding: "019-statrs-beta-cdf",
                points: expected::BETA_REG_LARGE,
            }],
        }
        .run(|p| Outcome::scalar(beta::beta_reg(p.scalar("a"), p.scalar("b"), p.scalar("x")))),
        Check {
            function: "statrs::function::beta::beta_reg",
            dataset: "betainc_switch",
            expected: &[Expected {
                finding: "019-statrs-beta-cdf",
                points: expected::BETA_REG_SWITCH,
            }],
        }
        .run(|p| Outcome::scalar(beta::beta_reg(p.scalar("a"), p.scalar("b"), p.scalar("x")))),
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
            function: "statrs::distribution::InverseGamma::cdf",
            dataset: "invgamma_cdf",
            expected: &[
                Expected {
                    finding: "012-statrs-gamma-cdf",
                    points: expected::INVGAMMA_CDF,
                },
                Expected {
                    finding: "024-statrs-inverse-gamma-cdf",
                    points: expected::INVGAMMA_TINY_X,
                },
            ],
        }
        .run(|p| Outcome::scalar(inverse_gamma(p.scalar("a")).cdf(p.scalar("x")))),
        Check {
            function: "statrs::distribution::InverseGamma::sf",
            dataset: "invgamma_sf",
            expected: &[
                Expected {
                    finding: "012-statrs-gamma-cdf",
                    points: expected::INVGAMMA_SF,
                },
                Expected {
                    finding: "024-statrs-inverse-gamma-cdf",
                    points: expected::INVGAMMA_TINY_X,
                },
            ],
        }
        .run(|p| Outcome::scalar(inverse_gamma(p.scalar("a")).sf(p.scalar("x")))),
        // Time-limited like `Gamma::inverse_cdf`: the series and continued
        // fraction inside `cdf` have no iteration limit.
        Check {
            function: "statrs::distribution::InverseGamma::inverse_cdf",
            dataset: "invgamma_ppf",
            expected: &[Expected {
                finding: "025-statrs-inverse-gamma-inverse-cdf",
                points: expected::INVGAMMA_PPF,
            }],
        }
        .run_with_time_limit(TIME_LIMIT, |p| {
            Outcome::scalar(inverse_gamma(p.scalar("a")).inverse_cdf(p.scalar("p")))
        }),
        Check {
            function: "statrs::distribution::Chi::cdf",
            dataset: "chi_cdf",
            expected: &[
                Expected {
                    finding: "012-statrs-gamma-cdf",
                    points: expected::CHI_CDF,
                },
                Expected {
                    finding: "035-statrs-chi-cdf",
                    points: expected::CHI_SQUARE_RANGE,
                },
                Expected {
                    finding: "035-statrs-chi-cdf",
                    points: expected::CHI_CDF_SUBNORMAL_SQUARE,
                },
            ],
        }
        .run(|p| Outcome::scalar(chi(p.scalar("k")).cdf(p.scalar("x")))),
        Check {
            function: "statrs::distribution::Chi::sf",
            dataset: "chi_sf",
            expected: &[
                Expected {
                    finding: "012-statrs-gamma-cdf",
                    points: expected::CHI_SF,
                },
                Expected {
                    finding: "035-statrs-chi-cdf",
                    points: expected::CHI_SQUARE_RANGE,
                },
            ],
        }
        .run(|p| Outcome::scalar(chi(p.scalar("k")).sf(p.scalar("x")))),
        // Time-limited like `InverseGamma::inverse_cdf`: the same search, and
        // `cdf`'s series and continued fraction have no iteration limit.
        Check {
            function: "statrs::distribution::Chi::inverse_cdf",
            dataset: "chi_ppf",
            expected: &[Expected {
                finding: "036-statrs-chi-inverse-cdf",
                points: expected::CHI_PPF,
            }],
        }
        .run_with_time_limit(TIME_LIMIT, |p| {
            Outcome::scalar(chi(p.scalar("k")).inverse_cdf(p.scalar("p")))
        }),
        Check {
            function: "statrs::distribution::Weibull::cdf",
            dataset: "weibull_cdf",
            expected: &[Expected {
                finding: "037-statrs-weibull-cdf",
                points: expected::WEIBULL_CDF,
            }],
        }
        .run(|p| Outcome::scalar(weibull(p.scalar("k"), p.scalar("lambda")).cdf(p.scalar("x")))),
        Check {
            function: "statrs::distribution::Weibull::sf",
            dataset: "weibull_sf",
            expected: &[Expected {
                finding: "037-statrs-weibull-cdf",
                points: expected::WEIBULL_SF,
            }],
        }
        .run(|p| Outcome::scalar(weibull(p.scalar("k"), p.scalar("lambda")).sf(p.scalar("x")))),
        Check {
            function: "statrs::distribution::Weibull::inverse_cdf",
            dataset: "weibull_ppf",
            expected: &[Expected {
                finding: "038-statrs-weibull-inverse-cdf",
                points: expected::WEIBULL_PPF,
            }],
        }
        .run(|p| {
            Outcome::scalar(weibull(p.scalar("k"), p.scalar("lambda")).inverse_cdf(p.scalar("p")))
        }),
        Check {
            function: "statrs::distribution::Pareto::cdf",
            dataset: "pareto_cdf",
            expected: &[Expected {
                finding: "040-statrs-pareto-cdf",
                points: expected::PARETO_CDF,
            }],
        }
        .run(|p| Outcome::scalar(pareto(p.scalar("alpha"), p.scalar("xm")).cdf(p.scalar("x")))),
        Check {
            function: "statrs::distribution::Pareto::sf",
            dataset: "pareto_sf",
            expected: &[Expected {
                finding: "040-statrs-pareto-cdf",
                points: expected::PARETO_SF,
            }],
        }
        .run(|p| Outcome::scalar(pareto(p.scalar("alpha"), p.scalar("xm")).sf(p.scalar("x")))),
        Check {
            function: "statrs::distribution::Pareto::inverse_cdf",
            dataset: "pareto_ppf",
            expected: &[Expected {
                finding: "041-statrs-pareto-inverse-cdf",
                points: expected::PARETO_PPF,
            }],
        }
        .run(|p| {
            Outcome::scalar(pareto(p.scalar("alpha"), p.scalar("xm")).inverse_cdf(p.scalar("p")))
        }),
        Check {
            function: "statrs::distribution::LogNormal::cdf",
            dataset: "lognorm_cdf",
            expected: &[Expected {
                finding: "039-statrs-log-normal",
                points: expected::LOGNORMAL_CDF,
            }],
        }
        .run(|p| Outcome::scalar(lognormal(p.scalar("mu"), p.scalar("sigma")).cdf(p.scalar("x")))),
        Check {
            function: "statrs::distribution::LogNormal::sf",
            dataset: "lognorm_sf",
            expected: &[Expected {
                finding: "039-statrs-log-normal",
                points: expected::LOGNORMAL_SF,
            }],
        }
        .run(|p| Outcome::scalar(lognormal(p.scalar("mu"), p.scalar("sigma")).sf(p.scalar("x")))),
        Check {
            function: "statrs::distribution::LogNormal::inverse_cdf",
            dataset: "lognorm_ppf",
            expected: &[],
        }
        .run(|p| {
            Outcome::scalar(lognormal(p.scalar("mu"), p.scalar("sigma")).inverse_cdf(p.scalar("p")))
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
        Check {
            function: "statrs::stats_tests::ttest_onesample::ttest_onesample",
            dataset: "ttest_1samp",
            expected: &[Expected {
                finding: "026-statrs-ttest-onesample",
                points: expected::TTEST_1SAMP,
            }],
        }
        .run(|p| {
            statistic_and_pvalue(ttest_onesample(
                p.inputs["x"].clone(),
                p.scalar("mu"),
                alternative(p.scalar("alternative")),
                NaNPolicy::Error,
            ))
        }),
        Check {
            function: "statrs::stats_tests::chisquare::chisquare",
            dataset: "chisquare",
            expected: &[Expected {
                finding: "027-statrs-chisquare",
                points: expected::CHISQUARE,
            }],
        }
        .run(|p| {
            statistic_and_pvalue(chisquare(
                &counts(p),
                Some(&p.inputs["f_exp"]),
                Some(ddof(p)),
            ))
        }),
        // The same function with its default, equal expected frequencies.
        Check {
            function: "statrs::stats_tests::chisquare::chisquare",
            dataset: "chisquare",
            expected: &[Expected {
                finding: "027-statrs-chisquare",
                points: expected::CHISQUARE_EQUAL,
            }],
        }
        .run_on(equal_expected(), |p| {
            statistic_and_pvalue(chisquare(&counts(p), None, Some(ddof(p))))
        }),
        Check {
            function: "statrs::stats_tests::f_oneway::f_oneway",
            dataset: "f_oneway",
            expected: &[Expected {
                finding: "028-statrs-f-oneway",
                points: expected::F_ONEWAY,
            }],
        }
        .run(|p| statistic_and_pvalue(f_oneway(groups(p), NaNPolicy::Error))),
        Check {
            function: "statrs::stats_tests::skewtest::skewtest",
            dataset: "skewtest",
            expected: &[
                Expected {
                    finding: "029-statrs-skewtest",
                    points: expected::SKEWTEST,
                },
                Expected {
                    finding: "009-statrs-normal",
                    points: expected::SKEWTEST_ERFC,
                },
            ],
        }
        .run(|p| {
            statistic_and_pvalue(skewtest(
                p.inputs["x"].clone(),
                alternative(p.scalar("alternative")),
                NaNPolicy::Error,
            ))
        }),
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
