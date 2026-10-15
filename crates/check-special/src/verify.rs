use std::time::Duration;

use harness::{Check, Expected, Outcome, finish};
use special::{Beta, Gamma};

use crate::expected;

/// How long one point of a time-limited check may take. A correct answer
/// takes microseconds; the limit only has to tell an endless loop apart on a
/// busy CI machine.
const TIME_LIMIT: Duration = Duration::from_secs(5);

/// I_x(a, b), with `ln_beta` precomputed as the documentation asks
/// ("`ln_beta` should be precomputed as `p.ln_beta(q)`").
fn inc_beta(a: f64, b: f64, x: f64) -> f64 {
    x.inc_beta(a, b, a.ln_beta(b))
}

/// P(a, x). special puts the argument first: `x.inc_gamma(a)`. It has no
/// upper function Q, so the `gamma_sf` data is not used.
fn inc_gamma(a: f64, x: f64) -> f64 {
    x.inc_gamma(a)
}

/// The x with I_x(a, b) = p, with `ln_beta` precomputed in the same way.
fn inv_inc_beta(a: f64, b: f64, p: f64) -> f64 {
    p.inv_inc_beta(a, b, a.ln_beta(b))
}

#[test]
fn special_matches_reference() {
    let results = [
        // First, and in child processes: `digamma` recurses once per unit
        // below 8 with no limit, so for large negative x it overflows the
        // stack, which aborts the process (finding 017).
        Check {
            function: "special::Gamma::digamma",
            dataset: "digamma",
            expected: &[Expected {
                finding: "017-special-digamma",
                points: expected::DIGAMMA,
            }],
        }
        .run_isolated("verify::special_matches_reference", TIME_LIMIT, |p| {
            Outcome::scalar(Gamma::digamma(p.scalar("x")))
        }),
        Check {
            function: "special::Gamma::trigamma",
            dataset: "trigamma",
            expected: &[Expected {
                finding: "018-special-trigamma",
                points: expected::TRIGAMMA,
            }],
        }
        .run(|p| Outcome::scalar(Gamma::trigamma(&p.scalar("x")))),
        Check {
            function: "special::Gamma::gamma",
            dataset: "gamma",
            expected: &[],
        }
        .run(|p| Outcome::scalar(Gamma::gamma(p.scalar("x")))),
        Check {
            function: "special::Gamma::ln_gamma",
            dataset: "gammaln",
            expected: &[],
        }
        .run(|p| Outcome::scalar(Gamma::ln_gamma(p.scalar("x")).0)),
        // The second value `ln_gamma` returns: the sign of gamma(x), +1 or -1.
        Check {
            function: "special::Gamma::ln_gamma (sign)",
            dataset: "gammasgn",
            expected: &[],
        }
        .run(|p| Outcome::scalar(f64::from(Gamma::ln_gamma(p.scalar("x")).1))),
        // Time-limited although no point is known to hang: neither the series
        // nor the continued fraction in `inc_gamma` has an iteration limit.
        Check {
            function: "special::Gamma::inc_gamma",
            dataset: "gammainc",
            expected: &[Expected {
                finding: "016-special-inc-gamma",
                points: expected::INC_GAMMA,
            }],
        }
        .run_with_time_limit(TIME_LIMIT, |p| {
            Outcome::scalar(inc_gamma(p.scalar("a"), p.scalar("x")))
        }),
        Check {
            function: "special::Gamma::inc_gamma",
            dataset: "gamma_cdf",
            expected: &[Expected {
                finding: "016-special-inc-gamma",
                points: expected::INC_GAMMA_CDF,
            }],
        }
        .run_with_time_limit(TIME_LIMIT, |p| {
            Outcome::scalar(inc_gamma(p.scalar("a"), p.scalar("x")))
        }),
        Check {
            function: "special::Gamma::inc_gamma",
            dataset: "gammainc_large",
            expected: &[Expected {
                finding: "016-special-inc-gamma",
                points: expected::INC_GAMMA_LARGE,
            }],
        }
        .run_with_time_limit(TIME_LIMIT, |p| {
            Outcome::scalar(inc_gamma(p.scalar("a"), p.scalar("x")))
        }),
        // Time-limited although no point is known to hang: the series in
        // `inc_beta` has no iteration limit.
        Check {
            function: "special::Beta::inc_beta",
            dataset: "betainc",
            expected: &[Expected {
                finding: "014-special-inc-beta",
                points: expected::INC_BETA,
            }],
        }
        .run_with_time_limit(TIME_LIMIT, |p| {
            Outcome::scalar(inc_beta(p.scalar("a"), p.scalar("b"), p.scalar("x")))
        }),
        Check {
            function: "special::Beta::inc_beta",
            dataset: "betainc_large",
            expected: &[Expected {
                finding: "014-special-inc-beta",
                points: expected::INC_BETA_LARGE,
            }],
        }
        .run_with_time_limit(TIME_LIMIT, |p| {
            Outcome::scalar(inc_beta(p.scalar("a"), p.scalar("b"), p.scalar("x")))
        }),
        Check {
            function: "special::Beta::inc_beta",
            dataset: "betainc_switch",
            expected: &[Expected {
                finding: "014-special-inc-beta",
                points: expected::INC_BETA_SWITCH,
            }],
        }
        .run_with_time_limit(TIME_LIMIT, |p| {
            Outcome::scalar(inc_beta(p.scalar("a"), p.scalar("b"), p.scalar("x")))
        }),
        // Last: it loops forever at some points (finding 015), and the
        // abandoned threads keep a CPU busy until the test ends.
        Check {
            function: "special::Beta::inv_inc_beta",
            dataset: "betaincinv",
            expected: &[Expected {
                finding: "015-special-inv-inc-beta",
                points: expected::INV_INC_BETA,
            }],
        }
        .run_with_time_limit(TIME_LIMIT, |p| {
            Outcome::scalar(inv_inc_beta(p.scalar("a"), p.scalar("b"), p.scalar("p")))
        }),
    ];
    finish("special", &results);
}
