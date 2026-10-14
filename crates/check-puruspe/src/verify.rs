use harness::{Check, Expected, Outcome, Point, finish};

use crate::expected;

/// `Y_n` through the integer-order function (the data has integer orders only).
fn yn(p: &Point) -> Outcome {
    Outcome::scalar(puruspe::Yn(p.scalar("nu") as u32, p.scalar("x")))
}

/// `Y_nu` through the real-order function, which also returns `J_nu`.
fn ynu(p: &Point) -> Outcome {
    Outcome::scalar(puruspe::Jnu_Ynu(p.scalar("nu"), p.scalar("x")).1)
}

/// `J_n` through the integer-order function.
fn jn(p: &Point) -> Outcome {
    Outcome::scalar(puruspe::Jn(p.scalar("nu") as u32, p.scalar("x")))
}

/// `J_nu` through the real-order function, which also returns `Y_nu`.
fn jnu(p: &Point) -> Outcome {
    Outcome::scalar(puruspe::Jnu_Ynu(p.scalar("nu"), p.scalar("x")).0)
}

/// `I_n` through the integer-order function.
fn in_(p: &Point) -> Outcome {
    Outcome::scalar(puruspe::In(p.scalar("nu") as u32, p.scalar("x")))
}

/// `I_nu` through the real-order function, which also returns `K_nu`.
fn inu(p: &Point) -> Outcome {
    Outcome::scalar(puruspe::Inu_Knu(p.scalar("nu"), p.scalar("x")).0)
}

#[test]
fn puruspe_matches_reference() {
    let results = [
        Check {
            function: "puruspe::erfcx",
            dataset: "erfcx",
            expected: &[],
        }
        .run(|p| Outcome::scalar(puruspe::erfcx(p.scalar("x")))),
        Check {
            function: "puruspe::gammp",
            dataset: "gammainc",
            expected: &[Expected {
                finding: "002-puruspe-gammp",
                points: expected::GAMMP,
            }],
        }
        .run(|p| Outcome::scalar(puruspe::gammp(p.scalar("a"), p.scalar("x")))),
        Check {
            function: "puruspe::betai",
            dataset: "betainc",
            expected: &[Expected {
                finding: "001-puruspe-betai",
                points: expected::BETAI,
            }],
        }
        .run(|p| Outcome::scalar(puruspe::betai(p.scalar("a"), p.scalar("b"), p.scalar("x")))),
        Check {
            function: "puruspe::betai",
            dataset: "betainc_large",
            expected: &[Expected {
                finding: "001-puruspe-betai",
                points: expected::BETAI_LARGE,
            }],
        }
        .run(|p| Outcome::scalar(puruspe::betai(p.scalar("a"), p.scalar("b"), p.scalar("x")))),
        Check {
            function: "puruspe::betai",
            dataset: "betainc_switch",
            expected: &[Expected {
                finding: "001-puruspe-betai",
                points: expected::BETAI_SWITCH,
            }],
        }
        .run(|p| Outcome::scalar(puruspe::betai(p.scalar("a"), p.scalar("b"), p.scalar("x")))),
        Check {
            function: "puruspe::Kn",
            dataset: "kv",
            expected: &[],
        }
        .run(|p| Outcome::scalar(puruspe::Kn(p.scalar("nu") as u32, p.scalar("x")))),
        Check {
            function: "puruspe::Yn",
            dataset: "yv",
            expected: &[],
        }
        .run(yn),
        Check {
            function: "puruspe::Yn",
            dataset: "yv_large",
            expected: &[],
        }
        .run(yn),
        Check {
            function: "puruspe::Jnu_Ynu",
            dataset: "yv",
            expected: &[],
        }
        .run(ynu),
        Check {
            function: "puruspe::Jnu_Ynu",
            dataset: "yv_large",
            expected: &[Expected {
                finding: "005-puruspe-jnu_ynu",
                points: expected::JNU_YNU_LARGE,
            }],
        }
        .run(ynu),
        Check {
            function: "puruspe::Jn",
            dataset: "jv",
            expected: &[Expected {
                finding: "006-puruspe-jn",
                points: expected::JN,
            }],
        }
        .run(jn),
        Check {
            function: "puruspe::Jn",
            dataset: "jv_large",
            expected: &[Expected {
                finding: "006-puruspe-jn",
                points: expected::JN_LARGE,
            }],
        }
        .run(jn),
        Check {
            function: "puruspe::Jnu_Ynu",
            dataset: "jv",
            expected: &[Expected {
                finding: "005-puruspe-jnu_ynu",
                points: expected::JNU_JV,
            }],
        }
        .run(jnu),
        Check {
            function: "puruspe::Jnu_Ynu",
            dataset: "jv_large",
            expected: &[Expected {
                finding: "005-puruspe-jnu_ynu",
                points: expected::JNU_JV_LARGE,
            }],
        }
        .run(jnu),
        Check {
            function: "puruspe::In",
            dataset: "iv",
            expected: &[Expected {
                finding: "007-puruspe-in",
                points: expected::IN,
            }],
        }
        .run(in_),
        Check {
            function: "puruspe::In",
            dataset: "iv_large",
            expected: &[Expected {
                finding: "007-puruspe-in",
                points: expected::IN_LARGE,
            }],
        }
        .run(in_),
        Check {
            function: "puruspe::Inu_Knu",
            dataset: "iv",
            expected: &[Expected {
                finding: "008-puruspe-inu_knu",
                points: expected::INU_IV,
            }],
        }
        .run(inu),
        Check {
            function: "puruspe::Inu_Knu",
            dataset: "iv_large",
            expected: &[Expected {
                finding: "008-puruspe-inu_knu",
                points: expected::INU_IV_LARGE,
            }],
        }
        .run(inu),
        Check {
            function: "puruspe::Inu_Knu",
            dataset: "iv_nonint",
            expected: &[Expected {
                finding: "008-puruspe-inu_knu",
                points: expected::INU_IV_NONINT,
            }],
        }
        .run(inu),
    ];
    finish("puruspe", &results);
}

/// The fix in `fixes/puruspe/`, on the same data as `puruspe::betai`. Its
/// remaining failures are the `ln_gamma` accuracy loss of finding 001, so
/// the expected lists alone would not notice a wrong answer at a listed
/// point; the bounds below (the verified accuracy documented in
/// `puruspe_fix`, per dataset) do.
#[test]
fn puruspe_fix_matches_reference() {
    let betai = |p: &Point| {
        Outcome::scalar(puruspe_fix::betai(
            p.scalar("a"),
            p.scalar("b"),
            p.scalar("x"),
        ))
    };
    let results = [
        Check {
            function: "puruspe_fix::betai",
            dataset: "betainc",
            expected: &[Expected {
                finding: "001-puruspe-betai",
                points: expected::FIX_BETAI,
            }],
        }
        .run(betai),
        Check {
            function: "puruspe_fix::betai",
            dataset: "betainc_large",
            expected: &[Expected {
                finding: "001-puruspe-betai",
                points: expected::FIX_BETAI_LARGE,
            }],
        }
        .run(betai),
        Check {
            function: "puruspe_fix::betai",
            dataset: "betainc_switch",
            expected: &[Expected {
                finding: "001-puruspe-betai",
                points: expected::FIX_BETAI_SWITCH,
            }],
        }
        .run(betai),
    ];
    finish("puruspe-fix", &results);
    for (result, bound) in results.iter().zip([2.5e-10, 5e-10, 3.5e-9]) {
        let worst = result.points.iter().map(|p| p.error).fold(0.0, f64::max);
        assert!(
            worst <= bound,
            "puruspe_fix::betai on {}: worst relative error {worst:e} > {bound:e}",
            result.check.dataset
        );
    }
}
