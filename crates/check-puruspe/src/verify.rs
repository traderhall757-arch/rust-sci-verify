use harness::{Check, Expected, Outcome, Point, finish};

use crate::expected;

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
            function: "puruspe::Kn",
            dataset: "kv",
            expected: &[],
        }
        .run(|p| Outcome::scalar(puruspe::Kn(p.scalar("nu") as u32, p.scalar("x")))),
    ];
    finish("puruspe", &results);
}

/// The fix in `fixes/puruspe/`, on the same data as `puruspe::betai`. Its
/// remaining failures are the `ln_gamma` accuracy loss of finding 001, so
/// the expected lists alone would not notice a wrong answer at a listed
/// point; the bound below (the verified accuracy documented in
/// `puruspe_fix`) does.
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
    ];
    finish("puruspe-fix", &results);
    let worst = results
        .iter()
        .flat_map(|r| &r.points)
        .map(|p| p.error)
        .fold(0.0, f64::max);
    assert!(
        worst <= 6e-10,
        "puruspe_fix::betai: worst relative error {worst:e}"
    );
}
