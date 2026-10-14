use harness::{Check, Expected, Outcome, finish};

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
    ];
    finish("puruspe", &results);
}
