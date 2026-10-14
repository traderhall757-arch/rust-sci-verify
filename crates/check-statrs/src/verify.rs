use harness::{Check, Expected, Outcome, finish};

use crate::expected;
use statrs::function::erf;

#[test]
fn statrs_matches_reference() {
    let results = [Check {
        function: "statrs::function::erf::erfc",
        dataset: "erfc",
        expected: &[Expected {
            finding: "004-statrs-erfc",
            points: expected::ERFC,
        }],
    }
    .run(|p| Outcome::scalar(erf::erfc(p.scalar("x"))))];
    finish("statrs", &results);
}
