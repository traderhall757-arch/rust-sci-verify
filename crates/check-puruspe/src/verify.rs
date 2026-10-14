use harness::{Check, Outcome, finish};

#[test]
fn puruspe_matches_reference() {
    let results = [Check {
        function: "puruspe::erfcx",
        dataset: "erfcx",
        expected: &[],
    }
    .run(|p| Outcome::scalar(puruspe::erfcx(p.scalar("x"))))];
    finish("puruspe", &results);
}
