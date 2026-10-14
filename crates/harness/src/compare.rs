//! Error metric and per-point classification (rules: `docs/conventions.md`).

use crate::check::Outcome;
use crate::data::{Point, Tag, Values};

/// Relative error of `value` against the correctly rounded `reference`.
///
/// `|value - reference| / max(|reference|, f64::MIN_POSITIVE)`: plain
/// relative error, except that below the smallest normal double the error
/// is measured relative to `f64::MIN_POSITIVE`, because subnormal results
/// carry fewer than 53 significant bits. Identical values (including equal
/// infinities and two NaNs) have error 0; any other non-finite pair has
/// infinite error.
#[must_use]
#[expect(
    clippy::float_cmp,
    reason = "exact equality (0.0 == -0.0, inf == inf) is intended"
)]
pub fn rel_error(value: f64, reference: f64) -> f64 {
    if value == reference || (value.is_nan() && reference.is_nan()) {
        return 0.0;
    }
    if !value.is_finite() || !reference.is_finite() {
        return f64::INFINITY;
    }
    (value - reference).abs() / reference.abs().max(f64::MIN_POSITIVE)
}

/// Largest [`rel_error`] over all outputs; infinite if shapes differ.
fn max_error(values: &[Vec<f64>], reference: &Values, outputs: &[String]) -> f64 {
    if values.len() != outputs.len() {
        return f64::INFINITY;
    }
    let mut worst = 0.0_f64;
    for (got, name) in values.iter().zip(outputs) {
        let want = &reference[name];
        if got.len() != want.len() {
            return f64::INFINITY;
        }
        for (g, w) in got.iter().zip(want) {
            worst = worst.max(rel_error(*g, *w));
        }
    }
    worst
}

/// How one point is judged.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Verdict {
    /// Crate and SciPy both agree with mpmath within tolerance.
    Pass,
    /// Crate agrees with mpmath; SciPy does not (SciPy's problem, not a failure).
    ScipyWrong,
    /// At a pole or overflow point the crate returned an accepted value
    /// (NaN, an infinity, or an error). Documented convention, not a bug.
    Convention,
    /// mpmath was not stable across precisions; the point is not judged.
    Skipped,
    /// Crate disagrees with mpmath; SciPy agrees.
    CrateWrong,
    /// Crate and SciPy both disagree with mpmath (reviewed by hand).
    BothWrong,
}

impl Verdict {
    /// True for verdicts that count against the crate.
    #[must_use]
    pub fn is_crate_failure(self) -> bool {
        matches!(self, Self::CrateWrong | Self::BothWrong)
    }
}

/// Accepted crate values at a tagged point: an error, NaN, or an infinity
/// (for overflow, the infinity must have the true value's sign).
fn accepted_at_tag(outcome: &Outcome, point: &Point, tag: Tag) -> bool {
    let values = match outcome {
        Outcome::Error(_) => return true,
        Outcome::Panic(_) => return false,
        Outcome::Value(v) => v,
    };
    let expected = point.mpmath.values().flatten();
    let got: Vec<f64> = values.iter().flatten().copied().collect();
    !got.is_empty()
        && got.len() == expected.clone().count()
        && got.iter().zip(expected).all(|(g, e)| {
            g.is_nan() || (g.is_infinite() && (tag == Tag::Pole || g.signum() == e.signum()))
        })
}

/// Classifies one point. Returns the verdict and the crate's error
/// (infinite for errors, panics and shape mismatches; NaN when not judged
/// numerically, i.e. at tagged points).
#[must_use]
pub fn classify(
    outcome: &Outcome,
    point: &Point,
    outputs: &[String],
    tolerance: f64,
) -> (Verdict, f64) {
    if point.tags.contains(&Tag::UnstableReference) {
        return (Verdict::Skipped, f64::NAN);
    }
    let tag = point
        .tags
        .iter()
        .copied()
        .find(|t| matches!(t, Tag::Pole | Tag::Overflow));
    if let Some(tag) = tag {
        let ok = accepted_at_tag(outcome, point, tag);
        return (
            if ok {
                Verdict::Convention
            } else {
                Verdict::CrateWrong
            },
            f64::NAN,
        );
    }
    let crate_err = match outcome {
        Outcome::Value(v) => max_error(v, &point.mpmath, outputs),
        Outcome::Error(_) | Outcome::Panic(_) => f64::INFINITY,
    };
    let scipy: Vec<Vec<f64>> = outputs.iter().map(|n| point.scipy[n].clone()).collect();
    let scipy_ok = max_error(&scipy, &point.mpmath, outputs) <= tolerance;
    let verdict = match (crate_err <= tolerance, scipy_ok) {
        (true, true) => Verdict::Pass,
        (true, false) => Verdict::ScipyWrong,
        (false, true) => Verdict::CrateWrong,
        (false, false) => Verdict::BothWrong,
    };
    (verdict, crate_err)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data::Dataset;
    use crate::data::tests::sample;

    fn names() -> Vec<String> {
        vec!["y".into()]
    }

    #[test]
    fn rel_error_cases() {
        assert_eq!(rel_error(1.0, 1.0).to_bits(), 0);
        assert_eq!(rel_error(f64::INFINITY, f64::INFINITY).to_bits(), 0);
        assert_eq!(rel_error(f64::NAN, f64::NAN).to_bits(), 0);
        assert_eq!(rel_error(0.0, -0.0).to_bits(), 0);
        assert!(rel_error(f64::NAN, 1.0).is_infinite());
        assert!(rel_error(1.0, f64::INFINITY).is_infinite());
        assert!((rel_error(1.1, 1.0) - 0.1).abs() < 1e-15);
        assert!((rel_error(-2.0, -1.0) - 1.0).abs() < 1e-15);
        // Subnormal reference: measured against the smallest normal double.
        let r = rel_error(0.0, 1e-310);
        assert!((r - 1e-310 / f64::MIN_POSITIVE).abs() < 1e-12);
    }

    fn point(i: usize) -> Point {
        Dataset::parse(&sample()).unwrap().points[i].clone()
    }

    #[test]
    fn untagged_points() {
        let p = point(0); // mpmath 0.1, scipy 0.1
        let n = names();
        assert_eq!(
            classify(&Outcome::scalar(0.1), &p, &n, 1e-13),
            (Verdict::Pass, 0.0)
        );
        let (v, e) = classify(&Outcome::scalar(0.2), &p, &n, 1e-13);
        assert_eq!(v, Verdict::CrateWrong);
        assert!((e - 1.0).abs() < 1e-12);
        assert_eq!(
            classify(&Outcome::Error("e".into()), &p, &n, 1e-13).0,
            Verdict::CrateWrong
        );
        assert_eq!(
            classify(&Outcome::Panic("p".into()), &p, &n, 1e-13).0,
            Verdict::CrateWrong
        );
        let mut bad_scipy = p.clone();
        bad_scipy.scipy.insert("y".into(), vec![f64::INFINITY]);
        assert_eq!(
            classify(&Outcome::scalar(0.1), &bad_scipy, &n, 1e-13).0,
            Verdict::ScipyWrong
        );
        assert_eq!(
            classify(&Outcome::scalar(5.0), &bad_scipy, &n, 1e-13).0,
            Verdict::BothWrong
        );
        // Shape mismatches count as infinitely wrong.
        assert!(
            classify(&Outcome::Value(vec![]), &p, &n, 1e-13)
                .1
                .is_infinite()
        );
        assert!(
            classify(&Outcome::Value(vec![vec![0.1, 0.1]]), &p, &n, 1e-13)
                .1
                .is_infinite()
        );
    }

    #[test]
    fn overflow_accepts_only_correctly_signed_infinity_nan_or_error() {
        let p = point(1); // overflow, mpmath +inf
        let n = names();
        let c = |o: Outcome| classify(&o, &p, &n, 1e-13).0;
        assert_eq!(c(Outcome::scalar(f64::INFINITY)), Verdict::Convention);
        assert_eq!(c(Outcome::scalar(f64::NAN)), Verdict::Convention);
        assert_eq!(c(Outcome::Error("overflow".into())), Verdict::Convention);
        assert_eq!(c(Outcome::scalar(f64::NEG_INFINITY)), Verdict::CrateWrong);
        assert_eq!(c(Outcome::scalar(f64::MAX)), Verdict::CrateWrong);
        assert_eq!(c(Outcome::Panic("boom".into())), Verdict::CrateWrong);
        assert_eq!(c(Outcome::Value(vec![])), Verdict::CrateWrong);
    }

    #[test]
    fn pole_accepts_either_infinity_but_never_a_finite_value() {
        let mut p = point(1);
        p.tags = vec![Tag::Pole];
        p.mpmath.insert("y".into(), vec![f64::NAN]);
        let n = names();
        let c = |o: Outcome| classify(&o, &p, &n, 1e-13).0;
        assert_eq!(c(Outcome::scalar(f64::NEG_INFINITY)), Verdict::Convention);
        assert_eq!(c(Outcome::scalar(f64::INFINITY)), Verdict::Convention);
        // statrs gamma(-6.0) = 5.94e12 must stay a finding.
        assert_eq!(c(Outcome::scalar(5.94e12)), Verdict::CrateWrong);
    }

    #[test]
    fn unstable_reference_is_skipped() {
        let mut p = point(0);
        p.tags = vec![Tag::UnstableReference];
        let (v, e) = classify(&Outcome::scalar(9.0), &p, &names(), 1e-13);
        assert_eq!(v, Verdict::Skipped);
        assert!(e.is_nan());
        assert!(!v.is_crate_failure());
        assert!(Verdict::BothWrong.is_crate_failure());
    }
}
