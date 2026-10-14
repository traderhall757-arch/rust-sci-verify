//! `puruspe::betai` from puruspe 0.4.4 with the fix in
//! `fixes/puruspe/0001-Fix-betai-for-large-a-and-b.patch` applied, so the
//! harness can verify the fix before it is offered upstream. Root cause and
//! evidence: `findings/001-puruspe-betai.md`.
//!
//! **Source and license.** The items below are copied verbatim from puruspe
//! 0.4.4 (<https://github.com/Axect/puruspe>, commit `f3eb94b`) after
//! applying the patch: `betai`, `betacf`, `betaiapprox` and `SWITCH` from
//! `src/beta.rs`; `EPS`, `FPMIN`, `Y` and `W` from `src/lib.rs`. puruspe is
//! MIT OR Apache-2.0, Copyright (c) 2021 Tae Geun Kim; its MIT license text
//! is in `LICENSE-MIT-puruspe`. `ln_gamma` is the pinned crate's own. The
//! algorithms follow Numerical Recipes, 3rd ed., section 6.4.
//! `scripts/check_fixes.py` (CI) checks that the copies equal the patched
//! upstream source, so the code tested here is the code in the patch.
//!
//! **Verified accuracy.** Relative error at most 5.15e-10 on the `betainc`
//! grid and 4.1e-10 on `betainc_large`, against mpmath (reports in
//! `reports/puruspe-fix-0.1.0.md`; CI asserts at most 6e-10). The fix does
//! not touch the remaining error, which comes from the prefactor
//! `exp(ln_gamma(a + b) - ln_gamma(a) - ln_gamma(b) + ...)`; see the finding.
//!
//! The copies keep upstream formatting and idioms, hence the `rustfmt::skip`
//! and the allowed lints.

#![allow(clippy::float_cmp, clippy::needless_range_loop)]

use puruspe::ln_gamma;

#[rustfmt::skip]
const SWITCH: usize = 3000;

#[rustfmt::skip]
const EPS: f64 = f64::EPSILON;
#[rustfmt::skip]
const FPMIN: f64 = f64::MIN_POSITIVE / EPS;
#[rustfmt::skip]
const Y: [f64; 18] = [
    0.0021695375159141994,
    0.011413521097787704,
    0.027972308950302116,
    0.051_727_015_600_492_42,
    0.082_502_225_484_340_94,
    0.12007019910960293,
    0.164_152_833_007_524_7,
    0.21442376986779355,
    0.27051082840644336,
    0.33199876341447887,
    0.39843234186401943,
    0.46931971407375483,
    0.544_136_055_566_579_7,
    0.622_327_452_880_310_8,
    0.703_315_004_655_971_7,
    0.786_499_107_683_134_5,
    0.871_263_896_190_615_2,
    0.956_981_801_526_291_4,
];
#[rustfmt::skip]
const W: [f64; 18] = [
    0.005_565_719_664_244_557,
    0.012_915_947_284_065_42,
    0.020181515297735382,
    0.027298621498568734,
    0.034_213_810_770_299_54,
    0.040_875_750_923_643_26,
    0.047_235_083_490_265_58,
    0.053_244_713_977_759_69,
    0.058_860_144_245_324_8,
    0.064_039_797_355_015_48,
    0.068_745_323_835_736_41,
    0.072_941_885_005_653_09,
    0.076_598_410_645_870_64,
    0.079_687_828_912_071_67,
    0.082_187_266_704_339_7,
    0.084_078_218_979_661_95,
    0.085_346_685_739_338_72,
    0.085_983_275_670_394_82,
];

/// The regularized incomplete beta function
/// I_x(a, b) = B(x; a, b) / B(a, b), for a > 0, b > 0 and 0 <= x <= 1.
///
/// # Panics
///
/// If a <= 0, b <= 0 or x is outside [0, 1] (as upstream).
#[rustfmt::skip]
pub fn betai(a: f64, b: f64, x: f64) -> f64 {
    assert!(a > 0f64 && b > 0f64, "Bad a or b in routine betai");
    assert!((0f64..=1f64).contains(&x), "Bad x in routine betai");
    if x == 0f64 || x == 1f64 {
        return x;
    }
    let switch = SWITCH as f64;
    if a > switch && b > switch {
        return betaiapprox(a, b, x);
    }
    let bt = (ln_gamma(a + b) - ln_gamma(a) - ln_gamma(b) + a * x.ln() + b * (1f64 - x).ln()).exp();
    if x < (a + 1f64) / (a + b + 2f64) {
        bt * betacf(a, b, x) / a
    } else {
        1f64 - bt * betacf(b, a, 1f64 - x) / b
    }
}

#[rustfmt::skip]
fn betacf(a: f64, b: f64, x: f64) -> f64 {
    let qab = a + b;
    let qap = a + 1f64;
    let qam = a - 1f64;
    let mut c = 1f64;
    let mut d = 1f64 - qab * x / qap;
    if d.abs() < FPMIN {
        d = FPMIN;
    }
    d = 1f64 / d;
    let mut h = d;
    for m in 1..10000 {
        let m = m as f64;
        let m2 = 2f64 * m;
        let mut aa = m * (b - m) * x / ((qam + m2) * (a + m2));
        d = 1f64 + aa * d;
        if d.abs() < FPMIN {
            d = FPMIN;
        }
        c = 1f64 + aa / c;
        if c.abs() < FPMIN {
            c = FPMIN;
        }
        d = 1f64 / d;
        h *= d * c;
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
        d = 1f64 + aa * d;
        if d.abs() < FPMIN {
            d = FPMIN;
        }
        c = 1f64 + aa / c;
        if c.abs() < FPMIN {
            c = FPMIN;
        }
        d = 1f64 / d;
        let del = d * c;
        h *= del;
        if (del - 1f64).abs() <= EPS {
            break;
        }
    }
    h
}

#[rustfmt::skip]
fn betaiapprox(a: f64, b: f64, x: f64) -> f64 {
    let a1 = a - 1f64;
    let b1 = b - 1f64;
    let mu = a / (a + b);
    let lnmu = mu.ln();
    let lnmuc = (1f64 - mu).ln();
    let mut t = (a * b / ((a + b).powi(2) * (a + b + 1f64))).sqrt();
    let xu = if x > a / (a + b) {
        if x >= 1f64 {
            return 1f64;
        }
        1f64.min((mu + 10f64 * t).max(x + 5f64 * t))
    } else {
        if x <= 0f64 {
            return 0f64;
        }
        0f64.max((mu - 10f64 * t).min(x - 5f64 * t))
    };
    let mut sum = 0f64;
    for j in 0..18 {
        t = x + (xu - x) * Y[j];
        sum += W[j] * (a1 * (t.ln() - lnmu) + b1 * ((1f64 - t).ln() - lnmuc)).exp();
    }
    let ans = sum
        * (xu - x)
        * (a1 * lnmu - ln_gamma(a) + b1 * lnmuc - ln_gamma(b) + ln_gamma(a + b)).exp();
    if x > mu {
        1f64 - ans
    } else {
        -ans
    }
}
