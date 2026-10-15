# Conventions: how every point is judged

The rules here are implemented in `crates/harness/src/compare.rs` (judging) and
`reference/rsv_reference/` (tagging). Changing either is a "Checks changed" item
in the pull request.

## The error metric

```
rel_error(value, reference) = |value - reference| / max(|reference|, 2.2250738585072014e-308)
```

This is ordinary relative error. Below the smallest normal double (about
2.2e-308) the error is measured relative to that number instead: subnormal results
carry fewer significant bits, so demanding full relative accuracy there would be
meaningless. Identical values, including two equal infinities, have error 0.

`reference` is the mpmath value **correctly rounded** to a double. A perfect
implementation has error at most half an ulp (about 1.1e-16).

## Tolerances

Each function has one tolerance, recorded with its reason in the header of its
data file. Every crate is held to the same tolerance for the same function.
Tolerances are set from the mathematics (conditioning and unavoidable rounding),
**before** looking at any crate's results, and are never loosened to make a
crate pass.

## Tags (assigned by rules, never by hand)

The generator assigns tags by fixed rules in code:

| Tag | Rule | mpmath value recorded |
|---|---|---|
| `pole` | the function has a pole or is undefined at this input (per-function rule on the inputs) | `"nan"` (not evaluated) |
| `overflow` | the true value is finite but larger than the largest double | `"inf"` or `"-inf"` |
| `unstable-reference` | mpmath at 50 and at 100 digits rounds to different doubles | the 50-digit value |

## Verdicts

| Verdict | Meaning | Fails CI? |
|---|---|---|
| **Pass** | crate and SciPy both within tolerance of mpmath | no |
| **SciPy wrong** | crate within tolerance, SciPy not. Recorded, not a crate failure | no |
| **Convention** | at a `pole` or `overflow` point, the crate returned an accepted value | no |
| **Skipped** | `unstable-reference`: the arbiter itself is in doubt | no |
| **Crate wrong** | crate outside tolerance, SciPy within | yes, unless listed as an expected failure with a finding |
| **Both wrong** | crate and SciPy both outside tolerance (reviewed by hand) | same as crate wrong |

### Accepted values at tagged points

A tagged point is a *convention difference* only if the crate returns one of:

- **NaN**;
- **an infinity**: either sign at a `pole`, but at an `overflow` point it must have the true value's sign;
- **an error** (e.g. `Err(..)`).

Any **finite** value at a pole or overflow point is **crate wrong**. Example:
statrs `gamma(-6.0)` returning `5.94e12` is a finding, not a convention. A
**panic** is never accepted: it is always crate wrong; nor is a crash (below).

### Calls that never return

A check whose crate function is known to loop forever at some inputs runs with
`Check::run_with_time_limit`: every point runs on its own thread, and a point
with no result when the limit (counted from the start of the check) expires
gets the outcome "no result within N s". Like a panic, that is always crate
wrong. The limit is far above the time a correct call takes (check-statrs uses
5 s for calls that take microseconds), so it only separates an endless loop
from a slow machine. Abandoned threads run until the test process exits, so
such a check comes last in its test.

### Calls that end the process

Some failures cannot be caught inside the process: a stack overflow aborts
it (special's `digamma` recurses once per unit below 8 and does this for
large negative x, finding 017). A check whose crate function can do this runs
with `Check::run_isolated`: its points are evaluated in child processes (the
same test executable, running only that check), one after another. If a
child dies during a point, the point gets the outcome "crash: stack
overflow" (or "crash: process ended (...)" for another cause), always crate
wrong, and a new child continues at the next point. A point with no result
within the time limit after the previous one gets "no result within N s"
and its child is killed. Such a check comes first in its test, so that the
checks before it, which every child skips, do not even load their data
again.

## Expected failures

A crate-wrong point does not fail CI only if the check crate lists it as an
expected failure that cites an existing `findings/NNN-crate-function.md`. Each
check crate keeps one exact list of input points per finding in
`src/expected.rs` (a `harness::Expected` names the finding and its points). CI
fails when:

- a crate-wrong point is not listed;
- a listed point no longer fails (fixed upstream: update the finding);
- a listed point cites a findings file that does not exist.

## Limits of the arbiter

Agreement between mpmath at 50 and at 100 digits catches loss of precision
inside mpmath, but it **does not catch systematic errors**: a wrong formula, a
wrong branch, or a bug in mpmath's algorithm gives the same wrong answer at
every precision. Defences:

- where we implement a formula ourselves on top of mpmath (e.g. the erfcx
  asymptotic series), the generator's tests check it against an independent
  mpmath route where both are valid;
- SciPy is an independent implementation, so a crate-vs-mpmath disagreement
  where SciPy agrees with the crate is reviewed by hand before it becomes a finding;
- findings quote the mathematical identity or published value they rely on.
