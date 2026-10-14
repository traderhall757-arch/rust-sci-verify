# Findings

One file per confirmed bug: `NNN-crate-function.md` (e.g. `001-puruspe-betai.md`).
A finding is confirmed only when the crate disagrees with **mpmath** beyond the
function's tolerance. SciPy alone is not enough.

## Template

```markdown
# NNN: crate::function short description

- **Crate / version:** puruspe 0.4.4
- **Severity:** high | medium | low (wrong answer silently > error/NaN > tiny inaccuracy)
- **Status:** draft | approved by owner | reported (link) | fixed upstream (version)
- **Affected inputs:** the input ranges, from the report

## Reproduction

(minimal Rust program, `cargo new` + one dependency)

## Expected

| Input | Crate | mpmath (true value) | SciPy |
|---|---|---|---|

## Notes

Root cause if known; the identity or source the expected value relies on.
```

Statuses move only forward. "Reported" requires the owner's approval for that
specific report.
