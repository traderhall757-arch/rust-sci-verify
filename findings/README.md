# Findings

One file per confirmed bug: `NNN-crate-function.md` (e.g. `001-puruspe-betai.md`).
A finding is confirmed only when the crate disagrees with **mpmath** beyond the
function's tolerance. SciPy alone is not enough.

## Template

```markdown
# NNN: crate::function short description

- **Crate / version:** puruspe 0.4.4
- **Severity:** high | medium | low (wrong answer silently > error/NaN > tiny inaccuracy)
- **Status:** draft | approved (date) | reported (link) | fixed upstream (version)
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

## Approval

Claude approves findings itself (owner, 2026-10-14), but only after every
step below is done and recorded in the finding. Until then the status stays
"draft". Findings 001 to 005 were approved by the owner before this rule.

1. **Reproduced** in a fresh `cargo new` project with the crate at the pinned
   version; the program in "Reproduction" is the one that was run.
2. **Arbiter checked**: mpmath agrees at 50 and 100 digits (the generator
   enforces this), and its route is tested against an independent one
   (identity, recurrence, expansion or finite sum) in
   `reference/tests/test_functions.py`.
3. **Second opinion** on the expected values in the "Expected" table: SciPy
   within tolerance, or, where SciPy is also wrong, a second mpmath route.
4. **Off-grid confirmation**: a random check or scan of the affected range,
   with method and seed in the finding, shows the bug is not an artifact of
   the grid and gives the affected range ("Affected inputs").
5. **Root cause** identified from the source, or the leads written down;
   severity set by the template's rule.

The approval line names the date and says "approved (Claude, checklist
complete)". The owner can reopen any finding. Approval never covers
contacting another project (rule 8 in `CLAUDE.md`).
