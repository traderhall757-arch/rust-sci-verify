"""Tests for scripts/check_fixes.py, on a small fake upstream crate and patch."""

from pathlib import Path

import check_fixes

UPSTREAM = """\
const K: f64 = 2.0;
const T: [f64; 2] = [
    1.0,
    2.0,
];

/// Doc.
pub fn f(x: f64) -> f64 {
    x * K
}

fn g(x: f64) -> f64 {
    x + 1.0
}
"""

PATCH = """\
--- a/src/lib.rs
+++ b/src/lib.rs
@@ -8,4 +8,4 @@
 pub fn f(x: f64) -> f64 {
-    x * K
+    x * K + 1.0
 }

"""

FIXED_F = "pub fn f(x: f64) -> f64 {\n    x * K + 1.0\n}"


def write(root: Path, rel: str, text: str) -> None:
    p = root / rel
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(text)


def setup(root: Path, copy: str, patch: str = PATCH) -> dict[str, Path]:
    """A workspace with one fix crate, and the fake upstream crate it patches."""
    write(root, "upstream/Cargo.toml", '[package]\nname = "up"\n')
    write(root, "upstream/src/lib.rs", UPSTREAM)
    write(root, "ws/Cargo.toml", '[workspace]\nmembers = ["crates/up-fix", "crates/other"]\n')
    write(root, "ws/crates/other/Cargo.toml", '[package]\nname = "other"\n')
    write(
        root,
        "ws/crates/up-fix/Cargo.toml",
        '[package]\nname = "up-fix"\n\n[package.metadata.rsv.fix]\nupstream = "up"\n'
        'patch = "fixes/up/0001.patch"\nitems = { "src/lib.rs" = ["f", "K", "T"] }\n',
    )
    write(root, "ws/crates/up-fix/src/lib.rs", copy)
    write(root, "ws/fixes/up/0001.patch", patch)
    return {"up": root / "upstream"}


GOOD_COPY = (
    "use x;\n\nconst K: f64 = 2.0;\n#[rustfmt::skip]\nconst T: [f64; 2] = [\n    1.0,\n"
    "    2.0,\n];\n\n/// Our doc.\n#[rustfmt::skip]\n" + FIXED_F + "\n"
)


def test_extract_item_fn_and_consts():
    assert check_fixes.extract_item(UPSTREAM, "K") == "const K: f64 = 2.0;"
    assert (
        check_fixes.extract_item(UPSTREAM, "T") == "const T: [f64; 2] = [\n    1.0,\n    2.0,\n];"
    )
    assert check_fixes.extract_item(UPSTREAM, "f") == "pub fn f(x: f64) -> f64 {\n    x * K\n}"
    assert check_fixes.extract_item(UPSTREAM, "g") == "fn g(x: f64) -> f64 {\n    x + 1.0\n}"
    assert check_fixes.extract_item("pub(crate) fn h() {\n}\n", "h") == "pub(crate) fn h() {\n}"


def test_extract_item_absent_or_unterminated():
    assert check_fixes.extract_item(UPSTREAM, "missing") is None
    # A name that only prefixes another item does not match it.
    assert check_fixes.extract_item(UPSTREAM, "fn") is None
    assert check_fixes.extract_item("fn f() {\n    1\n", "f") is None


def test_matching_copy_passes(tmp_path):
    ups = setup(tmp_path, GOOD_COPY)
    assert check_fixes.fix_crates(tmp_path / "ws")[0][0] == "up-fix"
    assert check_fixes.problems(tmp_path / "ws", ups, run_tests=False) == []


def test_copy_that_differs_from_the_patch_fails(tmp_path):
    ups = setup(tmp_path, GOOD_COPY.replace("x * K + 1.0", "x * K + 2.0"))
    assert check_fixes.problems(tmp_path / "ws", ups, run_tests=False) == [
        "up-fix: f differs from patched src/lib.rs"
    ]


def test_copy_of_the_unpatched_code_fails(tmp_path):
    ups = setup(tmp_path, GOOD_COPY.replace("x * K + 1.0", "x * K"))
    assert check_fixes.problems(tmp_path / "ws", ups, run_tests=False) == [
        "up-fix: f differs from patched src/lib.rs"
    ]


def test_missing_items_are_reported(tmp_path):
    ups = setup(tmp_path, GOOD_COPY.replace("const K: f64 = 2.0;\n", ""))
    fix_toml = tmp_path / "ws/crates/up-fix/Cargo.toml"
    fix_toml.write_text(fix_toml.read_text().replace('"T"]', '"T", "nope"]'))
    assert check_fixes.problems(tmp_path / "ws", ups, run_tests=False) == [
        "up-fix: K not found in the fix crate",
        "up-fix: nope not found in patched src/lib.rs",
    ]


def test_patch_that_does_not_apply_fails(tmp_path):
    ups = setup(tmp_path, GOOD_COPY, patch=PATCH.replace("x * K\n", "x * Q\n"))
    found = check_fixes.problems(tmp_path / "ws", ups, run_tests=False)
    assert len(found) == 1
    assert found[0].startswith("up-fix: fixes/up/0001.patch does not apply")


def test_unknown_upstream_fails(tmp_path):
    setup(tmp_path, GOOD_COPY)
    assert check_fixes.problems(tmp_path / "ws", {}, run_tests=False) == [
        "up-fix: upstream crate 'up' is not a locked dependency"
    ]


def test_patched_upstream_tests_run(tmp_path):
    # A real (tiny) crate: its test fails after the patch, so --run-tests reports it.
    ups = setup(tmp_path, GOOD_COPY)
    write(
        tmp_path,
        "upstream/Cargo.toml",
        '[package]\nname = "up"\nversion = "0.1.0"\nedition = "2021"\n',
    )
    write(
        tmp_path,
        "upstream/src/lib.rs",
        UPSTREAM + "\n#[test]\nfn t() {\n    assert_eq!(f(1.0), 2.0);\n}\n",
    )
    write(
        tmp_path,
        "upstream/Cargo.lock",
        'version = 4\n\n[[package]]\nname = "up"\nversion = "0.1.0"\n',
    )
    found = check_fixes.problems(tmp_path / "ws", ups, run_tests=True)
    assert len(found) == 1
    assert found[0].startswith("up-fix: patched upstream tests fail:")
