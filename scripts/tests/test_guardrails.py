"""Tests for the guardrail scripts, on small fake workspaces."""

import subprocess
from pathlib import Path

import check_api_tests
import check_orphans
import check_pins


def write(root: Path, rel: str, text: str) -> None:
    p = root / rel
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(text)


def crate(root: Path, name: str, role: str | None, src: str, deps: str = "") -> None:
    meta = f'\n[package.metadata.rsv]\nrole = "{role}"\n' if role else ""
    write(root, f"crates/{name}/Cargo.toml", f'[package]\nname = "{name}"\n{meta}{deps}')
    write(root, f"crates/{name}/src/lib.rs", src)


def workspace(root: Path, *members: str) -> None:
    listed = ", ".join(f'"crates/{m}"' for m in members)
    write(root, "Cargo.toml", f"[workspace]\nmembers = [{listed}]\n")


# check_orphans


def test_dep_info_path_for_libraries_and_test_binaries():
    deps = Path("/t/debug/deps")
    assert check_orphans.dep_info_path(deps / "libharness-ab12.rlib") == deps / "harness-ab12.d"
    assert check_orphans.dep_info_path(deps / "libharness-ab12.rmeta") == deps / "harness-ab12.d"
    assert check_orphans.dep_info_path(deps / "harness-ab12") == deps / "harness-ab12.d"
    assert check_orphans.dep_info_path(deps / "libfoo-1.so") == deps / "foo-1.d"


def test_parse_dep_info_and_orphans(tmp_path):
    write(tmp_path, "crates/a/src/lib.rs", "")
    write(tmp_path, "crates/a/src/used.rs", "")
    write(tmp_path, "crates/a/src/stray.rs", "")
    write(tmp_path, "crates/a/target/x.rs", "")
    text = (
        f"{tmp_path}/out.d: crates/a/src/lib.rs crates/a/src/used.rs\n\n"
        "crates/a/src/lib.rs:\n"
        f"{tmp_path}/crates/a/src/used.rs:\n"
        "Cargo.toml:\n"
    )
    compiled = check_orphans.parse_dep_info(text, tmp_path)
    assert compiled == {
        (tmp_path / "crates/a/src/lib.rs").resolve(),
        (tmp_path / "crates/a/src/used.rs").resolve(),
    }
    assert check_orphans.orphans(tmp_path, compiled) == [
        (tmp_path / "crates/a/src/stray.rs").resolve()
    ]


# check_api_tests


def test_api_rules(tmp_path):
    workspace(tmp_path, "impl", "check-impl", "tool", "norole")
    crate(
        tmp_path,
        "impl",
        "implementation",
        "pub fn fresnel() {}\npub(crate) fn hidden() {}\npub fn sici() {}\n",
    )
    crate(tmp_path, "check-impl", "check", 'const C: &str = "x"; // function: "impl::fresnel"\n')
    crate(tmp_path, "tool", "harness", "pub fn anything() {}\n")
    crate(tmp_path, "norole", None, "")
    found = check_api_tests.problems(tmp_path)
    # fresnel is checked; pub(crate) hidden is internal; sici is not checked.
    assert found == [
        "impl: pub fn sici has no reference check",
        "norole: [package.metadata.rsv] role must be one of ['check', 'harness', 'implementation']",
    ]


def test_check_crates_must_not_export(tmp_path):
    workspace(tmp_path, "check-x")
    crate(tmp_path, "check-x", "check", "pub fn helper() {}\n")
    assert check_api_tests.problems(tmp_path) == [
        "check-x: check crates must not define pub fn (helper)"
    ]


# check_pins


def pinned_repo(
    root: Path, toolchain="1.97.0", py_dep="scipy==1.18.1", target='"=0.4.4"', tool="x@1.0"
):
    write(root, "rust-toolchain.toml", f'[toolchain]\nchannel = "{toolchain}"\n')
    write(root, "Cargo.lock", "")
    write(root, "reference/uv.lock", "")
    write(
        root,
        "reference/pyproject.toml",
        f'[project]\ndependencies = ["{py_dep}"]\n[dependency-groups]\ndev = ["pytest==9.1.1"]\n',
    )
    workspace(root, "check-p")
    crate(
        root,
        "check-p",
        "check",
        "",
        f"[dev-dependencies]\nh = {{ workspace = true }}\np = {target}\n",
    )
    write(root, ".github/workflows/ci.yml", f"steps:\n  - with:\n      tool: {tool}\n")
    subprocess.run(["git", "init", "-q"], cwd=root, check=True)
    subprocess.run(["git", "add", "-A"], cwd=root, check=True)


def test_pins_clean(tmp_path):
    pinned_repo(tmp_path)
    assert check_pins.problems(tmp_path) == []


def test_pins_violations(tmp_path):
    pinned_repo(tmp_path, toolchain="stable", py_dep="scipy>=1", target='"0.4"', tool="x@1,y")
    subprocess.run(["git", "rm", "-q", "--cached", "Cargo.lock"], cwd=tmp_path, check=True)
    assert check_pins.problems(tmp_path) == [
        "rust-toolchain.toml: channel 'stable' is not an exact x.y.z release",
        "Cargo.lock: not committed",
        "reference/pyproject.toml: 'scipy>=1' is not pinned with ==",
        "crates/check-p: p = '0.4' is not pinned with =x.y.z",
        "ci.yml: install-action tool 'y' has no @version",
    ]
