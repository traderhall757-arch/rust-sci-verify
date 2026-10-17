"""Tests for select_datasets.py: which datasets a pull request regenerates."""

import subprocess
from pathlib import Path

import pytest
import select_datasets as sd

REPO = Path(__file__).resolve().parent.parent.parent
REAL = (REPO / sd.FUNCTIONS).read_text()

BASE = '''"""Docstring."""

import math

from .spec import Spec

_X = (1.0, 2.0)


def _shared() -> list[float]:
    return [math.sqrt(x) for x in _X]


def _a_grid():
    return [(x,) for x in _shared()]


def _b_grid():
    # b only
    return [(x + 1,) for x in _shared()]


def _c_grid():
    return [(3.0,)]


A = Spec(name="a", grid=_a_grid, tolerance=1e-13)
B = Spec(name="b", grid=_b_grid, tolerance=1e-13)
C = Spec(name="c", grid=_c_grid, tolerance=1e-13)

ALL = (A, B, C)
'''


def selection(head: str, base: str | None = BASE, changed=(sd.FUNCTIONS,), package=None):
    return sd.select(list(changed), base, head, package or {})


def regenerated(s: sd.Selection) -> str | list[str]:
    return "all" if s.everything else s.names


# analyse


def test_analyse_finds_datasets_in_order():
    r = sd.analyse(BASE)
    assert r.order == ["a", "b", "c"]
    assert len(set(r.datasets.values())) == 3


def test_real_functions_py_is_analysable_and_matches_the_data_files():
    r = sd.analyse(REAL)
    files = sorted(p.stem for p in (REPO / sd.DATA).glob("*.jsonl"))
    assert sorted(r.order) == files


@pytest.mark.parametrize(
    "source, reason",
    [
        ("x = globals()['A']", "globals"),
        ("def f():\n    return eval('A')", "eval"),
        ("def f(g):\n    return g.__globals__", "__globals__"),
        ("import sys\nALL = ()", "imports sys"),
        ("from . import functions\nALL = ()", "imports ."),
        ("def f(): pass\ndef f(): pass\nALL = ()", "twice"),
        ("A = Spec(name='a')\nALL = [A, A]", "listed twice"),
        ("A = Spec(name='a' + 'b')\nALL = (A,)", "no literal name"),
        ("A = Spec(name='a')\nALL = (A, make())", "not a name bound to Spec"),
        ("A = Spec(name='a')\nALL = tuple([A])", "not a literal"),
        ("A = 1", "no ALL"),
    ],
)
def test_analyse_rejects_what_it_cannot_follow(source, reason):
    with pytest.raises(sd.Unsupported, match=reason):
        sd.analyse(source)


# select


def test_comment_or_formatting_change_regenerates_nothing():
    head = BASE.replace("# b only", "# only b").replace("(x + 1,)", "( x+1 , )")
    s = selection(head)
    assert regenerated(s) == []
    assert all(line.startswith("skip ") for line in s.report)
    assert len(s.report) == 3  # every dataset is reported, skipped ones too


def test_one_spec_changed_regenerates_only_that_dataset():
    head = BASE.replace(
        'name="b", grid=_b_grid, tolerance=1e-13', 'name="b", grid=_b_grid, tolerance=1e-12'
    )
    s = selection(head)
    assert regenerated(s) == ["b"]
    assert "regenerate b: its Spec or code it uses changed" in s.report


def test_helper_change_regenerates_every_dataset_using_it():
    assert regenerated(selection(BASE.replace("x + 1", "x + 2"))) == ["b"]
    assert regenerated(selection(BASE.replace("math.sqrt(x)", "math.cbrt(x)"))) == ["a", "b"]
    assert regenerated(selection(BASE.replace("(1.0, 2.0)", "(1.0, 3.0)"))) == ["a", "b"]


def test_new_dataset_and_changed_data_file():
    head = BASE.replace("ALL = (A, B, C)", 'D = Spec(name="d", grid=_c_grid)\n\nALL = (A, B, C, D)')
    assert regenerated(selection(head)) == ["d"]
    s = selection(BASE, changed=[sd.DATA + "c.jsonl"])
    assert regenerated(s) == ["c"]
    assert "regenerate c: its data file changed" in s.report


def test_reordering_all_regenerates_nothing():
    # Each dataset runs in its own process, so order cannot matter.
    assert regenerated(selection(BASE.replace("(A, B, C)", "(C, A, B)"))) == []


@pytest.mark.parametrize(
    "head",
    [
        BASE.replace("import math", "import math\nimport itertools"),  # import
        BASE.replace("_X = (1.0, 2.0)", "_X = (1.0, 2.0)\nprint(_X)"),  # module statement
        BASE.replace("_X = (1.0, 2.0)", "_X = (1.0, 2.0)\n_Y = sorted(_X)"),  # call at import
        BASE.replace("def _c_grid():", "@decorate\ndef _c_grid():"),  # decorator call
        BASE.replace("_X = (1.0, 2.0)", "_X = (1.0, 2.0)\nx = y = 1"),  # binds two names
        BASE.replace("return [(3.0,)]", "return eval('[(3.0,)]')"),  # dynamic lookup
    ],
)
def test_shared_or_unanalysable_change_regenerates_all(head):
    assert regenerated(selection(head)) == "all"


def test_code_run_at_import_is_shared_with_what_it_reaches():
    # _Y is computed at import by calling _shared, so a change to _shared (or
    # anything it reaches) can affect every dataset.
    base = BASE.replace("ALL = (A, B, C)", "_Y = tuple(_shared())\n\nALL = (A, B, C)")
    head = base.replace("math.sqrt(x)", "math.cbrt(x)")
    assert regenerated(selection(head, base)) == "all"
    # A lambda body does not run at import; only its defaults do.
    base = BASE.replace("ALL = (A, B, C)", "_L = lambda t=1.0: _shared()\n\nALL = (A, B, C)")
    assert regenerated(selection(base.replace("x + 1", "x + 2"), base)) == ["b"]
    lam = base.replace("t=1.0", "t=float(1)")
    assert regenerated(selection(lam, base)) == "all"


@pytest.mark.parametrize(
    "path",
    [
        "reference/rsv_reference/generate.py",
        "reference/rsv_reference/exact.py",
        "reference/rsv_reference/spec.py",
        "reference/rsv_reference/__main__.py",
        "reference/pyproject.toml",
        "reference/uv.lock",
        "reference/.python-version",
        ".github/workflows/ci.yml",
        "scripts/select_datasets.py",
    ],
)
def test_generator_environment_or_workflow_change_regenerates_all(path):
    s = selection(BASE, changed=[path])
    assert regenerated(s) == "all"
    assert s.report == [f"regenerate all: {path} changed"]


def test_unrelated_changes_regenerate_nothing():
    changed = ["reference/tests/test_functions.py", "crates/harness/src/lib.rs", "docs/status.md"]
    assert regenerated(selection(BASE, changed=changed)) == []


@pytest.mark.parametrize(
    "spec_imports",
    [
        "from . import functions",
        "from .functions import ALL",
        "from rsv_reference import functions",
        "from rsv_reference.functions import ALL",
        "import rsv_reference.functions",
        "import rsv_reference",
        "from .exact import x",  # exact.py imports functions
    ],
)
def test_module_imported_by_functions_importing_it_back_regenerates_all(spec_imports):
    package = {
        sd.PACKAGE + "spec.py": f"import math\n{spec_imports}\n",
        sd.PACKAGE + "exact.py": "from .functions import ALL\n",
    }
    s = selection(BASE, package=package)
    assert regenerated(s) == "all"
    assert "imports functions.py, which imports it" in s.report[0]


def test_driver_may_import_functions():
    # functions.py does not import the driver, so the driver cannot be used
    # to look definitions up from inside it.
    package = {
        sd.PACKAGE + "spec.py": "import math\nfrom .exact import x\n",
        sd.PACKAGE + "exact.py": "import math\n",
        sd.PACKAGE + "cli.py": "from .functions import ALL\n",
        sd.PACKAGE + "__main__.py": "from .cli import main\n",
    }
    assert regenerated(selection(BASE, package=package)) == []


def test_real_package_does_not_import_functions_back():
    package = {
        f"{sd.PACKAGE}{p.name}": p.read_text()
        for p in (REPO / sd.PACKAGE).glob("*.py")
        if p.name != "functions.py"
    }
    assert "cli" in {p.split("/")[-1][:-3] for p in package}
    assert sd._imports_functions_back(REAL, package) is None


def test_new_or_unanalysable_base_regenerates_all():
    assert regenerated(selection(BASE, base=None)) == "all"
    assert regenerated(selection(BASE, base="ALL = make()")) == "all"


def test_real_file_one_spec_or_one_helper():
    erfcx = REAL.replace('name="erfcx",', 'name="erfcx",  # touched', 1)
    assert erfcx != REAL
    assert regenerated(selection(erfcx, REAL)) == []  # a comment only
    head = REAL.replace("costs up to x^2 * 2^-53", "costs at most x^2 * 2^-53", 1)
    assert head != REAL
    assert regenerated(selection(head, REAL)) == ["erfcx"]
    # erfc, the normal and log-normal distributions and the Levy cdf and quantile
    # are computed from _mp_erfcx.
    head = REAL.replace("def _mp_erfcx(x: mpmath.mpf)", "def _mp_erfcx(x: mpmath.mpf, _unused=0)")
    assert regenerated(selection(head, REAL)) == [
        "erfcx",
        "erfc",
        "norm_cdf",
        "norm_sf",
        "norm_ppf",
        "lognorm_cdf",
        "lognorm_sf",
        "lognorm_ppf",
        "levy_cdf",
        "levy_ppf",
    ]
    head = REAL.replace("range(-300, 301, 10)", "range(-300, 301, 20)")
    used = regenerated(selection(head, REAL))
    assert "erfcx" in used and "jv" not in used


# main, on a small git repository


def git(root: Path, *args: str) -> None:
    subprocess.run(["git", "-C", str(root), *args], check=True, capture_output=True)


def commit(root: Path, files: dict[str, str | None]) -> None:
    for rel, text in files.items():
        p = root / rel
        if text is None:
            p.unlink()
            continue
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(text)
    git(root, "add", "-A")
    git(root, "-c", "user.name=t", "-c", "user.email=t@t", "commit", "-qm", "c")


@pytest.fixture
def repo(tmp_path):
    git(tmp_path, "init", "-q")
    files = {sd.FUNCTIONS: BASE, sd.PACKAGE + "generate.py": "import math\n"}
    files |= {f"{sd.DATA}{n}.jsonl": "{}\n" for n in "abc"}
    commit(tmp_path, files)
    return tmp_path


def run(repo: Path, capsys, registry: str = "a b c") -> tuple[int, str, str]:
    code = sd.main(["--base", "HEAD~1", "--repo", str(repo), "--registry", registry])
    out = capsys.readouterr()
    return code, out.out.strip(), out.err


def test_main_prints_selection_and_reports_every_dataset(repo, capsys):
    commit(repo, {sd.FUNCTIONS: BASE.replace("x + 1", "x + 2")})
    code, out, err = run(repo, capsys)
    assert (code, out) == (0, "b")
    assert "regenerate b:" in err and "skip a:" in err and "skip c:" in err


def test_main_all_and_none(repo, capsys):
    commit(repo, {sd.PACKAGE + "generate.py": "import math  # changed\n"})
    assert run(repo, capsys)[:2] == (0, "all")
    commit(repo, {"README.md": "hi\n"})
    assert run(repo, capsys)[:2] == (0, "none")


def test_main_unanalysable_head_regenerates_all(repo, capsys):
    commit(repo, {sd.FUNCTIONS: BASE.replace("ALL = (A, B, C)", "ALL = tuple([A, B, C])")})
    code, out, err = run(repo, capsys)
    assert (code, out) == (0, "all")
    assert "not a literal" in err


def test_main_fails_if_source_and_imported_registry_differ(repo, capsys):
    code, _, err = run(repo, capsys, registry="a b c d")
    assert code == 1 and "datasets without a file ['d']" in err
    commit(repo, {f"{sd.DATA}d.jsonl": "{}\n"})
    code, _, err = run(repo, capsys, registry="a b c d")
    assert code == 1 and "differs from the imported registry" in err
    assert run(repo, capsys, registry="a b c")[0] == 1  # d.jsonl has no dataset


def test_main_fails_without_functions_py(repo, capsys):
    commit(repo, {sd.FUNCTIONS: None})
    code, _, err = run(repo, capsys)
    assert code == 1 and "does not exist" in err
