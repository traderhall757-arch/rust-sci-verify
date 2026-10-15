"""Choose the reference datasets a pull request must regenerate.

Usage: select_datasets.py --base REV [--head REV] --registry "NAME ..."

Prints `all`, `none` or the dataset names to stdout, and one line per dataset
saying whether it is regenerated and why to stderr. The CI job
`reference-fresh` passes the result to `python -m rsv_reference`; pushes to
`main` and the weekly run always regenerate everything.

Why skipping a dataset is safe. `python -m rsv_reference` generates every
dataset in its own fresh process, so a dataset's file depends only on

1. the Python environment and the generator (`reference/` outside
   `functions.py`, `data/` and `tests/`; the CI workflow; this script):
   any change regenerates everything;
2. the code of `functions.py` that runs in every process at import (imports,
   statements that bind no single name, and assignments that call anything
   but `Spec`): any change regenerates everything;
3. the definitions its `Spec` reaches through names, transitively.

A dataset is skipped only when (3) is unchanged, compared as syntax trees
(comments and formatting do not count), and its data file is unchanged too.
Anything the analysis cannot follow (dynamic lookups such as `globals()`,
names bound twice, unknown imports, a module imported by `functions.py` that
imports it back) regenerates everything instead of guessing.
"""

import argparse
import ast
import hashlib
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path

PACKAGE = "reference/rsv_reference/"
FUNCTIONS = PACKAGE + "functions.py"
DATA = "reference/data/"
TESTS = "reference/tests/"
# Changes to these can affect every dataset (besides reference/ itself).
EVERYTHING = (".github/workflows/ci.yml", "scripts/select_datasets.py")
# The only imports functions.py may have; anything else could reach code or
# state the analysis does not see.
IMPORTS = {"itertools", "math", "mpmath", "numpy", "scipy.special", "scipy.stats", ".spec"}
# Names and attributes that look up module contents by string.
DYNAMIC_NAMES = {"eval", "exec", "compile", "globals", "locals", "vars", "__import__", "__name__"}
DYNAMIC_ATTRS = {"__dict__", "__globals__", "__builtins__", "modules", "import_module"}
# The only call allowed in a definition that runs at import without making it
# shared: Spec is a frozen dataclass from spec.py.
IMPORT_TIME_CALLS = {"Spec"}


class Unsupported(Exception):
    """functions.py uses something the analysis cannot follow."""


@dataclass(frozen=True)
class Registry:
    """What functions.py defines, reduced to fingerprints."""

    shared: str  # code every dataset depends on
    datasets: dict[str, str]  # dataset name -> code only its Spec reaches
    order: list[str]  # dataset names in ALL order


@dataclass(frozen=True)
class Selection:
    everything: bool
    report: list[str]  # one line per dataset, or the reason for regenerating all
    names: list[str] = field(default_factory=list)  # if not everything


def _fingerprint(nodes: list[ast.AST]) -> str:
    h = hashlib.sha256()
    for node in nodes:
        h.update(ast.dump(node).encode())
        h.update(b"\0")
    return h.hexdigest()


def _import_time_parts(node: ast.stmt) -> list[ast.AST]:
    """The parts of a top-level definition evaluated when the module is imported."""
    if isinstance(node, ast.FunctionDef):
        a = node.args
        args = [*a.posonlyargs, *a.args, *a.kwonlyargs, a.vararg, a.kwarg]
        notes = [x.annotation for x in args if x is not None and x.annotation is not None]
        defaults = [*a.defaults, *(d for d in a.kw_defaults if d is not None)]
        return [*node.decorator_list, *defaults, *notes, *([node.returns] if node.returns else [])]
    if isinstance(node, ast.AnnAssign):
        return [node.annotation, node.value]
    return [node.value]


def _walk_import_time(node: ast.AST):
    """Like ast.walk, but skips lambda bodies (they run later, not at import)."""
    todo = [node]
    while todo:
        n = todo.pop()
        yield n
        if isinstance(n, ast.Lambda):
            todo.extend(n.args.defaults)
            todo.extend(d for d in n.args.kw_defaults if d is not None)
        else:
            todo.extend(ast.iter_child_nodes(n))


def _check_dynamic(tree: ast.Module) -> None:
    for n in ast.walk(tree):
        if isinstance(n, ast.Name) and n.id in DYNAMIC_NAMES:
            raise Unsupported(f"functions.py uses `{n.id}` (line {n.lineno})")
        if isinstance(n, ast.Attribute) and n.attr in DYNAMIC_ATTRS:
            raise Unsupported(f"functions.py uses `.{n.attr}` (line {n.lineno})")


def _import_name(node: ast.Import | ast.ImportFrom) -> list[str]:
    if isinstance(node, ast.Import):
        return [a.name for a in node.names]
    return ["." * node.level + (node.module or "")]


def analyse(source: str) -> Registry:
    """Split functions.py into shared code and per-dataset code."""
    tree = ast.parse(source)
    _check_dynamic(tree)
    shared: list[ast.AST] = []
    bindings: dict[str, ast.stmt] = {}
    for i, node in enumerate(tree.body):
        if i == 0 and isinstance(node, ast.Expr) and isinstance(node.value, ast.Constant):
            continue  # module docstring
        if isinstance(node, ast.Import | ast.ImportFrom):
            for module in _import_name(node):
                if module not in IMPORTS:
                    raise Unsupported(f"functions.py imports {module} (line {node.lineno})")
            shared.append(node)
            continue
        name = None
        if isinstance(node, ast.FunctionDef):
            name = node.name
        elif (
            isinstance(node, ast.Assign)
            and len(node.targets) == 1
            and isinstance(node.targets[0], ast.Name)
        ):
            name = node.targets[0].id
        elif (
            isinstance(node, ast.AnnAssign)
            and isinstance(node.target, ast.Name)
            and node.value is not None
        ):
            name = node.target.id
        if name is None:
            shared.append(node)
        elif name in bindings:
            raise Unsupported(f"functions.py binds `{name}` twice (line {node.lineno})")
        else:
            bindings[name] = node

    deps = {
        name: {n.id for n in ast.walk(node) if isinstance(n, ast.Name) and n.id in bindings}
        - {name}
        for name, node in bindings.items()
    }

    def closure(roots: set[str]) -> set[str]:
        seen: set[str] = set()
        todo = list(roots)
        while todo:
            name = todo.pop()
            if name not in seen:
                seen.add(name)
                todo.extend(deps[name])
        return seen

    # Definitions that call anything but Spec at import run code in every
    # process; they, and everything they reach, are shared.
    runs_at_import = set()
    for name, node in bindings.items():
        if isinstance(node, ast.FunctionDef) and node.decorator_list:
            runs_at_import.add(name)  # a decorator is called at import
        for part in _import_time_parts(node):
            for n in _walk_import_time(part):
                if isinstance(n, ast.Call) and not (
                    isinstance(n.func, ast.Name) and n.func.id in IMPORT_TIME_CALLS
                ):
                    runs_at_import.add(name)
    shared_names = closure(runs_at_import)
    shared.extend(bindings[name] for name in sorted(shared_names))

    if "ALL" not in bindings:
        raise Unsupported("functions.py defines no ALL")
    registry = bindings["ALL"].value
    if not isinstance(registry, ast.Tuple | ast.List):
        raise Unsupported("ALL is not a literal tuple or list")
    datasets: dict[str, str] = {}
    order: list[str] = []
    for element in registry.elts:
        binding = bindings.get(element.id) if isinstance(element, ast.Name) else None
        spec = binding.value if isinstance(binding, ast.Assign) else None
        if not (
            isinstance(spec, ast.Call)
            and isinstance(spec.func, ast.Name)
            and spec.func.id == "Spec"
        ):
            raise Unsupported(
                f"ALL entry `{ast.unparse(element)}` is not a name bound to Spec(...)"
            )
        names = [k.value for k in spec.keywords if k.arg == "name"]
        if len(names) != 1 or not (
            isinstance(names[0], ast.Constant) and isinstance(names[0].value, str)
        ):
            raise Unsupported(f"Spec `{element.id}` has no literal name=")
        dataset = names[0].value
        if dataset in datasets:
            raise Unsupported(f"dataset {dataset} is listed twice")
        reached = sorted(closure({element.id}) - shared_names)
        datasets[dataset] = _fingerprint([bindings[n] for n in reached])
        order.append(dataset)
    return Registry(_fingerprint(shared), datasets, order)


def _package_imports(text: str) -> set[str]:
    """The rsv_reference modules a module imports, by name ("spec", ...)."""
    found = set()
    for node in ast.walk(ast.parse(text)):
        if isinstance(node, ast.ImportFrom):
            module = node.module or ""
            if node.level == 0 and module.split(".")[0] != "rsv_reference":
                continue
            parts = module.split(".")[1:] if node.level == 0 else module.split(".")
            if parts and parts[0]:
                found.add(parts[0])
            else:
                found.update(a.name for a in node.names)  # from . import x
        elif isinstance(node, ast.Import):
            for a in node.names:
                parts = a.name.split(".")
                if parts[0] == "rsv_reference":
                    found.add(parts[1] if len(parts) > 1 else "functions")
    return found


def _imports_functions_back(functions: str, package: dict[str, str]) -> str | None:
    """A package module that functions.py imports (transitively) and that
    imports functions.py itself, giving its code a way to reach definitions by
    string; such access is invisible to the analysis. Modules functions.py
    does not import (the command-line driver) may import it freely."""
    modules = {Path(path).stem: text for path, text in package.items()}
    seen: set[str] = set()
    todo = sorted(_package_imports(functions))
    while todo:
        name = todo.pop()
        if name in seen or name not in modules:
            continue
        seen.add(name)
        imported = _package_imports(modules[name])
        if "functions" in imported:
            return f"{PACKAGE}{name}.py"
        todo.extend(sorted(imported))
    return None


def _affects_everything(path: str) -> bool:
    if path in EVERYTHING:
        return True
    if not path.startswith("reference/") or path == FUNCTIONS:
        return False
    return not path.startswith((DATA, TESTS))


def select(changed: list[str], base: str | None, head: str, package: dict[str, str]) -> Selection:
    """Decide what to regenerate.

    changed: paths that differ between base and head; base, head: functions.py
    at each (base is None if it did not exist); package: the other modules of
    rsv_reference at head, by path.
    """

    def everything(reason: str) -> Selection:
        return Selection(True, [f"regenerate all: {reason}"])

    for path in changed:
        if _affects_everything(path):
            return everything(f"{path} changed")
    if base is None:
        return everything(f"{FUNCTIONS} is new")
    try:
        cycle = _imports_functions_back(head, package)
        if cycle:
            return everything(f"{cycle} imports functions.py, which imports it")
        head_registry = analyse(head)
        base_registry = analyse(base)
    except Unsupported as e:
        return everything(str(e))
    if base_registry.shared != head_registry.shared:
        return everything(
            "shared code in functions.py changed (imports, module-level statements, "
            "or definitions that run at import)"
        )
    names, report = [], []
    for name in head_registry.order:
        reason = None
        if name not in base_registry.datasets:
            reason = "new dataset"
        elif base_registry.datasets[name] != head_registry.datasets[name]:
            reason = "its Spec or code it uses changed"
        elif f"{DATA}{name}.jsonl" in changed:
            reason = "its data file changed"
        if reason:
            names.append(name)
            report.append(f"regenerate {name}: {reason}")
        else:
            report.append(f"skip {name}: Spec, the code it uses and its data file are unchanged")
    return Selection(False, report, names)


def _git(repo: Path, *args: str) -> str:
    return subprocess.run(
        ["git", "-C", str(repo), *args], check=True, capture_output=True, text=True
    ).stdout


def _show(repo: Path, rev: str, path: str) -> str | None:
    r = subprocess.run(
        ["git", "-C", str(repo), "show", f"{rev}:{path}"],
        check=False,
        capture_output=True,
        text=True,
    )
    return r.stdout if r.returncode == 0 else None


def main(argv: list[str]) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--base", required=True, help="commit to compare with (the PR's base)")
    p.add_argument("--head", default="HEAD")
    p.add_argument("--repo", type=Path, default=Path(__file__).resolve().parent.parent)
    p.add_argument(
        "--registry",
        required=True,
        help="dataset names in ALL order, as the imported rsv_reference package reports them",
    )
    args = p.parse_args(argv)
    repo, head_rev, registry = args.repo, args.head, args.registry.split()

    head = _show(repo, head_rev, FUNCTIONS)
    if head is None:
        print(f"error: {FUNCTIONS} does not exist at {head_rev}", file=sys.stderr)
        return 1
    # Never skip silently: every data file must belong to a dataset and every
    # dataset must have a file, and the datasets the analysis sees must be
    # exactly the ones the package really generates.
    files = {
        Path(path).stem
        for path in _git(repo, "ls-tree", "--name-only", head_rev, DATA).split()
        if path.endswith(".jsonl")
    }
    if files != set(registry):
        print(
            "error: data files and datasets differ: "
            f"files without a dataset {sorted(files - set(registry))}, "
            f"datasets without a file {sorted(set(registry) - files)}",
            file=sys.stderr,
        )
        return 1
    try:
        order = analyse(head).order
    except Unsupported:
        order = registry  # select() regenerates everything in this case
    if order != registry:
        print(
            "error: ALL in functions.py, read as source, differs from the imported registry:\n"
            f"  source:   {' '.join(order)}\n  imported: {' '.join(registry)}",
            file=sys.stderr,
        )
        return 1

    changed = _git(repo, "diff", "--name-only", "--no-renames", args.base, head_rev).split()
    package = {
        path: _show(repo, head_rev, path) or ""
        for path in _git(repo, "ls-tree", "--name-only", head_rev, PACKAGE).split()
        if path.endswith(".py") and path != FUNCTIONS
    }
    selection = select(changed, _show(repo, args.base, FUNCTIONS), head, package)
    for line in selection.report:
        print(line, file=sys.stderr)
    if selection.everything:
        print("all")
    elif selection.names:
        print(" ".join(selection.names))
    else:
        print("none")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
