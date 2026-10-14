"""Tests for the weekly-comparison, pin-relaxing and mutation-summary scripts."""

import json
import subprocess
from pathlib import Path

import mutation_score
import relax_pins
import weekly_compare

HEADER = '{"schema":1}\n'


def test_relax_pins_only_touches_exact_pins():
    text = 'puruspe = "=0.4.4"\nharness.workspace = true\nserde = "1"\n'
    assert relax_pins.relax(text) == 'puruspe = "*"\nharness.workspace = true\nserde = "1"\n'


def report(dir_: Path, name: str, row: str) -> None:
    dir_.mkdir(exist_ok=True)
    (dir_ / name).write_text(
        f"# x\n\n| Function | ... |\n|---|\n| `p::f` | {row} |\n| `p::g` | same |\n"
    )


def test_report_changes_matches_crates_across_versions(tmp_path):
    old, new = tmp_path / "old", tmp_path / "new"
    report(old, "p-0.4.4.md", "361 | 317 | 0")
    report(new, "p-0.5.0.md", "361 | 316 | 1")
    lines = weekly_compare.report_changes(old, new)
    assert lines[0] == "### p: committed 0.4.4, latest 0.5.0"
    assert "| `p::f` | 361 | 317 | 0 | 361 | 316 | 1 |" in lines
    assert not any("p::g" in ln for ln in lines)
    report(new, "p-0.5.0.md", "361 | 317 | 0")
    assert weekly_compare.report_changes(old, new) == []
    (new / "notes.txt").write_text("ignored")
    report(new, "q-1.0.0.md", "1 | 1 | 0")
    assert weekly_compare.report_changes(old, new)[0] == "### q: committed none, latest 1.0.0"


def test_test_problems():
    assert weekly_compare.test_problems(0, "anything") == []
    log = "noise\nt::f[1.0]: unexpected CrateWrong\nt::f[2.0]: listed in 001 but no longer fails\n"
    out = weekly_compare.test_problems(101, log)
    assert out[0] == "`cargo test` failed (exit 101)."
    assert out[2:] == [
        "- `t::f[1.0]: unexpected CrateWrong`",
        "- `t::f[2.0]: listed in 001 but no longer fails`",
    ]
    build = weekly_compare.test_problems(101, "a\nerror[E0425]: cannot find function\n")
    assert build == [
        "`cargo test` failed (exit 101). Last lines:",
        "",
        "```",
        "a",
        "error[E0425]: cannot find function",
        "```",
    ]
    many = "\n".join(f"unexpected {i}" for i in range(35))
    assert weekly_compare.test_problems(1, many)[-1] == "- ... and 5 more"


def git(root: Path, *args: str) -> None:
    subprocess.run(["git", *args], cwd=root, check=True, capture_output=True)


def point(mp: float, sp: float) -> str:
    return json.dumps({"in": {"x": [1.0]}, "mpmath": {"y": [mp]}, "scipy": {"y": [sp]}}) + "\n"


def test_tool_versions_without_data(tmp_path):
    (tmp_path / "reference/data").mkdir(parents=True)
    assert weekly_compare.tool_versions(tmp_path) == []


def test_data_changes_against_head(tmp_path):
    data = tmp_path / "reference/data"
    data.mkdir(parents=True)
    (data / "a.jsonl").write_text(HEADER + point(1, 1) + point(2, 2))
    (data / "b.jsonl").write_text(HEADER + point(1, 1))
    git(tmp_path, "init", "-q")
    git(tmp_path, "add", "-A")
    git(tmp_path, "-c", "user.name=t", "-c", "user.email=t@t", "commit", "-qm", "x")
    assert weekly_compare.data_changes(tmp_path) == []
    assert weekly_compare.tool_versions(tmp_path)[0] == "- scipy: committed ?, latest ?"
    (data / "a.jsonl").write_text(HEADER + point(1, 1) + point(2, 3))
    (data / "b.jsonl").write_text(HEADER)
    (data / "c.jsonl").write_text(HEADER)
    (data / "a.jsonl").write_text(
        '{"scipy":{"version":"1.19.0"},"numpy":"2.6.0"}\n' + point(1, 1) + point(2, 3)
    )
    assert weekly_compare.tool_versions(tmp_path)[:2] == [
        "- scipy: committed ?, latest 1.19.0",
        "- numpy: committed ?, latest 2.6.0",
    ]
    assert weekly_compare.data_changes(tmp_path) == [
        "- `reference/data/a.jsonl`: 0 mpmath value(s) and 1 SciPy value(s) changed",
        "- `reference/data/b.jsonl`: point count 1 -> 0",
        "- `reference/data/c.jsonl`: new dataset",
    ]


def outcomes(caught: int, missed: int, timeout: int = 0) -> dict:
    missed_list = [
        {"scenario": {"Mutant": {"name": f"m{i}"}}, "summary": "MissedMutant"}
        for i in range(missed)
    ]
    return {
        "caught": caught,
        "missed": missed,
        "timeout": timeout,
        "unviable": 3,
        "outcomes": [{"scenario": "Baseline", "summary": "Success"}, *missed_list],
    }


def test_mutation_summary():
    text = mutation_score.summary(outcomes(3, 1))
    assert "**Mutation score: 75.0%** (3 caught of 4 testable mutants)" in text
    assert "| 3 | 1 | 0 | 3 |" in text
    assert text.endswith("- `m0`\n")
    assert "n/a" in mutation_score.summary(outcomes(0, 0))
    assert "did not catch" not in mutation_score.summary(outcomes(2, 0))


def test_mutation_main(tmp_path, capsys):
    assert mutation_score.main([str(tmp_path)]) == 1
    (tmp_path / "outcomes.json").write_text(json.dumps(outcomes(1, 0)))
    assert mutation_score.main([str(tmp_path)]) == 0
    assert "100.0%" in capsys.readouterr().out


def test_violations():
    assert mutation_score.violations(outcomes(3, 0), 0, None) == []
    assert mutation_score.violations(outcomes(3, 1), 0, None) == [
        "1 surviving mutant(s); at most 0 allowed"
    ]
    assert mutation_score.violations(outcomes(3, 0, timeout=1), 0, None) != []
    assert mutation_score.violations(outcomes(19, 1), None, 95) == []
    assert mutation_score.violations(outcomes(18, 2), None, 95) == [
        "mutation score 90.0% is below 95%"
    ]
    assert mutation_score.violations(outcomes(0, 0), 0, 95) == []


def test_main_enforces_rule(tmp_path, capsys):
    (tmp_path / "outcomes.json").write_text(json.dumps(outcomes(3, 1)))
    assert mutation_score.main([str(tmp_path)]) == 0
    assert mutation_score.main([str(tmp_path), "--max-survivors", "0"]) == 1
    assert "**FAIL:** 1 surviving" in capsys.readouterr().out


def test_exclusions_need_reasons(tmp_path):
    ok = 'exclude_re = [\n  # reason: equivalent mutant\n  "report.rs:15",\n]\n'
    bad = 'exclude_re = [\n  # just a note\n  "a",\n  # reason: fine\n  "b",\n  "c",\n]\n'
    assert mutation_score.unexplained_exclusions("exclude_re = []\n") == []
    assert mutation_score.unexplained_exclusions("# nothing\n") == []
    assert mutation_score.unexplained_exclusions(ok) == []
    assert mutation_score.unexplained_exclusions(bad) == ['"a"', '"c"']
    cfg = tmp_path / "mutants.toml"
    cfg.write_text(bad)
    assert mutation_score.main(["--check-config", str(cfg)]) == 1
    cfg.write_text(ok)
    assert mutation_score.main(["--check-config", str(cfg)]) == 0
