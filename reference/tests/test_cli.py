"""`python -m rsv_reference`: one fresh process per dataset, in parallel."""

from pathlib import Path

from rsv_reference.cli import DATA_DIR, main


def test_selected_datasets_equal_the_committed_files(tmp_path: Path, capsys):
    assert main(["gamma", "erfcx"], tmp_path) == 0
    assert sorted(p.name for p in tmp_path.iterdir()) == ["erfcx.jsonl", "gamma.jsonl"]
    for name in ("erfcx", "gamma"):
        assert (tmp_path / f"{name}.jsonl").read_bytes() == (
            DATA_DIR / f"{name}.jsonl"
        ).read_bytes()
    lines = capsys.readouterr().out.splitlines()
    # Printed in registry order, whatever order the processes finish in.
    assert [line.split(":")[0] for line in lines] == ["erfcx", "gamma"]
    assert lines[1].startswith("gamma: 901 points (") and lines[1].endswith(" s")


def test_unknown_name_is_an_error(tmp_path: Path, capsys):
    assert main(["erfcx", "nope"], tmp_path) == 2
    assert "unknown function(s): nope" in capsys.readouterr().err
    assert not list(tmp_path.iterdir())
