"""Package-level smoke tests that do not need a C++ build."""
from pathlib import Path

import pytest

import ratelib
from ratelib.demos import repo_root


def test_version():
    assert ratelib.__version__ == "0.1.0"


def test_repo_root_has_data():
    assert (repo_root() / "data" / "holidays" / "us_sifma.txt").exists()
    assert (repo_root() / "data" / "sobol" / "new-joe-kuo-6.21201.txt").exists()


def test_find_build_dir_env(monkeypatch, tmp_path: Path):
    monkeypatch.setenv("RL_BUILD_DIR", str(tmp_path))
    assert ratelib.find_build_dir() == tmp_path


def test_run_demo_reports_missing_build(monkeypatch, tmp_path: Path):
    monkeypatch.setenv("RL_BUILD_DIR", str(tmp_path))
    with pytest.raises(FileNotFoundError):
        ratelib.run_demo("rl_demo_schedule")


def test_run_analysis_needs_a_section_and_a_build(monkeypatch, tmp_path: Path):
    monkeypatch.setenv("RL_BUILD_DIR", str(tmp_path))
    with pytest.raises(ValueError):
        ratelib.run_analysis()
    with pytest.raises(FileNotFoundError):
        ratelib.run_analysis("xva")


def test_run_analysis_reuses_tables_already_written(monkeypatch, tmp_path: Path):
    monkeypatch.setenv("RL_BUILD_DIR", str(tmp_path))
    out = tmp_path / "notebook-results" / "xva"
    out.mkdir(parents=True)
    (out / "xva_table.csv").write_text("case,cva_bp\nnetted,15.2\n")
    tables = ratelib.run_analysis("xva")
    assert list(tables) == ["xva_table"]
    assert tables["xva_table"]["cva_bp"][0] == 15.2
