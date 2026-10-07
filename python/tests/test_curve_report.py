"""Curve report tests on a hand-written CSV (pure Python)."""
import pytest

from ratelib.market import load_reprice_csv

CSV = """id,type,curve,maturity,quote,model,residual,residual_bp,raw_residual
SOFR-OIS-1Y,ois_swap,USD-SOFR,2027-09-09,3.48,3.4800000000001,1e-13,1e-11,2e-15
SOFR-OIS-10Y,ois_swap,USD-SOFR,2036-09-09,3.45,3.45,0,0,0
SR3-Z6,future,USD-SOFR,2027-03-17,96.33,96.330002,2e-6,2e-4,1e-8
FF-SOFR-30Y,basis_swap,USD-SOFR,2056-09-11,-4.0,-3.9999,1e-4,1e-4,3e-7
"""


def test_load_and_stats():
    rep = load_reprice_csv(CSV, tolerance_bp=1e-8)
    assert len(rep.table) == 4
    assert rep.max_abs_bp == pytest.approx(2e-4)
    assert not rep.passed
    assert "FAIL" in rep.verdict()
    assert list(rep.worst(2)["id"]) == ["SR3-Z6", "FF-SOFR-30Y"]
    bc = rep.by_curve()
    assert bc.loc["USD-SOFR", "instruments"] == 4
    md = rep.to_markdown()
    assert md.startswith("| id | type |")
    assert "FF-SOFR-30Y" in md


def test_pass_with_loose_tolerance(tmp_path):
    rep = load_reprice_csv(CSV, tolerance_bp=1e-3)
    assert rep.passed
    rep.write(tmp_path / "report.md")
    assert (tmp_path / "report.md").read_text(encoding="utf-8").endswith("bp)\n")
    rep.write(tmp_path / "report.csv")
    again = load_reprice_csv(tmp_path / "report.csv")
    assert len(again.table) == 4


def test_missing_columns():
    with pytest.raises(ValueError):
        load_reprice_csv("id,quote\nA,1\n")
