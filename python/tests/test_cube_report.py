"""Cube report tests on hand-written CSV (pure Python)."""
import pytest

from ratelib.vol import CubeReport, load_arbitrage_csv, load_fit_csv
from ratelib.vol.cube_builder import USD_POLICY, SabrPolicy

FIT = """id,expiry,tenor,forward_pct,atm_bp,alpha,beta,rho,nu,rms_bp,max_bp,iterations,status,ok,flags
1Yx10Y,1Y,10Y,3.6012,88.01,0.04591,0.5,-0.1331,0.4127,2.1e-03,4.4e-03,7,CostSmall,1,
1Yx2Y,1Y,2Y,3.3804,95.83,0.05102,0.5,-0.1004,0.4501,1.9e-03,3.8e-03,6,CostSmall,1,
10Yx2Y,10Y,2Y,3.6421,80.31,0.04312,0.5,-0.2426,0.2337,2.5e-03,5.0e-03,9,CostSmall,1,
10Yx10Y,10Y,10Y,3.7610,72.48,0.03620,0.5,-0.2700,0.9999,4.0e-02,1.3e-01,40,MaxIterations,0,nu at bound;max residual 0.13 bp > 0.1 bp
"""

ARB = """id,expiry,tenor,forward_pct,atm_bp,k_min_pct,k_max_pct,min_density,first_negative_below_pct,first_negative_above_pct,call_spread_ok,nonpositive_vol,arbitrage_free
1Yx10Y,1Y,10Y,3.6012,88.01,0.08,7.12,1.2e-01,,,1,0,1
1Yx2Y,1Y,2Y,3.3804,95.83,0.01,7.21,3.1e-02,,,1,0,1
10Yx2Y,10Y,2Y,3.6421,80.31,0.01,13.8,-4.0e+00,0.412,,1,0,0
10Yx10Y,10Y,10Y,3.7610,72.48,0.01,12.9,-2.1e+04,2.101,,1,0,0
"""


def test_statistics_flags_and_verdict():
    rep = load_fit_csv(FIT)
    assert len(rep.fits) == 4
    assert rep.max_abs_bp == pytest.approx(0.13)
    assert not rep.passed
    assert list(rep.failed_nodes()["id"]) == ["10Yx10Y"]
    assert rep.worst(1)["id"].iloc[0] == "10Yx10Y"
    assert "FAIL" in rep.verdict() and "1 flagged" in rep.verdict()
    assert rep.fits.loc[0, "flags"] == ""


def test_pivot_is_in_axis_order_not_string_order():
    rep = load_fit_csv(FIT)
    piv = rep.pivot("max_bp")
    assert list(piv.index) == ["1Y", "10Y"]
    assert list(piv.columns) == ["2Y", "10Y"]
    assert piv.loc["10Y", "10Y"] == pytest.approx(0.13)


def test_negative_density_map():
    rep = load_fit_csv(FIT, arbitrage=load_arbitrage_csv(ARB))
    assert rep.arbitrage_nodes() == ["10Yx2Y", "10Yx10Y"]
    m = rep.negative_density_map()
    assert m.loc["10Y", "10Y"] == pytest.approx(2.101)
    assert m.isna().loc["1Y", "10Y"]
    assert "2 nodes with a negative density" in rep.verdict()
    with pytest.raises(ValueError):
        load_fit_csv(FIT).negative_density_map()


def test_pass_markdown_and_round_trip(tmp_path):
    good = "\n".join(FIT.splitlines()[:4]) + "\n"
    rep = load_fit_csv(good)
    assert rep.passed
    assert "PASS" in rep.verdict()
    md = rep.to_markdown()
    assert md.startswith("| id | forward_pct |")
    rep.write(tmp_path / "cube.md")
    assert (tmp_path / "cube.md").read_text(encoding="utf-8").strip().endswith("0 flagged")
    rep.write(tmp_path / "cube.csv")
    assert len(load_fit_csv(tmp_path / "cube.csv").fits) == 3


def test_missing_columns_and_policy_args():
    with pytest.raises(ValueError):
        load_fit_csv("id,max_bp\nA,1\n")
    with pytest.raises(ValueError):
        load_arbitrage_csv("id,arbitrage_free\nA,1\n")
    assert USD_POLICY == SabrPolicy(0.5, 0.0)
    assert SabrPolicy(0.0, 0.02).args() == ["0.0", "0.02"]
    assert isinstance(load_fit_csv(FIT), CubeReport)
