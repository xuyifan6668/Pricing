"""Vol screen validation tests (pure Python; run in CI without a C++ build)."""
import copy
import datetime as dt
import json

import pytest

from ratelib.demos import repo_root
from ratelib.vol import VolQuoteValidationError, load_vol_quotes, tenor_years, validate_vol_quotes
from ratelib.vol.vol_quotes import default_path, parse_vol_quotes, summary

PATH = repo_root() / "data" / "vol" / "usd_swaption_cube_2026-09-04.json"


def _doc():
    with open(PATH, encoding="utf-8") as f:
        return json.load(f)


def test_tenor_years_matches_cube_axis():
    assert tenor_years("1M") == pytest.approx(1 / 12)
    assert tenor_years("10Y") == 10.0
    assert tenor_years("2W") == pytest.approx(14 / 365)
    with pytest.raises(ValueError):
        tenor_years("0Y")
    with pytest.raises(ValueError):
        tenor_years("ten years")


def test_shipped_file_is_valid_and_complete():
    assert default_path() == PATH
    q = load_vol_quotes(PATH)
    assert q.as_of == dt.date(2026, 9, 4)
    assert q.index == "USD-SOFR"
    assert len(q.nodes) == 81
    assert q.expiries()[0] == "1M" and q.expiries()[-1] == "10Y"
    assert q.tenors()[0] == "1Y" and q.tenors()[-1] == "30Y"
    atm = q.atm_matrix()
    assert atm.shape == (9, 9)
    assert list(atm.index) == q.expiries()
    assert atm.loc["1M", "1Y"] > atm.loc["10Y", "30Y"]  # vols fall with expiry and tenor
    frame = q.to_frame()
    assert len(frame) == 81 * 7
    assert set(frame["offset_bp"]) == {-100.0, -50.0, -25.0, 0.0, 25.0, 50.0, 100.0}
    warnings = validate_vol_quotes(q, today=dt.date(2026, 9, 4))
    assert not any("days old" in w for w in warnings)
    assert not any("concave" in w for w in warnings)
    node = q.find("10Y", "30Y")
    assert node.id == "10Yx30Y"
    assert node.points()[3] == (0.0, node.atm_bp)
    assert "81 nodes" in summary(q)


def test_stale_and_concave_warnings():
    q = load_vol_quotes(PATH, validate=False)
    assert any("days old" in w for w in validate_vol_quotes(q, today=dt.date(2026, 10, 1)))
    d = _doc()
    node = d["nodes"][40]
    for p in node["smile"]:
        if abs(p["offset_bp"]) == 100:
            p["vol_bp"] = node["atm_bp"] - 5.0
    w = validate_vol_quotes(parse_vol_quotes(d), today=dt.date(2026, 9, 4))
    assert any("concave" in x for x in w)


@pytest.mark.parametrize(
    "mutate, fragment",
    [
        (lambda d: d.update(unit="lognormal_pct"), "not supported"),
        (lambda d: d["nodes"].append(copy.deepcopy(d["nodes"][0])), "duplicate node"),
        (lambda d: d["nodes"].pop(10), "grid has no node"),
        (lambda d: d["nodes"][0].update(atm_bp=0.0), "outside (0, 1000)"),
        (lambda d: d["nodes"][0]["smile"][0].update(vol_bp=1500.0), "outside (0, 1000)"),
        (lambda d: d["nodes"][0]["smile"][0].update(offset_bp=0), "offset 0"),
        (lambda d: d["nodes"][0]["smile"].reverse(), "not strictly increasing"),
        (lambda d: d["nodes"][0].update(forward_pct=350.0), "wrong unit"),
        (lambda d: d["nodes"][0].update(expiry="1Q"), "bad tenor"),
    ],
)
def test_hard_errors(mutate, fragment):
    d = _doc()
    mutate(d)
    with pytest.raises(VolQuoteValidationError) as e:
        validate_vol_quotes(parse_vol_quotes(d), today=dt.date(2026, 9, 4))
    assert fragment in str(e.value)


def test_missing_keys_are_structural_errors():
    d = _doc()
    del d["nodes"][0]["atm_bp"]
    with pytest.raises(VolQuoteValidationError):
        parse_vol_quotes(d)
