"""Quote validation tests (pure Python; run in CI without a C++ build)."""
import copy
import datetime as dt
import json

import pytest

from ratelib.demos import repo_root
from ratelib.market import QuoteValidationError, load_quotes, validate_quotes
from ratelib.market.quotes import parse_quotes, summary


def _doc():
    with open(repo_root() / "data" / "quotes" / "usd_2026-09-04.json", encoding="utf-8") as f:
        return json.load(f)


def test_shipped_file_is_valid():
    qs = load_quotes(repo_root() / "data" / "quotes" / "usd_2026-09-04.json")
    assert qs.as_of == dt.date(2026, 9, 4)
    assert len(qs.quotes) >= 35
    w = validate_quotes(qs, today=dt.date(2026, 9, 4))
    # The cross-assigned 25Y/30Y basis swaps are flagged as such (a warning, by design).
    assert any("FF-SOFR-30Y" in x and "cross-assignment" in x for x in w)
    assert not any("stale" in x or "days old" in x for x in w)
    assert qs.find("FF-SOFR-30Y").curve_name == "USD-SOFR"
    assert qs.find("FF-SOFR-1Y").curve_name == "USD-FEDFUNDS"
    assert "USD-SOFR" in summary(qs)


def test_stale_warning():
    qs = load_quotes(repo_root() / "data" / "quotes" / "usd_2026-09-04.json", validate=False)
    w = validate_quotes(qs, today=dt.date(2026, 10, 1))
    assert any("days old" in x for x in w)


@pytest.mark.parametrize(
    "mutate, fragment",
    [
        (lambda d: d["quotes"].append(copy.deepcopy(d["quotes"][0])), "duplicate id"),
        (lambda d: d["quotes"][0].update(type="swaption"), "unknown type"),
        (lambda d: d["quotes"][0].update(unit="ticks"), "unknown unit"),
        (lambda d: d["quotes"][0].update(value=385.0), "wrong unit"),
        (lambda d: next(q for q in d["quotes"] if q["type"] == "future").update(value=3.7, unit="pct"), "quoted as prices"),
        (lambda d: next(q for q in d["quotes"] if q["type"] == "future").update(start="2026-12-15"), "not an IMM date"),
        (lambda d: next(q for q in d["quotes"] if q["type"] == "basis_swap").update(reference_index="USD-FEDFUNDS"), "same index"),
        (lambda d: next(q for q in d["quotes"] if q["id"] == "SOFR-OIS-5Y").update(tenor="3Y"), "pin the same pillar"),
        (lambda d: next(q for q in d["quotes"] if q["id"] == "SR3-H7").update(start="2027-06-16", end="2027-09-15"), "strip gap"),
        (lambda d: d["quotes"][0].pop("tenor"), "needs a tenor"),
    ],
)
def test_hard_errors(mutate, fragment):
    d = _doc()
    mutate(d)
    qs = parse_quotes(d)
    with pytest.raises(QuoteValidationError) as e:
        validate_quotes(qs, today=dt.date(2026, 9, 4))
    assert fragment in str(e.value)


def test_all_errors_reported_at_once():
    d = _doc()
    d["quotes"][0]["type"] = "swaption"
    d["quotes"][1]["unit"] = "ticks"
    with pytest.raises(QuoteValidationError) as e:
        validate_quotes(parse_quotes(d), today=dt.date(2026, 9, 4))
    assert len(str(e.value).splitlines()) == 2


def test_rate_conversion():
    d = _doc()
    qs = parse_quotes(d)
    assert qs.find("SOFR-OIS-10Y").rate == pytest.approx(0.0345)
    assert qs.find("FF-SOFR-1Y").rate == pytest.approx(-2.0e-4)
    with pytest.raises(ValueError):
        _ = qs.find("SR3-Z6").rate
