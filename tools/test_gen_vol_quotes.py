"""Pin the vol screen generator: formula limits, the ATM inversion, the stressed node, and the committed file.

Run: uv run pytest tools/
"""
import json

import pytest

import gen_vol_quotes as g


def test_axis_years():
    assert g.years("1M") == pytest.approx(1 / 12)
    assert g.years("10Y") == 10.0
    assert g.years("1W") == pytest.approx(7 / 365)


def test_beta_zero_atm_limit():
    a, rho, nu, t = 0.0085, -0.4, 0.5, 3.0
    expected = a * (1 + (2 - 3 * rho * rho) * nu * nu * t / 24)
    assert g.hagan_normal_vol(a, 0.0, rho, nu, 0.0, 0.03, 0.03, t) == pytest.approx(expected, rel=1e-15)


def test_rho_zero_symmetry_and_nu_zero_limit():
    for beta in (0.0, 0.5, 1.0):
        a = {0.0: 0.009, 0.5: 0.047, 1.0: 0.25}[beta]
        assert g.hagan_normal_vol(a, beta, 0.0, 0.6, 0.0, 0.03, 0.042, 2.0) == pytest.approx(g.hagan_normal_vol(a, beta, 0.0, 0.6, 0.0, 0.042, 0.03, 2.0), rel=1e-14)
    f, k, t, a, b = 0.035, 0.025, 2.0, 0.05, 0.5
    c = 1 - b
    lead = a * c * (f - k) / (f**c - k**c)
    corr = -b * (2 - b) * a * a / (24 * (f * k) ** c)
    assert g.hagan_normal_vol(a, b, -0.3, 0.0, 0.0, f, k, t) == pytest.approx(lead * (1 + corr * t), rel=1e-13)


def test_continuity_through_the_money():
    for beta in (0.0, 0.5, 1.0):
        a = {0.0: 0.009, 0.5: 0.047, 1.0: 0.25}[beta]
        v0 = g.hagan_normal_vol(a, beta, -0.25, 0.5, 0.0, 0.036, 0.036, 1.0)
        assert abs(g.hagan_normal_vol(a, beta, -0.25, 0.5, 0.0, 0.036, 0.036 + 1e-9, 1.0) - v0) < 1e-9
        assert g.atm_normal_vol(a, beta, -0.25, 0.5, 0.0, 0.036, 1.0) == pytest.approx(v0, rel=1e-14)


def test_alpha_from_atm_inverts():
    for beta in (0.0, 0.5, 1.0):
        for rho in (-0.7, 0.0, 0.5):
            for nu in (0.0, 0.4, 1.0):
                a = g.alpha_from_atm(0.0085, 0.036, 5.0, beta, rho, nu, 0.0)
                assert g.atm_normal_vol(a, beta, rho, nu, 0.0, 0.036, 5.0) == pytest.approx(0.0085, rel=1e-13)
    with pytest.raises(ValueError):
        g.alpha_from_atm(0.0150, 0.03, 30.0, 1.0, 0.0, 0.0, 0.0)


def test_document_shape_and_stressed_node():
    doc = g.build_document()
    assert doc["asOf"] == "2026-09-04" and doc["unit"] == "normal_bp"
    assert len(doc["nodes"]) == 81
    for n in doc["nodes"]:
        assert [p["offset_bp"] for p in n["smile"]] == g.OFFSETS_BP
        assert 2.0 < n["forward_pct"] < 6.0
        assert 40.0 < n["atm_bp"] < 200.0
        assert all(n["forward_pct"] / 100 + p["offset_bp"] * 1e-4 > 0 for p in n["smile"])
    stressed = g.node_parameters(*g.STRESSED)
    assert stressed["nu"] == g.STRESSED_NU and stressed["rho"] == g.STRESSED_RHO
    benign = g.node_parameters("10Y", "20Y")
    assert benign["nu"] < 0.3


def test_committed_file_matches_the_generator():
    committed = json.loads(g.OUT.read_text(encoding="utf-8"))
    assert committed == g.build_document(), "data/vol file is stale: run `uv run tools/gen_vol_quotes.py`"


def test_every_smile_is_a_smile_not_a_frown():
    # The wings average above the ATM vol on every node (SABR with nu > 0 is convex near the money). The sign of
    # the skew is NOT pinned: at short expiries the beta = 0.5 backbone and a large nu lift the payer wing.
    for n in g.build_document()["nodes"]:
        vols = {p["offset_bp"]: p["vol_bp"] for p in n["smile"]}
        assert 0.5 * (vols[-100] + vols[100]) > n["atm_bp"], (n["expiry"], n["tenor"])
        assert 0.5 * (vols[-25] + vols[25]) > n["atm_bp"] - 0.01, (n["expiry"], n["tenor"])
