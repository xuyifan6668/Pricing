#!/usr/bin/env python3
"""Generate the synthetic USD swaption vol screen data/vol/usd_swaption_cube_2026-09-04.json.

Run:  uv run tools/gen_vol_quotes.py          (byte-identical output if the surface is unchanged)
Test: uv run pytest tools/

SIMPLIFICATION: synthetic marks, not vendor data. The screen is generated from smooth forward, ATM-vol, rho and
nu surfaces in a plausible September-2026 USD configuration (forwards 3.3-3.8%, ATM normal vols 70-106 bp, negative
skew growing with expiry, vol-of-vol decaying with expiry), beta = 0.5, no shift. Quotes are Hagan normal vols
(the same formula as cpp/vol/include/rl/vol/detail/sabr.ipp) rounded to 0.01 bp, forwards to 0.0001%, so a
correct per-node fit leaves residuals of a few thousandths of a bp: honest, not zero.

One node is deliberately stressed: 10Yx30Y has nu = 0.70 and rho = -0.45. It fits to the same precision as every
other node, and its implied density goes negative about 170 bp below the forward, inside one standard deviation.
That is the demonstration that a perfect fit and an arbitrage-free smile are different things
(docs/RESULTS.md section 3; notebook 07).

beta and the shift are NOT written into the quote file: they are model policy, chosen by the cube builder
(python/ratelib/vol/cube_builder.py) and passed to the fit.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

AS_OF = "2026-09-04"
EXPIRIES = ["1M", "3M", "6M", "1Y", "2Y", "3Y", "5Y", "7Y", "10Y"]
TENORS = ["1Y", "2Y", "3Y", "5Y", "7Y", "10Y", "15Y", "20Y", "30Y"]
OFFSETS_BP = [-100, -50, -25, 25, 50, 100]
BETA = 0.5
SHIFT = 0.0
STRESSED = ("10Y", "30Y")
STRESSED_RHO = -0.45
STRESSED_NU = 0.70
OUT = Path(__file__).resolve().parents[1] / "data" / "vol" / "usd_swaption_cube_2026-09-04.json"


def years(label: str) -> float:
    """Axis year fraction, as market::cubeAxisYears: months/12, days/365."""
    n, unit = int(label[:-1]), label[-1].upper()
    return {"D": n / 365.0, "W": 7.0 * n / 365.0, "M": n / 12.0, "Y": float(n)}[unit]


# ---------------------------------------------------------------------------------------------------------------
# Generating surfaces
# ---------------------------------------------------------------------------------------------------------------
def forward(e: float, t: float) -> float:
    return 0.0330 + 0.0010 * math.log1p(e) / math.log(11.0) + 0.0045 * (1.0 - math.exp(-(e + t) / 12.0)) - 0.0012 * math.exp(-t / 2.0)


def atm_vol(e: float, t: float) -> float:
    return 0.0068 + 0.0030 * math.exp(-e / 2.5) + 0.0012 * math.exp(-t / 4.0) - 0.0004 * math.log1p(t) / math.log(31.0) + 0.0006 * (1.0 - math.exp(-e / 3.0))


def rho_surface(e: float, t: float) -> float:
    return -0.05 - 0.20 * (1.0 - math.exp(-e / 4.0)) - 0.05 * (1.0 - math.exp(-t / 10.0))


def nu_surface(e: float, t: float) -> float:
    return 0.18 + 0.50 * math.exp(-e / 1.2) + 0.08 * math.exp(-t / 5.0)


# ---------------------------------------------------------------------------------------------------------------
# Hagan normal vol (B.69) with the integral-form zeta — a transcription of cpp/vol/include/rl/vol/detail/sabr.ipp
# ---------------------------------------------------------------------------------------------------------------
def _zeta_over_xhat(z: float, rho: float) -> float:
    if abs(z) < 1e-6:
        return 1.0 - 0.5 * rho * z + (2.0 - 3.0 * rho * rho) * z * z / 12.0
    sd = math.sqrt(1.0 - 2.0 * rho * z + z * z)
    if z >= 0.0:
        xhat = math.log1p(z * (z - 2.0 * rho + sd + 1.0) / ((sd + 1.0) * (1.0 - rho)))
    else:
        xhat = -math.log1p(z * (z - 2.0 * rho - sd - 1.0) / ((sd + 1.0) * (1.0 + rho)))
    return z / xhat


def hagan_normal_vol(alpha: float, beta: float, rho: float, nu: float, shift: float, f_: float, k_: float, t: float) -> float:
    if beta == 0.0:
        zeta = nu / alpha * (f_ - k_)
        return alpha * _zeta_over_xhat(zeta, rho) * (1.0 + (2.0 - 3.0 * rho * rho) * nu * nu / 24.0 * t)
    f, k = f_ + shift, k_ + shift
    if f <= 0.0 or k <= 0.0:
        raise ValueError("forward + shift and strike + shift must be positive when beta > 0")
    c = 1.0 - beta
    r = math.log(f / k)
    near_atm = abs(r) < 1e-8
    if beta == 1.0:
        ratio = 1.0 + r * (0.5 + r / 6.0) if near_atm else math.expm1(r) / r
        factor0, q = alpha * k * ratio, r
    else:
        ratio = (1.0 + r * (0.5 + r / 6.0)) / (1.0 + c * r * (0.5 + c * r / 6.0)) if near_atm else c * math.expm1(r) / math.expm1(c * r)
        factor0 = alpha * k**beta * ratio
        q = k**c * math.expm1(c * r) / c
    zeta = nu / alpha * q
    fav_c = math.exp(c * 0.5 * (math.log(f) + math.log(k)))
    corr = -beta * (2.0 - beta) * alpha * alpha / (24.0 * fav_c * fav_c) + rho * alpha * nu * beta / (4.0 * fav_c) + (2.0 - 3.0 * rho * rho) * nu * nu / 24.0
    return factor0 * _zeta_over_xhat(zeta, rho) * (1.0 + corr * t)


def atm_normal_vol(alpha: float, beta: float, rho: float, nu: float, shift: float, f_: float, t: float) -> float:
    if beta == 0.0:
        return alpha * (1.0 + (2.0 - 3.0 * rho * rho) * nu * nu * t / 24.0)
    f = f_ + shift
    fc = f ** (1.0 - beta)
    corr = -beta * (2.0 - beta) * alpha * alpha / (24.0 * fc * fc) + rho * alpha * nu * beta / (4.0 * fc) + (2.0 - 3.0 * rho * rho) * nu * nu / 24.0
    return alpha * f**beta * (1.0 + corr * t)


def alpha_from_atm(sigma: float, f_: float, t: float, beta: float, rho: float, nu: float, shift: float) -> float:
    """Smallest positive root of the ATM cubic, by bisection between 0 and its first stationary point (West 2005)."""
    f = f_ + shift
    fb = 1.0 if beta == 0.0 else f**beta
    fc = 1.0 if beta == 0.0 else f ** (1.0 - beta)
    a1 = fb * (1.0 + (2.0 - 3.0 * rho * rho) * nu * nu * t / 24.0)
    a2 = 0.0 if beta == 0.0 else fb * rho * nu * beta * t / (4.0 * fc)
    a3 = 0.0 if beta == 0.0 else -fb * beta * (2.0 - beta) * t / (24.0 * fc * fc)
    if a3 == 0.0:
        return sigma / a1
    disc = math.sqrt(4.0 * a2 * a2 - 12.0 * a3 * a1)
    a_star = 2.0 * a1 / (-2.0 * a2 + disc) if a2 < 0.0 else (-2.0 * a2 - disc) / (6.0 * a3)
    g = lambda a: ((a3 * a + a2) * a + a1) * a - sigma  # noqa: E731
    if g(a_star) < 0.0:
        raise ValueError("ATM vol unattainable")
    lo, hi = 0.0, a_star
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if g(mid) < 0.0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


# ---------------------------------------------------------------------------------------------------------------
# The screen
# ---------------------------------------------------------------------------------------------------------------
def node_parameters(expiry: str, tenor: str) -> dict:
    """The generating parameters of one node, after rounding the forward and ATM vol as quoted."""
    e, t = years(expiry), years(tenor)
    fwd = round(forward(e, t) * 100.0, 4) / 100.0
    atm = round(atm_vol(e, t) * 1e4, 2) / 1e4
    rho, nu = (STRESSED_RHO, STRESSED_NU) if (expiry, tenor) == STRESSED else (rho_surface(e, t), nu_surface(e, t))
    alpha = alpha_from_atm(atm, fwd, e, BETA, rho, nu, SHIFT)
    return {"expiry": expiry, "tenor": tenor, "expiry_years": e, "forward": fwd, "atm": atm, "alpha": alpha, "beta": BETA, "rho": rho, "nu": nu, "shift": SHIFT}


def build_document() -> dict:
    nodes = []
    for expiry in EXPIRIES:
        for tenor in TENORS:
            p = node_parameters(expiry, tenor)
            smile = []
            for off in OFFSETS_BP:
                v = hagan_normal_vol(p["alpha"], BETA, p["rho"], p["nu"], SHIFT, p["forward"], p["forward"] + off * 1e-4, p["expiry_years"])
                smile.append({"offset_bp": off, "vol_bp": round(v * 1e4, 2)})
            nodes.append({"expiry": expiry, "tenor": tenor, "forward_pct": round(p["forward"] * 100.0, 4), "atm_bp": round(p["atm"] * 1e4, 2), "smile": smile})
    return {
        "_comment": (
            "SIMPLIFICATION: synthetic USD swaption normal vols (bp) for as-of 2026-09-04, generated by tools/gen_vol_quotes.py from smooth "
            "SABR surfaces (beta 0.5, no shift) with Hagan's normal-vol formula, rounded to 0.01 bp; not vendor data. Each node carries its "
            "forward swap rate (production derives it from the curve). Schema: cpp/market/include/rl/market/vol_quotes.hpp. Node 10Yx30Y is "
            "deliberately stressed (nu 0.70, rho -0.45): it fits as well as every other node and its density goes negative about 170 bp "
            "below the forward (docs/RESULTS.md section 3, notebook 07)."
        ),
        "asOf": AS_OF,
        "currency": "USD",
        "index": "USD-SOFR",
        "unit": "normal_bp",
        "nodes": nodes,
    }


def render(doc: dict) -> str:
    return json.dumps(doc, indent=2) + "\n"


def main() -> None:
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(render(build_document()), encoding="utf-8", newline="\n")
    print(f"wrote {OUT} ({len(EXPIRIES) * len(TENORS)} nodes)")


if __name__ == "__main__":
    main()
