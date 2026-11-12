#!/usr/bin/env python3
"""Generate the synthetic USD swaption vol screen data/vol/usd_swaption_cube_2026-09-04.json.

Run:  uv run tools/gen_vol_quotes.py          (byte-identical output if the surface is unchanged)
Test: uv run pytest tools/

SIMPLIFICATION: synthetic marks, not vendor data. The screen is generated in a plausible September-2026 USD
configuration: smooth forward, rho and nu surfaces (forwards 3.3-3.8%, negative skew growing with expiry, vol-of-vol
decaying with expiry), beta = 0.5, no shift, and ATM normal vols (63-113 bp) implied by a two-factor Gaussian model
with the parameters below. The ATM matrix therefore has the tenor structure of a model in which rates of different
maturities are imperfectly correlated, and a two-factor calibration to it has something to find. Quotes are Hagan normal vols
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


# The ATM surface comes from a two-factor Gaussian model, so that the screen carries a correlation structure for a
# two-factor calibration to find. Brigo-Mercurio G2++ notation: dx = -a x dt + sigma dW1, dy = -b y dt + eta dW2,
# d<W1, W2> = rho dt, with both volatilities scaled by s(u) = 1 + c exp(-u / lam) (higher volatility in the near term).
G2_A, G2_B = 0.40, 0.03          # a fast and a slow factor
G2_SIGMA, G2_ETA = 0.0123, 0.0097
G2_RHO = -0.70                   # negatively correlated: the curve steepens and flattens as well as shifting
G2_SCALE_C, G2_SCALE_LAMBDA = 0.46, 1.3


def _bond_loading(k: float, u: float) -> float:
    return u if abs(k * u) < 1e-9 else -math.expm1(-k * u) / k


def _decay_integral(k: float, m: float, e: float) -> float:
    """Integral over u in [0, e] of exp(-m u) exp(-k (e - u))."""
    return e * math.exp(-k * e) if abs(k - m) < 1e-9 else (math.exp(-m * e) - math.exp(-k * e)) / (k - m)


def atm_vol(e: float, t: float) -> float:
    """Normal vol of the forward swap rate in the two-factor model, with the swap's bond weights frozen at today's
    values on a flat curve at the node's forward (annual fixed coupons). Frozen weights are the standard market
    approximation; an exact two-factor pricer reproduces these vols to a few tenths of a basis point.

    The swap rate's loading on factor i at time u < e is sigma_i s(u) exp(-k_i (e - u)) G_i, where
    G_i = [P_n B_i(n) + S sum_j P_j B_i(j)] / A collects the bonds of the swap; its variance to expiry is the
    integral of the squared, correlated sum of the two loadings."""
    f = forward(e, t)
    n = int(round(t))
    bonds = [(1.0 + f) ** (-j) for j in range(n + 1)]
    annuity = sum(bonds[1:])
    speeds, vols = (G2_A, G2_B), (G2_SIGMA, G2_ETA)
    g = [(bonds[n] * _bond_loading(k, n) + f * sum(bonds[j] * _bond_loading(k, j) for j in range(1, n + 1))) / annuity for k in speeds]
    c, lam = G2_SCALE_C, G2_SCALE_LAMBDA
    variance = 0.0
    for i in range(2):
        for j in range(2):
            kk = speeds[i] + speeds[j]
            scale = _decay_integral(kk, 0.0, e) + 2.0 * c * _decay_integral(kk, 1.0 / lam, e) + c * c * _decay_integral(kk, 2.0 / lam, e)
            variance += (1.0 if i == j else G2_RHO) * vols[i] * vols[j] * g[i] * g[j] * scale
    return math.sqrt(variance / e)


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
            "SIMPLIFICATION: synthetic USD swaption normal vols (bp) for as-of 2026-09-04, generated by tools/gen_vol_quotes.py: ATM vols from a "
            "two-factor Gaussian model (a 0.40, b 0.03, rho -0.70), smiles from smooth SABR rho and nu surfaces (beta 0.5, no shift) with "
            "Hagan's normal-vol formula, rounded to 0.01 bp; not vendor data. Each node carries its "
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
