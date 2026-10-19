"""Cube assembly orchestration: policy -> per-node C++ fits -> cube JSON -> report.

What lives here and why: the SABR policy (beta and the shift per currency), the choice of
which nodes to fit and in what order, where the cube is written, and the report. All of that changes during
research; the fit kernel itself is C++ (`vol::fitSabrSmile`, driven through `rl_demo_vol` until the pybind11 module
arrives ).

DESIGN: the policy is a frozen dataclass passed explicitly to every call. A module-level "current beta" would make
two notebooks that ran the same code disagree, and the snapshot hash would not see why.
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from ratelib.demos import run_demo, run_demo_text
from ratelib.vol.report import CubeReport, load_arbitrage_csv, load_fit_csv


@dataclass(frozen=True)
class SabrPolicy:
    """Desk policy for a currency's SABR cube: never calibrated, never differentiated."""

    beta: float = 0.5
    shift: float = 0.0

    def args(self) -> list[str]:
        return [repr(self.beta), repr(self.shift)]


# USD, v1: beta = 0.5 and no shift. Strikes on the screen are all positive; the arbitrage scan shows what that costs
# at low strikes on long expiries (notebook 07).
USD_POLICY = SabrPolicy(beta=0.5, shift=0.0)


def fit_nodes(policy: SabrPolicy = USD_POLICY, build_dir: Path | None = None) -> pd.DataFrame:
    """One row per node: fitted parameters, rms/max residual in bp, iterations, status, flags."""
    return run_demo("rl_demo_vol", "fit", "all", *policy.args(), build_dir=build_dir)


def fit_residuals(policy: SabrPolicy = USD_POLICY, build_dir: Path | None = None) -> pd.DataFrame:
    """One row per quote: market vol, model vol, residual in bp."""
    return run_demo("rl_demo_vol", "residuals", *policy.args(), build_dir=build_dir)


def arbitrage_scan(policy: SabrPolicy = USD_POLICY, n_std_dev: float = 4.0, points: int = 801, build_dir: Path | None = None) -> pd.DataFrame:
    """One row per node: the first strike below/above the forward with a negative density, over F ± n_std_dev·σ√T."""
    text = run_demo_text("rl_demo_vol", "arbitrage", *policy.args(), n_std_dev, points, build_dir=build_dir)
    return load_arbitrage_csv(text)


def build_cube(policy: SabrPolicy = USD_POLICY, out_path: Path | None = None, build_dir: Path | None = None) -> dict:
    """Fit every node and return the cube document (market::toJson schema); write it if `out_path` is given."""
    text = run_demo_text("rl_demo_vol", "cube", *policy.args(), build_dir=build_dir)
    doc = json.loads(text)
    if out_path is not None:
        Path(out_path).write_text(json.dumps(doc, indent=2) + "\n", encoding="utf-8")
    return doc


def query_smile(cube_path: Path, expiry_years: float, tenor_years: float, forward: float, k_min: float, k_max: float, n: int = 101,
                build_dir: Path | None = None) -> pd.DataFrame:
    """Smile at an arbitrary (expiry, tenor) from a fitted cube file, at the caller's forward (decimal inputs)."""
    return run_demo("rl_demo_vol", "query", cube_path, expiry_years, tenor_years, forward, k_min, k_max, n, build_dir=build_dir)


@dataclass
class CubeBuild:
    policy: SabrPolicy
    report: CubeReport
    cube: dict
    cube_path: Path | None


def build(policy: SabrPolicy = USD_POLICY, out_dir: Path | None = None, build_dir: Path | None = None, tolerance_bp: float = 0.1) -> CubeBuild:
    """The morning cube build: fit, assemble, scan for arbitrage, report. Writes cube JSON and the report if `out_dir`."""
    fits = fit_nodes(policy, build_dir)
    report = load_fit_csv(fits.to_csv(index=False), tolerance_bp=tolerance_bp, residuals=fit_residuals(policy, build_dir),
                          arbitrage=arbitrage_scan(policy, build_dir=build_dir))
    cube_path = None
    if out_dir is not None:
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        cube_path = out_dir / "usd_swaption_cube_fitted.json"
        report.write(out_dir / "cube_report.md")
    cube = build_cube(policy, cube_path, build_dir)
    return CubeBuild(policy=policy, report=report, cube=cube, cube_path=cube_path)
