"""Cube calibration report: the second table of the calibration report.

Consumes the CSV printed by `rl_demo_vol fit all` (one row per node), optionally `rl_demo_vol residuals` (one row
per quote) and `rl_demo_vol arbitrage` (one row per node), and renders: the residual table, the worst nodes, an
expiry x tenor heat-map of max residuals, the negative-density map, and a one-line verdict. Nothing here fits.
"""
from __future__ import annotations

import io
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from ratelib.vol.vol_quotes import tenor_years

FIT_COLUMNS = ["id", "expiry", "tenor", "forward_pct", "atm_bp", "alpha", "beta", "rho", "nu", "rms_bp", "max_bp", "iterations", "status", "ok", "flags"]
ARBITRAGE_COLUMNS = ["id", "expiry", "tenor", "first_negative_below_pct", "first_negative_above_pct", "call_spread_ok", "nonpositive_vol", "arbitrage_free"]


def _read(text_or_path: str | Path) -> pd.DataFrame:
    src = text_or_path
    if isinstance(src, Path) or (isinstance(src, str) and "\n" not in src and Path(src).exists()):
        return pd.read_csv(src)
    return pd.read_csv(io.StringIO(str(src)))


def _axis_pivot(df: pd.DataFrame, value: str) -> pd.DataFrame:
    exps = sorted(df["expiry"].unique(), key=tenor_years)
    tens = sorted(df["tenor"].unique(), key=tenor_years)
    return df.pivot(index="expiry", columns="tenor", values=value).reindex(index=exps, columns=tens)


@dataclass
class CubeReport:
    fits: pd.DataFrame
    residuals: pd.DataFrame | None = None
    arbitrage: pd.DataFrame | None = None
    tolerance_bp: float = 0.1

    @property
    def max_abs_bp(self) -> float:
        return float(self.fits["max_bp"].abs().max()) if len(self.fits) else 0.0

    @property
    def mean_rms_bp(self) -> float:
        return float(self.fits["rms_bp"].mean()) if len(self.fits) else 0.0

    def failed_nodes(self) -> pd.DataFrame:
        bad = (self.fits["ok"] != 1) | (self.fits["flags"] != "") | (self.fits["max_bp"].abs() >= self.tolerance_bp)
        return self.fits[bad]

    @property
    def passed(self) -> bool:
        return len(self.fits) > 0 and self.failed_nodes().empty

    def worst(self, n: int = 5) -> pd.DataFrame:
        return self.fits.reindex(self.fits["max_bp"].abs().sort_values(ascending=False).index).head(n)

    def pivot(self, value: str = "max_bp") -> pd.DataFrame:
        """Expiry x tenor matrix of a per-node column, axes in year order."""
        return _axis_pivot(self.fits, value)

    def arbitrage_nodes(self) -> list[str]:
        if self.arbitrage is None:
            return []
        return list(self.arbitrage.loc[self.arbitrage["arbitrage_free"] != 1, "id"])

    def negative_density_map(self) -> pd.DataFrame:
        """Expiry x tenor matrix of the first strike (pct) below the forward where the density goes negative; NaN = clean."""
        if self.arbitrage is None:
            raise ValueError("no arbitrage scan attached (run `rl_demo_vol arbitrage` and pass it in)")
        return _axis_pivot(self.arbitrage, "first_negative_below_pct")

    def verdict(self) -> str:
        status = "PASS" if self.passed else "FAIL"
        text = (
            f"{status}: {len(self.fits)} nodes, max |residual| {self.max_abs_bp:.3e} bp, mean rms {self.mean_rms_bp:.3e} bp "
            f"(tolerance {self.tolerance_bp:g} bp), {len(self.failed_nodes())} flagged"
        )
        if self.arbitrage is not None:
            text += f"; {len(self.arbitrage_nodes())} nodes with a negative density inside the scanned range"
        return text

    def to_markdown(self) -> str:
        cols = ["id", "forward_pct", "atm_bp", "alpha", "rho", "nu", "rms_bp", "max_bp", "iterations", "status", "flags"]
        t = self.fits[cols].copy()
        for c, fmt in [("forward_pct", "%.4f"), ("atm_bp", "%.2f"), ("alpha", "%.5f"), ("rho", "%.4f"), ("nu", "%.4f"), ("rms_bp", "%.2e"), ("max_bp", "%.2e")]:
            t[c] = t[c].map(lambda v, f=fmt: f % v)
        header = "| " + " | ".join(cols) + " |\n|" + "|".join("---" for _ in cols) + "|\n"
        rows = "".join("| " + " | ".join(str(v) for v in r) + " |\n" for r in t.itertuples(index=False))
        return header + rows + f"\n{self.verdict()}\n"

    def write(self, path: str | Path) -> None:
        p = Path(path)
        if p.suffix == ".md":
            p.write_text(self.to_markdown(), encoding="utf-8")
        else:
            self.fits.to_csv(p, index=False)


def load_fit_csv(text_or_path: str | Path, tolerance_bp: float = 0.1, residuals: pd.DataFrame | None = None, arbitrage: pd.DataFrame | None = None) -> CubeReport:
    """Parse `rl_demo_vol fit all` output from a path or the CSV text itself."""
    df = _read(text_or_path)
    missing = [c for c in FIT_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"fit CSV is missing columns: {missing}")
    df["flags"] = df["flags"].fillna("").astype(str)
    return CubeReport(fits=df, residuals=residuals, arbitrage=arbitrage, tolerance_bp=tolerance_bp)


def load_arbitrage_csv(text_or_path: str | Path) -> pd.DataFrame:
    """Parse `rl_demo_vol arbitrage` output; empty strike cells become NaN."""
    df = _read(text_or_path)
    missing = [c for c in ARBITRAGE_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"arbitrage CSV is missing columns: {missing}")
    return df
