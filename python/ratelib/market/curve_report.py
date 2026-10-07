"""Curve repricing report: the first table of the calibration report.

Consumes the CSV printed by `rl_demo_bootstrap reprice` (or `RepricingReport::toCsv()` from C++):
id,type,curve,maturity,quote,model,residual,residual_bp,raw_residual — and renders it as a DataFrame,
a Markdown table and a one-line verdict. Nothing here prices anything.
"""
from __future__ import annotations

import io
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

REQUIRED_COLUMNS = ["id", "type", "curve", "maturity", "quote", "model", "residual", "residual_bp"]


@dataclass
class CurveReport:
    table: pd.DataFrame
    tolerance_bp: float = 1e-8

    @property
    def max_abs_bp(self) -> float:
        return float(self.table["residual_bp"].abs().max()) if len(self.table) else 0.0

    @property
    def mean_abs_bp(self) -> float:
        return float(self.table["residual_bp"].abs().mean()) if len(self.table) else 0.0

    @property
    def passed(self) -> bool:
        return self.max_abs_bp < self.tolerance_bp

    def worst(self, n: int = 5) -> pd.DataFrame:
        return self.table.reindex(self.table["residual_bp"].abs().sort_values(ascending=False).index).head(n)

    def by_curve(self) -> pd.DataFrame:
        g = self.table.groupby("curve")["residual_bp"]
        return pd.DataFrame({"instruments": g.size(), "max_abs_bp": g.apply(lambda s: s.abs().max()), "mean_abs_bp": g.apply(lambda s: s.abs().mean())})

    def verdict(self) -> str:
        status = "PASS" if self.passed else "FAIL"
        return f"{status}: {len(self.table)} instruments, max |residual| {self.max_abs_bp:.3e} bp, mean {self.mean_abs_bp:.3e} bp (tolerance {self.tolerance_bp:g} bp)"

    def to_markdown(self, float_fmt: str = "%.6g") -> str:
        cols = ["id", "type", "curve", "maturity", "quote", "model", "residual_bp"]
        t = self.table[cols].copy()
        t["quote"] = t["quote"].map(lambda v: float_fmt % v)
        t["model"] = t["model"].map(lambda v: float_fmt % v)
        t["residual_bp"] = t["residual_bp"].map(lambda v: "%.2e" % v)
        header = "| " + " | ".join(cols) + " |\n|" + "|".join("---" for _ in cols) + "|\n"
        rows = "".join("| " + " | ".join(str(v) for v in r) + " |\n" for r in t.itertuples(index=False))
        return header + rows + f"\n{self.verdict()}\n"

    def write(self, path: str | Path) -> None:
        p = Path(path)
        if p.suffix == ".md":
            p.write_text(self.to_markdown(), encoding="utf-8")
        else:
            self.table.to_csv(p, index=False)


def load_reprice_csv(text_or_path: str | Path, tolerance_bp: float = 1e-8) -> CurveReport:
    """Parse the reprice CSV from a path or from the CSV text itself."""
    src = text_or_path
    if isinstance(src, Path) or (isinstance(src, str) and "\n" not in src and Path(src).exists()):
        df = pd.read_csv(src)
    else:
        df = pd.read_csv(io.StringIO(str(src)))
    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"reprice CSV is missing columns: {missing}")
    df["maturity"] = pd.to_datetime(df["maturity"]).dt.date
    return CurveReport(table=df, tolerance_bp=tolerance_bp)


def from_dataframe(df: pd.DataFrame, tolerance_bp: float = 1e-8) -> CurveReport:
    return load_reprice_csv(df.to_csv(index=False), tolerance_bp)
