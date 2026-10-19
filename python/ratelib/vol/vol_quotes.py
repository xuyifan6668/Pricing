"""Swaption vol screen validation: the checks a vol desk runs before the morning cube fit.

Mirrors the schema of cpp/market/include/rl/market/vol_quotes.hpp. The C++ loader validates structure (a complete
grid, positive vols, strictly increasing offsets); this module validates *market sense* and says what the fit will
struggle with: stale as-of, inconsistent offset sets, ATM term-structure jumps, concave smiles SABR cannot produce,
non-positive strikes that need beta = 0 or a shift.
"""
from __future__ import annotations

import datetime as dt
import json
import re
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

_TENOR_RE = re.compile(r"^(\d+)([DWMY])$", re.IGNORECASE)


class VolQuoteValidationError(ValueError):
    """Raised with every hard error found, one per line."""


def tenor_years(label: str) -> float:
    """Axis year fraction, exactly as market::cubeAxisYears: months/12 for M and Y, days/365 for D and W."""
    m = _TENOR_RE.match(label)
    if not m or int(m.group(1)) <= 0:
        raise ValueError(f"bad tenor label '{label}'")
    n, unit = int(m.group(1)), m.group(2).upper()
    return {"D": n / 365.0, "W": 7.0 * n / 365.0, "M": n / 12.0, "Y": float(n)}[unit]


@dataclass
class SmilePoint:
    offset_bp: float
    vol_bp: float


@dataclass
class VolNode:
    expiry: str
    tenor: str
    forward_pct: float
    atm_bp: float
    smile: list[SmilePoint] = field(default_factory=list)

    @property
    def id(self) -> str:
        return f"{self.expiry}x{self.tenor}"

    @property
    def expiry_years(self) -> float:
        return tenor_years(self.expiry)

    @property
    def tenor_years(self) -> float:
        return tenor_years(self.tenor)

    def points(self) -> list[tuple[float, float]]:
        """(offset_bp, vol_bp) including ATM at offset 0, sorted by offset."""
        return sorted([(0.0, self.atm_bp)] + [(p.offset_bp, p.vol_bp) for p in self.smile])

    def strikes_pct(self) -> list[float]:
        return [self.forward_pct + off / 100.0 for off, _ in self.points()]


@dataclass
class VolQuotes:
    as_of: dt.date
    currency: str
    index: str
    unit: str
    nodes: list[VolNode] = field(default_factory=list)

    def expiries(self) -> list[str]:
        return sorted({n.expiry for n in self.nodes}, key=tenor_years)

    def tenors(self) -> list[str]:
        return sorted({n.tenor for n in self.nodes}, key=tenor_years)

    def find(self, expiry: str, tenor: str) -> VolNode:
        for n in self.nodes:
            if n.expiry == expiry and n.tenor == tenor:
                return n
        raise KeyError(f"{expiry}x{tenor}")

    def atm_matrix(self) -> pd.DataFrame:
        """ATM normal vols in bp, expiries down, tenors across, in axis order."""
        df = pd.DataFrame([(n.expiry, n.tenor, n.atm_bp) for n in self.nodes], columns=["expiry", "tenor", "atm_bp"])
        return df.pivot(index="expiry", columns="tenor", values="atm_bp").reindex(index=self.expiries(), columns=self.tenors())

    def to_frame(self) -> pd.DataFrame:
        """Long format, one row per quoted strike (ATM included)."""
        rows = []
        for n in self.nodes:
            for off, vol in n.points():
                rows.append((n.id, n.expiry, n.tenor, n.expiry_years, n.tenor_years, n.forward_pct, off, n.forward_pct + off / 100.0, vol))
        return pd.DataFrame(rows, columns=["id", "expiry", "tenor", "expiry_years", "tenor_years", "forward_pct", "offset_bp", "strike_pct", "vol_bp"])


def parse_vol_quotes(doc: dict) -> VolQuotes:
    """Build VolQuotes from the JSON document (structure only; see validate_vol_quotes)."""
    for key in ("asOf", "index", "nodes"):
        if key not in doc:
            raise VolQuoteValidationError(f"vol quote file needs '{key}'")
    nodes = []
    for e in doc["nodes"]:
        for key in ("expiry", "tenor", "forward_pct", "atm_bp"):
            if key not in e:
                raise VolQuoteValidationError(f"vol node {e.get('expiry', '?')}x{e.get('tenor', '?')}: missing '{key}'")
        smile = [SmilePoint(float(p["offset_bp"]), float(p["vol_bp"])) for p in e.get("smile", [])]
        nodes.append(VolNode(str(e["expiry"]).upper(), str(e["tenor"]).upper(), float(e["forward_pct"]), float(e["atm_bp"]), smile))
    return VolQuotes(dt.date.fromisoformat(doc["asOf"]), str(doc.get("currency", "USD")), str(doc["index"]), str(doc.get("unit", "normal_bp")), nodes)


def load_vol_quotes(path: str | Path, validate: bool = True) -> VolQuotes:
    with open(path, encoding="utf-8") as f:
        q = parse_vol_quotes(json.load(f))
    if validate:
        validate_vol_quotes(q)
    return q


def validate_vol_quotes(q: VolQuotes, *, today: dt.date | None = None, max_age_days: int = 5, max_atm_jump: float = 0.25) -> list[str]:
    """Return warnings; raise VolQuoteValidationError listing every hard error.

    Hard errors: unit other than normal_bp, bad tenor labels, duplicate nodes, holes in the expiry x tenor grid,
    ATM or smile vols outside (0, 1000) bp, forwards outside (-2%, 25%), offsets that are zero, duplicated or not
    increasing.
    Warnings: stale as-of; nodes quoting different offset sets; an ATM vol that moves by more than `max_atm_jump`
    (relative) between adjacent expiries at the same tenor; a concave smile (both 100 bp wings below ATM: no SABR
    with nu >= 0 fits it); a quoted strike <= 0 (needs beta = 0 or a shift under SABR).
    """
    errors: list[str] = []
    warnings: list[str] = []
    if q.unit != "normal_bp":
        errors.append(f"unit '{q.unit}' is not supported (normal_bp only in v1)")
    seen: set[str] = set()
    for n in q.nodes:
        try:
            tenor_years(n.expiry)
            tenor_years(n.tenor)
        except ValueError as exc:
            errors.append(f"{n.id}: {exc}")
            continue
        if n.id in seen:
            errors.append(f"{n.id}: duplicate node")
        seen.add(n.id)
        if not 0.0 < n.atm_bp < 1000.0:
            errors.append(f"{n.id}: ATM vol {n.atm_bp} bp outside (0, 1000) — wrong unit?")
        if not -2.0 < n.forward_pct < 25.0:
            errors.append(f"{n.id}: forward {n.forward_pct}% outside (-2%, 25%) — wrong unit?")
        offs = [p.offset_bp for p in n.smile]
        if any(o == 0.0 for o in offs):
            errors.append(f"{n.id}: smile offset 0 (the ATM point is atm_bp)")
        if any(b <= a for a, b in zip(offs, offs[1:])):
            errors.append(f"{n.id}: smile offsets not strictly increasing")
        for p in n.smile:
            if not 0.0 < p.vol_bp < 1000.0:
                errors.append(f"{n.id}: smile vol {p.vol_bp} bp at {p.offset_bp:+g} bp outside (0, 1000)")
    if not errors:
        valid = [n for n in q.nodes]
        for e in q.expiries():
            for t in q.tenors():
                if f"{e}x{t}" not in seen:
                    errors.append(f"grid has no node {e}x{t}")
    if errors:
        raise VolQuoteValidationError("\n".join(errors))

    today = today or dt.date.today()
    if (today - q.as_of).days > max_age_days:
        warnings.append(f"as-of {q.as_of} is {(today - q.as_of).days} days old")
    offset_sets = {tuple(p.offset_bp for p in n.smile) for n in valid}
    if len(offset_sets) > 1:
        warnings.append(f"nodes quote {len(offset_sets)} different offset sets")
    exps = q.expiries()
    for t in q.tenors():
        for a, b in zip(exps, exps[1:]):
            va, vb = q.find(a, t).atm_bp, q.find(b, t).atm_bp
            if abs(vb / va - 1.0) > max_atm_jump:
                warnings.append(f"ATM vol jumps {vb / va - 1.0:+.0%} from {a}x{t} to {b}x{t}")
    for n in valid:
        wings = {p.offset_bp: p.vol_bp for p in n.smile}
        lo, hi = min(wings, default=0.0), max(wings, default=0.0)
        if lo < 0.0 < hi and wings[lo] < n.atm_bp and wings[hi] < n.atm_bp:
            warnings.append(f"{n.id}: concave smile (both wings below ATM) — SABR cannot fit it")
        if any(s <= 0.0 for s in n.strikes_pct()):
            warnings.append(f"{n.id}: quoted strike <= 0 — needs beta = 0 or a shift")
    return warnings


def summary(q: VolQuotes) -> str:
    """One paragraph for the notebook."""
    atm = q.atm_matrix()
    return (
        f"as-of {q.as_of} {q.currency} {q.index}: {len(q.nodes)} nodes, {len(q.expiries())} expiries x {len(q.tenors())} tenors, "
        f"ATM {atm.min().min():.1f}-{atm.max().max():.1f} bp, offsets {sorted({p.offset_bp for n in q.nodes for p in n.smile})}"
    )


def default_path() -> Path:
    from ratelib.demos import repo_root

    return repo_root() / "data" / "vol" / "usd_swaption_cube_2026-09-04.json"
