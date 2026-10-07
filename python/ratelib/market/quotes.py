"""Quote-file validation: the checks a desk runs before a curve build, in the language they get rewritten in
weekly (docs/DESIGN.md §3: "market data assembly / validation: Python").

Mirrors the schema of cpp/market/include/rl/market/quote.hpp. The C++ loader validates structure; this
module validates *market sense*: monotone maturities per type, plausible levels, a consistent IMM strip,
stale as-of dates, and the curve assignment that makes the global solve necessary.
"""
from __future__ import annotations

import datetime as dt
import json
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

QUOTE_TYPES = {"deposit", "fra", "future", "ois_swap", "basis_swap", "turn"}
UNITS = {"decimal", "pct", "bp", "price"}
_TENOR_RE = re.compile(r"^(\d+)([DWMY])$", re.IGNORECASE)
_TENOR_DAYS = {"D": 1.0, "W": 7.0, "M": 30.4375, "Y": 365.25}


class QuoteValidationError(ValueError):
    """Raised with every problem found, one per line, so the whole file is reviewed at once."""


@dataclass
class Quote:
    id: str
    type: str
    index: str
    value: float
    unit: str = "pct"
    reference_index: str = ""
    curve: str = ""
    tenor: str | None = None
    start: dt.date | None = None
    end: dt.date | None = None
    convexity_bp: float = 0.0

    @property
    def rate(self) -> float:
        """Decimal rate (raises for price quotes)."""
        if self.unit == "price":
            raise ValueError(f"{self.id}: price quote has no rate")
        return self.value / {"decimal": 1.0, "pct": 100.0, "bp": 1e4}[self.unit]

    @property
    def curve_name(self) -> str:
        return self.curve or self.index

    def tenor_days(self) -> float | None:
        if self.tenor is None:
            return None
        m = _TENOR_RE.match(self.tenor)
        if not m:
            return None
        return int(m.group(1)) * _TENOR_DAYS[m.group(2).upper()]

    def maturity_days(self, as_of: dt.date) -> float | None:
        """Approximate days from as-of to the pillar, for ordering checks."""
        if self.end is not None:
            return float((self.end - as_of).days)
        return self.tenor_days()


@dataclass
class QuoteSet:
    as_of: dt.date
    currency: str
    quotes: list[Quote] = field(default_factory=list)

    def of_type(self, t: str) -> list[Quote]:
        return [q for q in self.quotes if q.type == t]

    def by_curve(self) -> dict[str, list[Quote]]:
        out: dict[str, list[Quote]] = {}
        for q in self.quotes:
            out.setdefault(q.curve_name, []).append(q)
        return out

    def find(self, quote_id: str) -> Quote:
        for q in self.quotes:
            if q.id == quote_id:
                return q
        raise KeyError(quote_id)


def _parse_date(s: str | None) -> dt.date | None:
    return None if s is None else dt.date.fromisoformat(s)


def parse_quotes(doc: dict) -> QuoteSet:
    """Build a QuoteSet from the JSON document (no market-sense validation; see validate_quotes)."""
    if "asOf" not in doc or "quotes" not in doc:
        raise QuoteValidationError("quote file needs 'asOf' and 'quotes'")
    quotes = []
    for e in doc["quotes"]:
        for key in ("id", "type", "value"):
            if key not in e:
                raise QuoteValidationError(f"quote {e.get('id', '?')}: missing '{key}'")
        quotes.append(
            Quote(
                id=str(e["id"]),
                type=str(e["type"]).lower(),
                index=str(e.get("index", "")),
                value=float(e["value"]),
                unit=str(e.get("unit", "pct")).lower(),
                reference_index=str(e.get("reference_index", "")),
                curve=str(e.get("curve", "")),
                tenor=e.get("tenor"),
                start=_parse_date(e.get("start")),
                end=_parse_date(e.get("end")),
                convexity_bp=float(e.get("convexity_bp", 0.0)),
            )
        )
    return QuoteSet(as_of=dt.date.fromisoformat(doc["asOf"]), currency=str(doc.get("currency", "USD")), quotes=quotes)


def load_quotes(path: str | Path, validate: bool = True) -> QuoteSet:
    with open(path, encoding="utf-8") as f:
        qs = parse_quotes(json.load(f))
    if validate:
        validate_quotes(qs)
    return qs


def _third_wednesday(year: int, month: int) -> dt.date:
    d = dt.date(year, month, 1)
    return d + dt.timedelta(days=(2 - d.weekday()) % 7 + 14)


def validate_quotes(qs: QuoteSet, *, today: dt.date | None = None, max_age_days: int = 5) -> list[str]:
    """Return the list of warnings; raise QuoteValidationError on hard errors.

    Hard errors: duplicate ids, unknown type/unit, missing fields per type, rates outside (-2%, 25%),
    futures prices outside (90, 100.5), basis spreads beyond ±100 bp, non-increasing maturities within a
    (curve, type) group, a futures strip that is not consecutive quarterly IMM periods, two instruments
    pinning the same (curve, maturity).
    Warnings: stale as-of (older than max_age_days relative to `today`), a discount curve with fewer than
    five pillars, a projection curve pinned by fewer than three instruments, no turn quote for a curve that
    spans a year end, a basis swap assigned to the reference curve (the cross-assignment that couples
    curves) with no comment.
    """
    errors: list[str] = []
    warnings: list[str] = []
    seen: set[str] = set()
    for q in qs.quotes:
        if q.id in seen:
            errors.append(f"{q.id}: duplicate id")
        seen.add(q.id)
        if q.type not in QUOTE_TYPES:
            errors.append(f"{q.id}: unknown type '{q.type}'")
            continue
        if q.unit not in UNITS:
            errors.append(f"{q.id}: unknown unit '{q.unit}'")
            continue
        if not q.index:
            errors.append(f"{q.id}: missing index")
        if q.type in ("deposit", "ois_swap", "basis_swap") and q.tenor is None and (q.start is None or q.end is None):
            errors.append(f"{q.id}: needs a tenor or start/end")
        if q.type in ("fra", "future", "turn") and (q.start is None or q.end is None):
            errors.append(f"{q.id}: needs start and end")
        if q.start is not None and q.end is not None and q.end <= q.start:
            errors.append(f"{q.id}: end must follow start")
        if q.tenor is not None and not _TENOR_RE.match(q.tenor):
            errors.append(f"{q.id}: bad tenor '{q.tenor}'")
        if q.type == "basis_swap":
            if not q.reference_index:
                errors.append(f"{q.id}: basis swap needs reference_index")
            elif q.reference_index == q.index:
                errors.append(f"{q.id}: basis swap legs on the same index")
        if q.type == "future":
            if q.unit != "price":
                errors.append(f"{q.id}: futures are quoted as prices")
            elif not 90.0 < q.value < 100.5:
                errors.append(f"{q.id}: futures price {q.value} out of range")
            if q.start is not None and q.start != _third_wednesday(q.start.year, q.start.month):
                errors.append(f"{q.id}: start {q.start} is not an IMM date")
        elif q.unit == "price":
            errors.append(f"{q.id}: only futures are quoted as prices")
        else:
            r = q.rate
            if q.type == "basis_swap" or q.type == "turn":
                if abs(r) > 0.01:
                    errors.append(f"{q.id}: spread/jump {r * 1e4:.1f} bp beyond ±100 bp — wrong unit?")
            elif not -0.02 < r < 0.25:
                errors.append(f"{q.id}: rate {r:.4%} outside (-2%, 25%) — wrong unit?")

    # Maturity ordering within (curve, type), and unique pillars per curve.
    groups: dict[tuple[str, str], list[Quote]] = {}
    for q in qs.quotes:
        if q.type == "turn":
            continue
        groups.setdefault((q.curve_name, q.type), []).append(q)
    for (curve, t), qlist in groups.items():
        days = [q.maturity_days(qs.as_of) for q in qlist]
        for a, b, qa, qb in zip(days, days[1:], qlist, qlist[1:]):
            if a is not None and b is not None and b <= a:
                errors.append(f"{qb.id}: maturity not after {qa.id} within {curve}/{t}")
    for curve, qlist in qs.by_curve().items():
        pillars: dict[float, str] = {}
        for q in qlist:
            if q.type == "turn":
                continue
            d = q.maturity_days(qs.as_of)
            if d is None:
                continue
            if d in pillars:
                errors.append(f"{q.id} and {pillars[d]} pin the same pillar on {curve}")
            pillars[d] = q.id
    # Futures strip: consecutive quarterly IMM periods.
    futs = sorted(qs.of_type("future"), key=lambda q: q.start or qs.as_of)
    for a, b in zip(futs, futs[1:]):
        if a.end != b.start:
            errors.append(f"{b.id}: futures strip gap — {a.id} ends {a.end}, {b.id} starts {b.start}")

    if errors:
        raise QuoteValidationError("\n".join(errors))

    # Warnings.
    today = today or dt.date.today()
    if (today - qs.as_of).days > max_age_days:
        warnings.append(f"as-of {qs.as_of} is {(today - qs.as_of).days} days old")
    for curve, qlist in qs.by_curve().items():
        n = sum(1 for q in qlist if q.type != "turn")
        if n < 3:
            warnings.append(f"curve {curve} is pinned by only {n} instrument(s)")
    for q in qs.quotes:
        if q.type == "basis_swap" and q.curve and q.curve == q.reference_index:
            warnings.append(f"{q.id}: assigned to the reference curve {q.curve} (cross-assignment; this couples the curves)")
    has_turn = any(q.type == "turn" for q in qs.quotes)
    if not has_turn:
        warnings.append("no turn-of-year quote")
    return warnings


def summary(qs: QuoteSet) -> str:
    """One line per (curve, type): count and maturity range. For the notebook."""
    lines = [f"as-of {qs.as_of} {qs.currency}: {len(qs.quotes)} quotes"]
    for curve, qlist in qs.by_curve().items():
        by_type: dict[str, list[Quote]] = {}
        for q in qlist:
            by_type.setdefault(q.type, []).append(q)
        for t, ql in by_type.items():
            lines.append(f"  {curve:14s} {t:10s} {len(ql):3d}  {ql[0].tenor or ql[0].start} … {ql[-1].tenor or ql[-1].end}")
    return "\n".join(lines)


def iter_ids(quotes: Iterable[Quote]) -> list[str]:
    return [q.id for q in quotes]
