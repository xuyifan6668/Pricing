#!/usr/bin/env python3
"""Reprice the benchmark case of `rl_analysis ... benchmark` with QuantLib and compare.

    uv run --with QuantLib python tools/benchmark_quantlib.py [--results build/results]

Every other check in this repository compares the library with itself: closed form against quadrature,
PDE against Monte Carlo. Those share one author, one curve stack and one reading of the conventions, so
a mistake common to all of them would pass. This script is the outside view. The C++ side hands over only
what defines the case: discount factors on a weekly grid, the swap's accrual dates, the Hull-White
parameters. QuantLib then builds its own swap, computes its own accruals, forwards and annuity, and
prices with its own engines:

    swap value, par rate      DiscountingSwapEngine
    European swaptions        JamshidianSwaptionEngine (Hull-White closed form)
    Bermudan swaption         TreeSwaptionEngine (trinomial lattice) and FdHullWhiteSwaptionEngine

What is shared is the model (Hull-White with the same two numbers), the curve and the dates. What is
not shared is every line of code. The floating leg is built in QuantLib as a term-rate leg on the same
curve; a compounded overnight coupon and a term coupon over the same period have the same value on one
curve, so the contract is the same.

Writes <results>/benchmark_quantlib.csv. QuantLib is not a dependency of the package: it is needed only
to regenerate this comparison.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def _read(path: Path) -> list[dict[str, str]]:
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def compare(results: Path, tree_steps: int = 2000, fd_time: int = 1000, fd_space: int = 1000) -> list[dict[str, object]]:
    """One row per compared quantity: case, quantity, this library, QuantLib, difference in bp of notional."""
    import QuantLib as ql

    def date(iso: str) -> "ql.Date":
        y, m, d = (int(x) for x in iso.split("-"))
        return ql.Date(d, m, y)

    curve_rows = _read(results / "benchmark_curve.csv")
    today = date(curve_rows[0]["date"])
    ql.Settings.instance().evaluationDate = today
    curve = ql.DiscountCurve([date(r["date"]) for r in curve_rows], [float(r["discount_factor"]) for r in curve_rows], ql.Actual365Fixed())
    curve.enableExtrapolation()
    handle = ql.YieldTermStructureHandle(curve)
    calendar, act360 = ql.NullCalendar(), ql.Actual360()
    # A 12-month term index on the discount curve itself, with no fixing lag and no date adjustment.
    index = ql.IborIndex("bench", ql.Period(12, ql.Months), 0, ql.USDCurrency(), calendar, ql.Unadjusted, False, act360, handle)

    legs = _read(results / "benchmark_legs.csv")
    europeans = _read(results / "benchmark_europeans.csv")
    out: list[dict[str, object]] = []

    def row(case: str, quantity: str, ours: float, theirs: float, scale: float = 1e4) -> None:
        out.append({"case": case, "quantity": quantity, "ratelib": ours * scale, "quantlib": theirs * scale, "difference": (ours - theirs) * scale})

    for case in _read(results / "benchmark_cases.csv"):
        cid, payer = case["case"], case["payer"] == "1"
        rate, kappa, sigma = float(case["fixed_rate"]), float(case["kappa"]), float(case["sigma"])
        fixed = [r for r in legs if r["case"] == cid and r["leg"] == "fixed"]
        floating = [r for r in legs if r["case"] == cid and r["leg"] == "float"]

        def schedule(rows: list[dict[str, str]], first: int = 0) -> "ql.Schedule":
            dates = [date(rows[first]["accrual_start"])] + [date(r["accrual_end"]) for r in rows[first:]]
            return ql.Schedule(dates, calendar, ql.Unadjusted)

        def swap(first: int = 0) -> "ql.VanillaSwap":
            kind = ql.VanillaSwap.Payer if payer else ql.VanillaSwap.Receiver
            s = ql.VanillaSwap(kind, 1.0, schedule(fixed, first), rate, act360, schedule(floating, first), index, 0.0, act360)
            s.setPricingEngine(ql.DiscountingSwapEngine(handle))
            return s

        label = f"{'payer' if payer else 'receiver'} {rate:.2%}"
        whole = swap()
        row(label, "swap value", float(case["swap_npv"]), whole.NPV())
        row(label, "par rate", float(case["par_rate"]), whole.fairRate())

        model = ql.HullWhite(handle, kappa, sigma)
        exercise_dates = []
        for e in (r for r in europeans if r["case"] == cid):
            first = next(i for i, r in enumerate(fixed) if r["accrual_start"] == e["exercise_date"])
            exercise_dates.append(date(e["exercise_date"]))
            option = ql.Swaption(swap(first), ql.EuropeanExercise(exercise_dates[-1]))
            option.setPricingEngine(ql.JamshidianSwaptionEngine(model))
            row(label, f"European, exercise {e['exercise_date']}", float(e["european"]), option.NPV())

        bermudan = ql.Swaption(whole, ql.BermudanExercise(exercise_dates))
        ours = float(case["bermudan_pde"])
        bermudan.setPricingEngine(ql.TreeSwaptionEngine(model, tree_steps))
        row(label, f"Bermudan, PDE against QuantLib lattice ({tree_steps} steps)", ours, bermudan.NPV())
        bermudan.setPricingEngine(ql.FdHullWhiteSwaptionEngine(model, fd_time, fd_space))
        row(label, f"Bermudan, PDE against QuantLib finite differences ({fd_time}x{fd_space})", ours, bermudan.NPV())
    return out


def write(results: Path, **kwargs: int) -> list[dict[str, object]]:
    rows = compare(results, **kwargs)
    with open(results / "benchmark_quantlib.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["case", "quantity", "ratelib", "quantlib", "difference"])
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.8f}" if isinstance(v, float) else v) for k, v in r.items()})
    return rows


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", default=str(ROOT / "build" / "results"))
    args = ap.parse_args()
    rows = write(Path(args.results))
    worst: dict[str, float] = {}
    for r in rows:
        kind = str(r["quantity"]).split(",")[0]
        worst[kind] = max(worst.get(kind, 0.0), abs(float(r["difference"])))  # type: ignore[arg-type]
    for kind, w in worst.items():
        print(f"{kind:12s} worst difference {w:.2e} (bp of notional; par rate in bp of rate)")


if __name__ == "__main__":
    main()
