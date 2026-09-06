# rl::core

Dates, calendars, day counts, schedules, currencies and index definitions. No market data, no
numerics, no templates on the scalar type: everything here is passive and is never differentiated.
Implements: (Hinnant date algorithms, ISDA 2006 day counts, anchor-based schedule
generation with stub detection); rationale in `// DESIGN:` comments.

## Contract

- `Date` is a serial day with epoch 1899-12-30 (Excel's), valid for 0001-01-01 … 9999-12-31.
  Construction validates the civil triple and throws `rl::InputError`.
- `Calendar` is a name plus a sorted, de-duplicated holiday table loaded from
  `data/holidays/<file>.txt`. Weekends are Saturday and Sunday for every calendar in v1
  (`// SIMPLIFICATION:` Middle-East weekends are not modelled). Holiday tables cover 2000–2075;
  outside that range a calendar silently degrades to weekends-only.
- `DayCount::yearFraction(dc, start, end)` is antisymmetric: `τ(a, b) = −τ(b, a)`.
- `Schedule::generate` returns strictly increasing unadjusted dates with `front() == effective`,
  `back() == termination`; adjusted dates are non-decreasing; `payment ≥ end`.
- `Index` is a value type compared by name. It describes conventions only.

## Invariants (all tested in `test/`)

- `civilFromDays(daysFromCivil(y, m, d)) == {y, m, d}` for every valid triple.
- `weekday(serial + 7k) == weekday(serial)`.
- `addMonths` clamps and is applied from the anchor, never iteratively.
- `adjust(d, bdc) == d` when `d` is a business day, for every convention.
- ACT/ACT ISDA telescopes across arbitrary sub-periods.

## Threading

Immutable after construction; `Calendar::get` caches behind a mutex and returns shared handles.
Everything else is a value type. Safe to share across the parallel risk driver.

## Data

`tools/gen_holidays.py` regenerates `data/holidays/*.txt`; `tools/test_gen_holidays.py` pins the
generator to independently known dates. The C++ loader never computes a holiday.
