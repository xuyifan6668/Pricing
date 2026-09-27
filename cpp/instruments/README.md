# rl::instruments

Pure trade descriptions: no market data, no model, not templated on the scalar.

## Contract

- `Conventions::usdSofrOis()` etc. are data: change a field, not code, to change a market convention.
  USD SOFR OIS: annual/annual ACT/360, MF on USD-SIFMA, T+2 on Fedwire, 2-day payment delay, single
  period up to 1Y. FF–SOFR basis: quarterly, spread on the FF leg (`// SIMPLIFICATION:` dealers also
  quote annual-pay). SR3 futures: IMM reference period, ACT/360, price = 100 − 100·R.
- `thirdWednesday`, `nextQuarterlyImm`, `immMonthCode` implement the CME quarterly cycle.
- Builders (`OisSwapBuilder`, `BasisSwapBuilder`, `DepositBuilder`, `FraBuilder`, `FutureBuilder`) turn
  conventions + a few numbers into dated trades; `makeFixedLeg` / `makeOvernightLeg` / `makeTermLeg`
  materialise cashflows from a `core::Schedule`. Notional defaults to 1e6; the leg's `Side` carries the sign.
- `cashflow.hpp`: `std::variant<FixedCoupon, OvernightCoupon, TermCoupon, NotionalExchange>` of
  public-member aggregates with stored `tau` and `core::Index` by value. `leg.hpp`: eager vector,
  one side, one currency, ordered by payment date. `swap.hpp`, `basis_swap.hpp`: opposite sides, one
  currency, legs owned. `deposit.hpp`, `fra.hpp`, `future.hpp`: flat records. `swaption.hpp`,
  `bermudan.hpp`: own their underlying `Swap`; validated expiry / exercise dates.

## Invariants (tested)

- IMM dates are Wednesdays between the 15th and 21st; `nextQuarterlyImm` is strictly after its argument.
- The 10y SOFR swap built from 2026-09-08 reproduces the hand-written schedule date for date; tenors ≤ 1Y
  give one period; 18M gives a short front stub plus one annual period.
- A leg's cashflows are ordered by payment date; a swap's legs have opposite sides and one currency.

## Threading

Value types, immutable after construction. Builders are single-use objects; nothing here is shared.
