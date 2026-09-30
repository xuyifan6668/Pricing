# rl::linear

Leg PV, the overnight compounding kernel, annuity, par rate / par spread, forwards and DV01s off a
`CurveSet<T>` (Andersen & Piterbarg vol. 1 §5–6, ISDA 2021 §7). Called ~1e6 times per
bootstrap solve.

## Contract

- Discounting off `curves.discount()`, projection off `curves.projection(index.name)`.
- A cashflow paying on or before the reference date is worth zero.
- `overnightRate`: realised fixings compound (or average, per `Index::compounding`) day by day up to the
  reference date, the projection curve's DF ratio covers the rest; lookback is applied as an observation
  shift; lockout freezes the last k fixings. Missing fixings throw `InputError` with index, date and coupon.
  `// SIMPLIFICATION:` averaged legs use the continuous average of the projected part (error ≈ r²/720
  relative, far below a tenth of a bp); lockout is not applied to projected days.
- `pv(Leg)` applies `sign()`; `pv(Swap)`/`pv(BasisSwap)` are signed sums of legs.
- `parRate = −sign_fix · PV_float / A_fix`; `parSpread = s − PV / (sign_leg · A_leg)` (uniform spread required).
- `fixedRateDv01 = sign_fix · A · 1e-4`; `parallelDv01` is a one-sided +1bp bump of every zero rate.
- Templates are defined in `detail/linear.ipp` and instantiated for `double` and `rl::aad::Real`; add a
  line to `RL_LINEAR_INSTANTIATE` when you add a function.

## Invariants (tested)

- Fixed coupon PV = N·r·τ·P(pay) (hand example); settled cashflows are zero.
- A spread-free OIS leg paid on accrual end telescopes to `P(start) − P(end)` on any curve shape.
- `pv(swap at parRate) = 0` to 1e-14 of notional; `A·K = −sign·PV_float`.
- Compounded SOFR vs averaged FF on identical curves differ by the compounding basis `≈ r²τ/2` (≈ 1.9 bp
  on a quarter at 3.85%).
- Realised-fixing coupons match a hand-compounded product; lookback/lockout use the documented days.
- `Real` and `double` agree bit for bit and the `Real` result is active.

## Threading

Pure functions of immutable inputs; safe to call from any number of threads.
