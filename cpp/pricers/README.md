# pricers

The Bermudan swaption priced three ways under LGM, swap exposure, and risk.

| Piece | File |
|---|---|
| The Bermudan as exercise rights + model; calibration to its own strip | `bermudan.hpp` |
| Crank–Nicolson finite differences, optional Rannacher start-up, exercise boundary | `bermudan_pde.hpp` |
| Longstaff–Schwartz lower bound; Andersen–Broadie upper bound | `bermudan_mc.hpp` |
| Swap exposure profile: EE (closed form and Monte Carlo), PFE, ENE | `exposure.hpp` |
| Par-quote delta by adjoint and by bump; Bermudan vega by calibration instrument | `risk.hpp` |

Everything is in deflated units: under LGM those are martingales in a driftless Gaussian state, so the
PDE is the heat equation in the variance clock ζ, Monte Carlo steps are exact, and no discounting appears
in any pricer.

## Results (docs/RESULTS.md §5–8)

- Lower bound ≤ PDE ≤ upper bound on every trade, duality gap 0.2–0.7 bp of notional.
- PDE: second order, 2 ms at 400 × 400, confirmed to 0.001 bp by exact-transition integration.
- Rannacher start-up: smoother gamma, slightly worse price. Both are measured and tested.
- Switch value rises from 98 to 174 bp as mean reversion goes from −2% to 15%.

## Tested

One exercise date = the European; ≥ every European and ≤ their sum; monotone in exercise rights and in
strike; deep in-the-money limit; mean reversion leaves the Europeans alone and raises the switch value;
convergence order; exercise boundary above the strike for a payer and equal to it at the last date;
in-sample bias; a worse rule gives a lower lower bound and a wider gap; exposure closed form = Monte
Carlo; adjoint delta = bumped delta; vega bounded by each European's own vega.
