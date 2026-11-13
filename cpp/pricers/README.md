# pricers

The Bermudan swaption priced three ways under LGM, swap exposure, and risk.

| Piece | File |
|---|---|
| The Bermudan as exercise rights + model; calibration to its own strip | `bermudan.hpp` |
| Crank–Nicolson finite differences, optional Rannacher start-up, exercise boundary | `bermudan_pde.hpp` |
| Longstaff–Schwartz lower bound; Andersen–Broadie upper bound | `bermudan_mc.hpp` |
| The Bermudan under G2++: exact two-state simulation, regression on both states | `bermudan_g2pp.hpp` |
| CMS spread option under G2++: quadrature, with the two rates' vols and correlation | `cms_spread.hpp` |
| Cheyette by Monte Carlo: Europeans, strip calibration, the Bermudan | `cheyette_mc.hpp` |
| Delta-hedging simulation; sensitivities and P&L explain | `hedging.hpp` |
| The Bermudan with the model calibrated at chosen strikes (`calibrateBermudanAtStrikes`) | `bermudan.hpp` |
| Swap exposure profile: EE (closed form and Monte Carlo), PFE, ENE | `exposure.hpp` |
| Par-quote delta by adjoint and by bump; Bermudan vega by calibration instrument | `risk.hpp` |

Everything is in deflated units: under LGM those are martingales in a driftless Gaussian state, so the
PDE is the heat equation in the variance clock ζ, Monte Carlo steps are exact, and no discounting appears
in any pricer.

## Results (docs/RESULTS.md §5–9, §12)

- Lower bound ≤ PDE ≤ upper bound on every trade, duality gap 0.2–0.7 bp of notional.
- PDE: second order, 2 ms at 400 × 400, confirmed to 0.001 bp by exact-transition integration.
- Rannacher start-up: smoother gamma, slightly worse price. Both are measured and tested.
- Switch value rises from 92 to 173 bp as mean reversion goes from −2% to 15%.
- Under G2++ (§9) factor correlation moves the price by 73 bp with every European repriced; with the
  second factor off the two-factor Monte Carlo reproduces the PDE.
- §12: calibrated at the exercise boundary the Bermudan moves 16 bp; Cheyette with the market's skew moves it
  1.2 bp. Delta hedge noise falls like √interval (278 bp unhedged, 7.4 bp daily); the P&L explain leaves under
  0.6 bp on day-sized moves.
- CMS spread option (§9.4), both swaption vols matched: 10 bp under one factor, 35–70 bp under two.

## Tested

One exercise date = the European; ≥ every European and ≤ their sum; monotone in exercise rights and in
strike; deep in-the-money limit; mean reversion leaves the Europeans alone and raises the switch value;
convergence order; exercise boundary above the strike for a payer and equal to it at the last date;
in-sample bias; a worse rule gives a lower lower bound and a wider gap; exposure closed form = Monte
Carlo; adjoint delta = bumped delta; vega bounded by each European's own vega. Two-factor: one exercise
date = the European; second factor off = the one-factor PDE; ≥ best European and ≤ their sum; monotone in
exercise rights; reproducible.
