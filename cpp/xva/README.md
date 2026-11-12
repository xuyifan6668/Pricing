# xva

Counterparty credit and funding adjustments on top of the pricing models.

| Piece | File |
|---|---|
| Credit curve: hazard rates bootstrapped from CDS spreads | `credit.hpp` |
| Exposure engine: a netting set of swaps on simulated paths, with collateral; one- and two-factor | `exposure_engine.hpp` |
| CVA, DVA, FCA, FBA from exposure; one swap's exposure in closed form | `xva.hpp` |
| Exposure of a Bermudan swaption: an option until exercised, a swap afterwards; one- and two-factor | `bermudan_exposure.hpp` |
| CS01 ladder; adjoint delta of a swap's CVA to the curve quotes | `cva_risk.hpp` |

Everything is in deflated units, as in the pricers: the expectation of a deflated value is its price today,
so expected exposure is already discounted and CVA is exposure weighted by default probabilities.

## Results (docs/RESULTS.md §11)

- One swap: closed-form CVA 30.62 bp, engine 30.53 bp.
- Eight-swap book: CVA 97.5 bp gross, 15.2 bp netted, 0.88 bp under a zero-threshold CSA with a 10-day
  margin period. Collateralised exposure spikes on coupon dates.
- Bermudan: CVA 38.5 bp against 30.6 bp for the underlying swap; DVA 0.6 bp against 15.1 bp.
- CVA delta to 37 curve quotes by adjoint in 2 ms, about 140 times faster than bumping.
- Under two factors the calibrated 10-year swap's CVA does not move with factor correlation; the book's
  netted CVA moves 15.0–15.6 bp and its gross CVA 91–103 bp.

## Tested

Bootstrap reprices every CDS; credit triangle; an arbitrageable spread curve is rejected. Engine = closed
form on coupon dates for payer and receiver; EPE + ENE = today's value of the remaining swap; a payer and
its receiver net to zero; identical trades do not net; netted ≤ gross. Collateral: residual exposure is
positive, grows with the margin period like its square root and with the threshold, and falls to zero at
maturity; an unreachable threshold is no agreement; a coupon inside the margin period leaves a spike. One
side's CVA exposure is the other's DVA exposure. Bermudan: grid and simulation reproduce the PDE price; the
first-date exposure is the price and is never negative; one exercise date is the European; deep in the
money is the swap's closed form. Adjoint CVA delta = central bump to 2e-10; the CS01 ladder adds up to a
parallel shift. Two factors: the engine reproduces the two-factor swaption strip on coupon dates; with the
second factor off it is the one-factor engine; factor correlation moves a curve-slope position and not the
calibrated swap; the regression-based Bermudan profile matches the one-factor grid with the second factor off.
