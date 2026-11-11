# xva

Counterparty credit and funding adjustments on top of the pricing models.

| Piece | File |
|---|---|
| Credit curve: hazard rates bootstrapped from CDS spreads | `credit.hpp` |
| Exposure engine: a netting set of swaps on simulated paths, with collateral | `exposure_engine.hpp` |
| CVA, DVA, FCA, FBA from exposure; one swap's exposure in closed form | `xva.hpp` |
| Exposure of a Bermudan swaption: an option until exercised, a swap afterwards | `bermudan_exposure.hpp` |
| CS01 ladder; adjoint delta of a swap's CVA to the curve quotes | `cva_risk.hpp` |

Everything is in deflated units, as in the pricers: the expectation of a deflated value is its price today,
so expected exposure is already discounted and CVA is exposure weighted by default probabilities.

## Results (docs/RESULTS.md §11)

- One swap: closed-form CVA 29.55 bp, engine 29.46 bp.
- Eight-swap book: CVA 93.9 bp gross, 14.9 bp netted, 0.89 bp under a zero-threshold CSA with a 10-day
  margin period. Collateralised exposure spikes on coupon dates.
- Bermudan: CVA 37.7 bp against 29.5 bp for the underlying swap; DVA 0.6 bp against 14.5 bp.
- CVA delta to 37 curve quotes by adjoint in 2 ms, about 140 times faster than bumping.

## Tested

Bootstrap reprices every CDS; credit triangle; an arbitrageable spread curve is rejected. Engine = closed
form on coupon dates for payer and receiver; EPE + ENE = today's value of the remaining swap; a payer and
its receiver net to zero; identical trades do not net; netted ≤ gross. Collateral: residual exposure is
positive, grows with the margin period like its square root and with the threshold, and falls to zero at
maturity; an unreachable threshold is no agreement; a coupon inside the margin period leaves a spike. One
side's CVA exposure is the other's DVA exposure. Bermudan: grid and simulation reproduce the PDE price; the
first-date exposure is the price and is never negative; one exercise date is the European; deep in the
money is the swap's closed form. Adjoint CVA delta = central bump to 2e-10; the CS01 ladder adds up to a
parallel shift.
