# models

The LGM one-factor model, European swaptions under it, and calibration to a co-terminal strip.

| Piece | File |
|---|---|
| The model: H(t), ζ(t), deflated bond, Hull–White mapping | `lgm.hpp` |
| A swap as a list of zero-coupon bonds; exercise rights | `flows.hpp` |
| Swaption price under LGM (Jamshidian); market (Bachelier) price; implied vol | `lgm_swaption.hpp` |
| Bootstrap of ζ to a swaption strip | `lgm_calibration.hpp` |
| The USD market used everywhere: bootstrapped curves + fitted cube | `usd_market.hpp` |

## The model in four lines

State x with x₀ = 0 and x_t ~ N(0, ζ(t)), no drift under the model's numeraire. Deflated bond:
`P(t,T;x)/N(t,x) = D(T) exp(−H(T) x − ½ H(T)² ζ(t))`. H is set by the mean reversion κ; ζ is the
accumulated variance. The discount curve D is an input, so the curve is fitted by construction.

## Tested

- Martingale property: E[deflated bond] = D(T) for every κ, t, T.
- The bond decomposition reprices a swap off the curve to 1e-13 of notional.
- Jamshidian closed form = numerical integration of the payoff; put–call parity; zero-variance limit.
- A one-period swaption reproduces the Hull–White zero-bond option formula.
- The strip reprices to 1e-6 bp for κ from −1% to 10%; κ changes ζ but not the calibrated prices.
- An inconsistent strip (one needing variance to fall) is flagged and left visibly mispriced.
