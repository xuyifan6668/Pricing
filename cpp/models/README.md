# models

The LGM one-factor and G2++ two-factor models, European swaptions under them, and their calibration.

| Piece | File |
|---|---|
| The model: H(t), ζ(t), deflated bond, Hull–White mapping | `lgm.hpp` |
| A swap as a list of zero-coupon bonds; exercise rights | `flows.hpp` |
| Swaption price under LGM (Jamshidian); market (Bachelier) price; implied vol | `lgm_swaption.hpp` |
| Bootstrap of ζ to a swaption strip; mean reversion marked to a second strip | `lgm_calibration.hpp` |
| The USD market used everywhere: bootstrapped curves + fitted cube | `usd_market.hpp` |
| G2++: covariance of the two states, deflated bond, forward-rate correlation | `g2pp.hpp` |
| Swaption price under G2++: closed form in one direction, quadrature in the other | `g2pp_swaption.hpp` |
| Cheyette (quasi-Gaussian) model with local volatility: two states, closed-form bonds | `cheyette.hpp` |
| G2++ shape from the ATM matrix; scale bootstrapped to a strip; two swaptions matched at a given ρ | `g2pp_calibration.hpp` |

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

## The second factor

Two driftless Gaussian states with covariance C(t); the deflated bond is `D(T) exp(−H(T)·x − ½ H(T)ᵀ C(t) H(T))`.
The shape (two mean reversions, σ₂/σ₁, ρ) sets how rates of different maturities and different dates are
correlated; a common volatility scale s(t) sets the level and is bootstrapped to a strip like ζ.

Tested: the deflated bond is a martingale; the bond option reproduces Brigo–Mercurio's formula; the swaption
formula equals two-dimensional integration, has converged in the quadrature order and satisfies put–call
parity; with the second factor off, or with two perfectly correlated factors of one speed, prices equal the
one-factor model's; the fit recovers a two-factor model's own vols; the scale reprices a strip to 1e-6 bp
for any shape and reduces to the one-factor calibration; an inconsistent strip is flagged.
