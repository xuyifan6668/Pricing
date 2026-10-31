# Design

Short on purpose. The library is as much software as the pricing needs and no more.

## Shape

```
core → math → aad → market → instruments → linear → bootstrap → vol → models → pricers
```

Each module depends only on modules to its left. The first eight are the curve and volatility layer; the
last two are the model and the products priced under it.

Two styles of code, deliberately:

- **Curve layer (`math` to `bootstrap`, `vol` formulas): templated on the scalar type.** Every function
  that turns curve nodes into a price compiles for `double` and for `aad::Real`, the tape-recording type.
  That one decision is what makes the bootstrap Jacobian and curve risk available by adjoint
  differentiation with no second implementation to keep in sync.
- **Model and pricer layer (`models`, `pricers`): plain `double`, free functions, small structs.** A
  calibration is a root search and a pricer is a loop; neither needs to be differentiated through, because
  risk on model prices is taken by bumping the market and re-calibrating.

## Decisions that matter

**Curves are bootstrapped globally, not one at a time.** The SOFR long end is quoted through Fed
Funds–SOFR basis swaps while Fed Funds cashflows are discounted on SOFR. One Newton solve over all nodes of
both curves reprices everything; a sequential pass leaves 5 bp (docs/RESULTS.md §1).

**Residuals are in rate units.** `model quote − market quote` for every instrument. A 1-week deposit and a
30-year basis swap are then directly comparable, the convergence number is the number in the repricing
report, and ∂F/∂q = −I.

**The bootstrap Jacobian comes from the adjoint tape.** The residuals are recorded once with the nodes as
inputs and swept once per row. It is exact for either interpolation scheme. A hand-derived Jacobian through
monotone-convex interpolation has to carry the dependence of each segment's breakpoint on the nodes; that
term is easy to drop and dropping it leaves every price right.

**Risk crosses the bootstrap by the implicit-function theorem, not by taping Newton.** At the solution
F(z(q), q) = 0, so dz/dq = −(∂F/∂z)⁻¹ ∂F/∂q. Node risk becomes par-quote risk with one transposed linear
solve. This differentiates the equation rather than the algorithm, so it does not depend on how the solver
got there.

**The AAD scalar has no ordering operators.** A branch on an active value must be written on its plain
value, which makes every "which regime am I in" decision visible in the code. Equality exists because
Eigen's product kernels need it.

**Normal vol is the working convention.** The cube stores normal vols and SABR parameters; the LGM
calibration targets normal vols. Lognormal quotes are converted on the way in, exactly, by price matching.

**SABR: β and the shift are policy, α comes from the ATM quote, ρ and ν are fitted.** β and ρ both produce
skew; fitting both is ill-conditioned. Solving α from the ATM vol inside every evaluation means each node
reprices ATM exactly and the search is two-dimensional.

**LGM, not Hull–White, as the parameterisation.** Same model. In LGM the deflated bond is a single
exponential in a driftless Gaussian state, so: the model fits the discount curve by construction, a
swaption depends on the volatility only through ζ at its expiry (calibration is a sequence of 1-D root
searches), the PDE is the heat equation, and Monte Carlo steps are exact.

**A swap is handed to the model as a list of zero-coupon bonds.** Fixed coupons are bonds; a compounded
overnight coupon over [s, e] is +1 at s and −1 at e. Everything the model prices is then Σ aₖ P(t, Tₖ), and
the decomposition is tested to reprice the swap off today's curve exactly.

**Calibration is to the trade's own co-terminal strip at the trade's strike.** Those are the Europeans the
Bermudan contains. An inconsistent strip (one that would need variance to fall) is reported and left
visibly mispriced, not smoothed over.

**Monte Carlo reports two numbers, never one.** An out-of-sample lower bound and a dual upper bound. The
gap between them measures the quality of the exercise rule without reference to any other pricer.

**No fast-math.** Reassociation would make results depend on compiler and inlining; determinism is worth
more than the speed here.

## What was deliberately left out

Python bindings, an Excel add-in, a parallel risk engine, a generic product/model registry, payoff
scripting. None of them would show anything about a model or a product. Risk on the Bermudan is bump and
re-calibrate; it is 37 PDE prices at 2 ms each.

## Error handling

Exceptions (`rl::InputError`, `rl::NumericalError`, `rl::CalibrationError`) with a context chain, thrown at
validation time. Inner loops do not throw. A calibration that fails says which instrument and why.
