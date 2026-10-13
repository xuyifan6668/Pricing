# rl::vol

Swaption volatility machinery: Bachelier and shifted Black-76 with greeks and implied vols, Hagan SABR
in its normal and lognormal forms with a shift, α from the ATM constraint, exact gradients, cube queries, the
per-node SABR calibration problem, and a static-arbitrage diagnostic. Runs ~1e2 times per cube build, once per
strip swaption per LGM calibration, and ~1e4 times per cube-vega ladder.

The judgement calls are in `normalToShiftedLognormalVol` and `shiftedLognormalToNormalVol` in
`src/conversion.cpp`, and `SabrNodeProblem::scaling`, `sabrInitialGuess` and `fitSabrSmile` in `src/sabr_fit.cpp`.
Everything else is closed forms with the source cited above the code.

## Contract

- Prices are undiscounted, per unit annuity. A payer swaption is a call on the forward swap rate, a receiver a put.
- Active set: the forward, the vols, α, ρ and ν are `T`; strike, expiry, β and shift are passive `double`.
- `bachelierPrice` accepts any sign of F and K. `blackPrice` requires F + s > 0 (`InputError`); K + s ≤ 0 prices
  intrinsic for the call and zero for the put, independent of the vol.
- Implied vols are `double`-only and solve on the out-of-the-money time value: an exact bracket and Newton on
  log price for Bachelier, a geometric bracket, Newton on log price and a Brent fallback for Black. They return 0 at
  intrinsic and throw below intrinsic, at the no-arbitrage bound, or when K + s ≤ 0.
- `sabrNormalVol` is Hagan et al. (2002) eq. (B.69) with ζ in integral form. The β = 0 path takes no logs or powers,
  so negative rates work. For β > 0 it may return a non-positive vol deep in the low-strike wing: that is the
  expansion failing, reported by `arbitrage.hpp`, not an exception.
- `alphaFromAtmNormalVol` returns the smallest positive root of the ATM cubic and then takes one Newton step in `T`
 . It throws `NumericalError` when the ATM vol exceeds the cubic's maximum.
- `sabrNormalVolGradient` and `sabrNormalVolGradientAtmConstrained` are exact: a `math::Dual` through the same code.
- `cubeParams` interpolates the cube (`market::SwaptionCube`) and re-solves α at the caller's forward, so the cube
  reprices the interpolated ATM vol at any (expiry, tenor, forward).
- `SabrNodeProblem` satisfies `math::CalibrationProblem`. Its residuals are raw, in bp of normal vol; the weights are
  `scaling(z)`. With `alphaFromAtm` the parameter vector is (ρ, ν).
- `checkDensity` scans away from the forward and reports the nearest strike on each side where the density is below
  −1e-6 or the vol is non-positive. The second difference uses h = 0.5 bp.
- Templates are declared in the headers, defined in `detail/*.ipp`, and instantiated for `double` in `src/vol.cpp`
  and for `rl::aad::Real` in `src/instantiate_real.cpp`. Add a line to the matching `RL_VOL_INSTANTIATE_*` macro
  when you add a template.
- `vol` links `core`, `math` and `market` only, never `linear`: forwards are arguments.

## Invariants (tested)

- Put–call parity to 1e-15 under both models, for negative forwards and with shifts. Haug's Black-76 reference
  value 1.7011. A shift is exactly a displaced forward and strike.
- Greeks match central differences. Implied vols round-trip to 1e-12 on OTM options out to ±6 standard deviations
  (Bachelier) and ±3 (Black).
- SABR limits: the β = 0 ATM closed form, the ν = 0 local-vol term with its α² correction, continuity and the skew
  through the money for β ∈ {0, 0.5, 1}, and continuity across the ζ/x̂(ζ) series threshold. With ρ = 0 the smile is
  symmetric in F ↔ K. The β = 1, ν = 0 normal form matches the exact Black-to-normal conversion to 1e-5. The
  lognormal and normal forms agree within 0.1 bp at 1Y.
- Gradients match finite differences, and the ATM-constrained gradient vanishes at the money. `alphaFromAtm`
  inverts to 1e-14 or throws only when a brute-force scan confirms the vol is unreachable.
- `Real` results are bit-identical to `double` and active, including through the ATM constraint.
- The cube reproduces node smiles, and reprices ATM off-grid at any forward.
- The Bachelier density integrates to 1 with mean F. The benign 1Y×10Y node is arbitrage-free over ±4 standard
  deviations. The stressed 10Y×30Y node goes negative near 2.09%, within one coarse grid step of the fine-grid answer.
- Also tested: the conversion round trips, the weighting, the initial guess, synthetic recovery to
  1e-6, the rounded-quote fit under 0.1 bp, frown failure detection, and the 81-node USD cube.

## Threading

Pure functions of their inputs. `SabrNodeProblem` is immutable after construction, so nodes can be fitted in
parallel with no shared state. The `Real` instantiations write to the calling thread's tape.
