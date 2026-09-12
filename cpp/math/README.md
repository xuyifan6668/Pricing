# rl::math

The numerical toolbox: `Scalar` concept, interpolation strategies, root finding, Levenberg–Marquardt,
normal distribution, Gauss–Hermite quadrature, tridiagonal solve, MT19937-64 and the seed policy,
Sobol (Joe–Kuo) and the Brownian bridge.

Sources are cited in the code. Curves use `Linear<T>` on log-discount or `MonotoneConvex<T>`
(Hagan–West); `NaturalCubic<T>` (vol cube
only) is implemented.

## Contract

- Routines templated on `T` (`interpolation`, `normal`, `tridiagonal`, `GaussHermite::integrate`)
  compile and give **bit-identical values** for `T = double` and `T = rl::aad::Real`. Templates are
  declared in `.hpp`, defined in `include/rl/math/detail/*.ipp`, and explicitly instantiated for
  `double` in `src/*.cpp` and for `Real` in `cpp/aad/src/instantiate_math.cpp` (target
  `rl::math_real`). Consumers that price under both scalars link `rl::math_real`.
- Solvers (`brent`, `newtonSafeguarded`, `levenbergMarquardt`) and random-number machinery
  (`Rng`, `Sobol`, `BrownianBridge`) are `double` only. Adjoints never tape through a solver; they
  cross it via the implicit-function theorem.
- `rl::aad::Real` has **no relational operators**. Branch on `value(x)`; use `rl::math::max/min`.
- `Dual<T, N>` (`dual.hpp`) is a forward-mode dual over `double` or `Real` for small local Jacobians (the SABR
  gradient). It satisfies `Scalar` and, like `Real`, has no relational operators. It is not the risk engine.
- `CalibrationProblem` (`calibration_problem.hpp`) is the concept every calibration satisfies: the bootstrap,
  the SABR node, the LGM strip. Residuals are raw; the scaling is separate.
- `Rng::uniform01()` is strictly inside (0, 1); `Rng::gaussian()` is the inverse cdf of that
  uniform, so the stream is bit-identical across MSVC, libstdc++ and libc++.
- `Sobol` point 0 is the origin (Joe–Kuo convention). Pass `skip ≥ 1` before feeding an inverse cdf.

## Invariants (all tested)

- Interpolators reproduce nodes; `MonotoneConvex` satisfies `∫ f = Δy` on every interval and is
  continuous; `NaturalCubic` has zero second derivative at both ends.
- `cdf(inverseCdf(p)) = p` to 1e-15 relative on `[1e-300, 1 − 1e-16]`.
- Gauss–Hermite is exact for polynomials of degree `2n − 1`; `Σ w = 1`.
- `multiply(solve(d)) == d`; Thomas matches Eigen's dense solve on dominant systems.
- The first `2^m` Sobol points are stratified per coordinate and form a `(0, m, 2)`-net in the
  first two dimensions.
- The Brownian bridge's linear map has covariance `min(tᵢ, tⱼ)` exactly.
- `Dual` tangents match closed forms and central differences for every elementary function, with identical values
  under `double` and `Real`.

## Threading

Everything is immutable after construction except `Rng`, `Sobol` (one per thread/path stream) and
the thread-local tape that `Real` writes to.
