# bootstrap

The global multi-curve solve, its Jacobian, and the transform from node risk to par-quote risk.

| Piece | File |
|---|---|
| Quotes → dated instruments, one pillar each | `instrument_set.hpp`, `pillars.hpp` |
| The calibration problem: parameters, residuals, quotes | `bootstrap_problem.hpp` |
| Residual per instrument type (rate form) | `detail/residuals.ipp` |
| Initial guess, damped Newton, sequential bootstrap for comparison | `global_bootstrap.cpp` |
| Jacobian by adjoint sweep; ∂F/∂q; implicit-function par risk | `jacobian.cpp` |
| Repricing report | `report.hpp` |

## Why global

The 25Y and 30Y SOFR points are quoted through Fed Funds–SOFR basis swaps, and Fed Funds cashflows are
discounted on SOFR. Each curve needs the other, so all nodes of both curves are solved at once.

## Results (docs/RESULTS.md §1–2)

- 37 instruments repriced to 1e-10 bp in 4 Newton iterations, under linear log-discount and monotone
  convex interpolation.
- One sequential pass leaves 5 bp on the long basis swaps.
- Par risk by the implicit-function theorem matches bump-and-rebuild to 1e-6 of DV01.

## Tested

Repricing of every instrument; quadratic convergence from a flat guess and from zero; residuals
commensurable across instruments (1 bp of quote reads as 1 bp); adjoint Jacobian against finite
differences; the Jacobian's sparsity (SOFR-only instruments do not see Fed Funds nodes); bit-identical
results on repeat; the turn-of-year jump is two days wide.
