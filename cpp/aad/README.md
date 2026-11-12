# aad

Tape-based adjoint algorithmic differentiation: the active scalar `Real` and its `Tape`.

- `Real` holds a value and a tape node. Every arithmetic operator and elementary function records its
  local partial derivatives when it is evaluated. Operations between constants record nothing.
- `Tape::propagate()` walks the record backwards: `adjoint[parent] += partial × adjoint[node]`. One sweep
  gives the derivative of one output with respect to every input.
- `Real` has `==` but no ordering operators. A branch on an active value is written on `value(x)`, so every
  regime decision that the tape will freeze is visible in the code.
- The tape is thread-local. `mark()` / `rewindTo()` let a nested recording (the bootstrap Jacobian) clean
  up after itself without disturbing an outer one.
- `eigen_traits.hpp` lets `Eigen::Matrix<Real>` do products and sums. Linear solves stay in `double`;
  derivatives cross a solve by the implicit-function theorem (see `bootstrap`).

Used for: the bootstrap Jacobian (`bootstrap::dFdz`) and swap par-quote risk
(`pricers::swapParDeltaAdjoint`), about 80× faster than bumping on a 200-swap book.

Tests: every operator's adjoint against central differences; a 1000-input gradient from one sweep;
identical values under `double` and `Real`.
