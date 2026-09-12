#pragma once
// Standard normal distribution: density, cdf, inverse cdf. Used by Bachelier/Black, the
// MC gaussian draws and the LGM analytic swaption — all hot.

#include "rl/math/scalar.hpp"

namespace rl::math {

/// φ(x) = exp(−x²/2)/√(2π)
template <Scalar T>
T normalPdf(const T& x);

/// Φ(x) = ½ erfc(−x/√2). Accurate in both tails (erfc does not cancel for x ≪ 0 the way ½(1+erf) does).
template <Scalar T>
T normalCdf(const T& x);

/// Φ⁻¹(p) for p ∈ (0, 1). Wichura (1988) algorithm AS241, PPND16 (relative accuracy ≈ 1e-16).
template <Scalar T>
T normalInverseCdf(const T& p);

} // namespace rl::math
