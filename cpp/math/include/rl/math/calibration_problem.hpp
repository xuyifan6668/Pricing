#pragma once
// The calibration-problem abstraction.
//
// Every calibration in this library has the same shape: a parameter vector z, a residual vector F(z; q)
// that vanishes (or is minimised) when the model reprices the instruments quoted at q, a scaling that makes
// residuals commensurable, and instrument ids for the report. Instances:
//    bootstrap::BootstrapProblem  z = all free nodes (+ turn jumps) of all curves; square; Newton.
//    vol::SabrNodeProblem         z = (α, ρ, ν) or (ρ, ν) with α from ATM; m > n; Levenberg–Marquardt.
//    models::LgmStripProblem      z = ζ increments; one swaption each; sequential 1-D roots.
//
// DESIGN: the concept lives in math, next to the double-only solvers that consume it, and not in bootstrap
// where it was first written. vol and models must
// be able to satisfy it without linking the bootstrap; bootstrap re-exports it with a using-declaration so
// no code changed meaning.
//
// DESIGN: residuals(z) are *raw* (quote units: rate, spread, bp of vol) and scaling(z) is separate, so the
// report shows what the desk reads while the solver sees commensurable numbers. A problem that folds the
// weights into its residuals makes "the 7y is off by 0.01 bp" unanswerable from the solver output.

#include "rl/math/lm.hpp" // math::Vector

#include <concepts>
#include <cstddef>
#include <string>
#include <vector>

namespace rl::math {

/// What every solver may assume about a calibration problem.
template <class P>
concept CalibrationProblem = requires(const P& p, const Vector& z) {
    { p.parameterCount() } -> std::convertible_to<std::size_t>;
    { p.residualCount() } -> std::convertible_to<std::size_t>;
    { p.residuals(z) } -> std::same_as<Vector>;
    { p.scaling(z) } -> std::same_as<Vector>;
    { p.instrumentIds() } -> std::convertible_to<const std::vector<std::string>&>;
};

} // namespace rl::math
