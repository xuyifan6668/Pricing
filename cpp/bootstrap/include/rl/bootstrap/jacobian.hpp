#pragma once
// The bootstrap Jacobian and the zero-to-par risk transform. Sensitivities (adjoint or bumped) stop at
// the curve nodes z; the implicit-function theorem carries them to the par quotes q a trader hedges with.
//
// Notation: n instruments, n parameters. F(z; q) ∈ ℝⁿ is the residual vector, z the node ordinates
// y = −log P at the pillars followed by turn jumps, curve by curve; q the quotes (decimal rates;
// futures in price).

#include "rl/bootstrap/bootstrap_problem.hpp"
#include "rl/math/lm.hpp"

#include <functional>

namespace rl::bootstrap {

/**
 * ∂F/∂z, n × n, by adjoint differentiation: the residuals are recorded once on the tape with the nodes
 * as inputs, and one reverse sweep per residual gives its row. Exact to rounding under either
 * interpolation scheme, and exactly sparse: a residual that never read a node has a zero there.
 *
 * Structure: close to lower-triangular within each curve's block (a 5y swap does not see the 7y node),
 * with the basis-swap and Fed Funds rows populated in BOTH curve blocks. That coupling is the picture of
 * why a sequential bootstrap fails.
 */
math::Matrix dFdz(const BootstrapProblem& problem, const math::Vector& z);

/// ∂F/∂q, n × n. With rate-form residuals this is −I (−1/100 for a future, quoted in price). It is
/// diagonal only because each quote is a contract term of exactly one instrument.
math::Matrix dFdq(const BootstrapProblem& problem, const math::Vector& z);

/**
 * Zero risk → par risk by the implicit-function theorem.
 *
 * At the solution F(z*(q); q) = 0 for all q, so ∂F/∂z · dz/dq + ∂F/∂q = 0 and, for a price P(z) with
 * node risk g = ∂P/∂z,
 *     ∂P/∂q = −gᵀ (∂F/∂z)⁻¹ ∂F/∂q.
 * Computed as one transposed solve (∂F/∂z)ᵀ w = g followed by −(∂F/∂q)ᵀ w: one solve per price, not
 * one per quote, and no inverse.
 *
 * The theorem differentiates the equation, not the algorithm, so the result is as accurate as the
 * bootstrap is converged. Units: currency per unit of quote; multiply by 1e-4 for "per basis point".
 */
math::Vector zeroToParRisk(const math::Vector& dPdz, const math::Matrix& dFdzAtSolution, const math::Matrix& dFdqAtSolution);

// --- helpers ---------------------------------------------------------------------------------------

/// A PV as a function of the curve set (what a pricer is, from the bootstrap's point of view).
using PvFunction = std::function<double(const market::CurveSet<double>&)>;

/// ∂P/∂z by central differences over the node ordinates (a cross-check for the adjoint route).
math::Vector zeroRiskByBump(const BootstrapProblem& problem, const math::Vector& z, const PvFunction& pv, double step = 1e-6);

/// dFdz, dFdq and zeroToParRisk in one call; returns ∂P/∂q per unit of quote, aligned with problem.instrumentIds().
math::Vector parRisk(const BootstrapProblem& problem, const math::Vector& z, const math::Vector& dPdz);

/// The brute-force reference: for each quote i, bump by `bump` (quote units), re-solve, reprice, difference.
/// O(n) bootstraps; test and demo only.
math::Vector parRiskByRebuild(const BootstrapProblem& problem, const math::Vector& z, const PvFunction& pv, double bump = 1e-4);

} // namespace rl::bootstrap
