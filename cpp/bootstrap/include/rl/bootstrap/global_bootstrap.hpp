#pragma once
// The global multi-curve solve: Newton on all free nodes of all curves at once, with a step clamp, a
// line search, and loud failure. `solveSequential` is the textbook curve-by-curve bootstrap, kept to
// show why it cannot reprice a coupled instrument set.
//
// The solver works in double. Risk does not differentiate through the Newton iterations: par risk
// crosses the solve by the implicit-function theorem (jacobian.hpp).

#include "rl/bootstrap/bootstrap_problem.hpp"
#include "rl/math/lm.hpp" // math::Vector / math::Matrix

#include <string>
#include <vector>

namespace rl::bootstrap {

struct BootstrapOptions {
    double toleranceBp = 1e-10;     ///< converged when max_i |scale_i · r_i| < toleranceBp
    int maxIterations = 50;
    bool lineSearch = true;         ///< backtracking on ‖S F‖₂ (halve the step until it decreases; at most 20 halvings)
    double maxStep = 0.05;          ///< clamp each Newton step in y = −log P to ±maxStep (5% in "rate × time")
    bool analyticJacobian = true;   ///< adjoint dFdz(); false → finiteDifferenceJacobian() (ten times slower; a cross-check)
    double fdStep = 1e-7;           ///< central-difference step in y for the FD Jacobian
    bool rescaleEachIteration = true; ///< recompute residualScale at every iterate, or freeze at the guess
    bool verbose = false;           ///< push one line per iteration into BootstrapResult::log
};

struct BootstrapResult {
    math::Vector z;                 ///< solution (or last iterate if !converged)
    math::Vector residuals;         ///< raw residuals at z
    math::Vector scaledResiduals;   ///< in bp of quote
    double maxAbsBp = 0.0;          ///< ‖scaledResiduals‖∞
    int iterations = 0;
    bool converged = false;
    std::vector<std::string> log;   ///< per-iteration diagnostics when verbose; always the failure reason if !converged
};

/// Initial guess: each curve flat at the rate of its shortest instrument (y_i = r₀ t_i). A curve pinned
/// only by basis swaps starts from the reference curve's front rate plus the shortest spread.
/// Production would seed from yesterday's solution instead, which has the same layout: pass it as z.
math::Vector initialGuess(const BootstrapProblem& problem);

/**
 * Damped Newton on F(z) = 0 with scaling S (residuals in basis points of quote):
 *
 *   solve   (S J_k) Δ = −S F(z_k)       J_k = ∂F/∂z at z_k (adjoint, or finite differences)
 *   clamp   Δ_i ← sign(Δ_i) · min(|Δ_i|, maxStep)
 *   search  λ ∈ {1, ½, ¼, …}: accept the first λ with ‖S F(z_k + λΔ)‖₂ < ‖S F(z_k)‖₂
 *   stop    when max_i |S_ii F_i| < toleranceBp, or maxIterations, or no λ reduces the residual
 *
 * Converges quadratically: four iterations from the flat guess on the USD set. A singular Jacobian
 * throws NumericalError naming the node no instrument pins; non-convergence returns converged = false
 * with the reason in `log`.
 */
BootstrapResult solve(const BootstrapProblem& problem, const BootstrapOptions& options = {});

/**
 * One pass of the sequential bootstrap: for each curve in turn, solve that curve's instruments over that
 * curve's nodes with every other curve held where it currently is.
 *
 * On the USD set this reprices every SOFR swap and future and every Fed Funds instrument, but the 25Y
 * and 30Y basis swaps that pin the SOFR long end were solved against the *guessed* Fed Funds curve and
 * are left about 5 bp off. Iterating the pass would converge slowly; the global Newton is the right tool.
 */
BootstrapResult solveSequential(const BootstrapProblem& problem, const BootstrapOptions& options = {});

// --- helpers ------------------------------------------------------------------

/// Central-difference Jacobian of the raw residual vector, ∂F_i/∂z_j ≈ (F_i(z + h e_j) − F_i(z − h e_j))/(2h).
/// O(n) residual evaluations of O(n) each. For testing dFdz and as the analyticJacobian = false fallback.
math::Matrix finiteDifferenceJacobian(const BootstrapProblem& problem, const math::Vector& z, double step = 1e-7);

/// Scaled residual vector S·F(z) in bp, using residualScale: the convergence quantity.
math::Vector scaledResiduals(const BootstrapProblem& problem, const math::Vector& z);

} // namespace rl::bootstrap
