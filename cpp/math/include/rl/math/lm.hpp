#pragma once
// Levenberg–Marquardt least squares over Eigen, double only (see rootfind.hpp for why solvers are
// never templated on the AAD scalar). Used by the SABR per-node fit and available to the
// global bootstrap as a fallback when Newton fails.

#include <Eigen/Dense>

#include <functional>

namespace rl::math {

using Vector = Eigen::VectorXd;
using Matrix = Eigen::MatrixXd;

struct LmOptions {
    double lambda0 = 1e-3;         ///< Initial damping.
    double lambdaUp = 10.0;        ///< Multiply on a rejected step.
    double lambdaDown = 0.1;       ///< Multiply on an accepted step.
    double gradTol = 1e-10;        ///< Stop when ‖Jᵀr‖∞ ≤ gradTol.
    double stepTol = 1e-12;        ///< Stop when ‖Δx‖ ≤ stepTol (1 + ‖x‖).
    double costTol = 1e-16;        ///< Stop when the relative cost decrease is below this.
    int maxIter = 200;
    bool useMarquardtScaling = true; ///< Damp with diag(JᵀJ) rather than the identity.
    Vector lowerBounds;            ///< Optional box constraints (empty = none); enforced by clamping.
    Vector upperBounds;
};

enum class LmStatus { Converged, GradientSmall, StepSmall, CostSmall, MaxIterations, Failed };
const char* toString(LmStatus s);

struct LmResult {
    Vector x;
    Vector residuals;
    double cost = 0.0;    ///< ½ ‖r‖²
    int iterations = 0;
    int evaluations = 0;
    LmStatus status = LmStatus::Failed;
};

using ResidualFn = std::function<Vector(const Vector& x)>;
using JacobianFn = std::function<Matrix(const Vector& x)>; // m × n, J_ij = ∂r_i/∂x_j

/// Minimise ½ ‖r(x)‖² from x0. If `jac` is empty, a central finite-difference Jacobian is used
/// (step h_j = 1e-6 (1 + |x_j|)).
LmResult levenbergMarquardt(const ResidualFn& residuals, const JacobianFn& jac, const Vector& x0, LmOptions opt = {});

/// Central-difference Jacobian helper (exposed for tests and for checking analytic Jacobians).
Matrix finiteDifferenceJacobian(const ResidualFn& residuals, const Vector& x, double relStep = 1e-6);

} // namespace rl::math
