#pragma once
// One-dimensional root finding, on double only.
//
// Why not templated on T: adjoints never tape through a solver's iterations.
// If x* solves f(x*; θ) = 0 then dx*/dθ = −(∂f/∂θ)/(∂f/∂x) at x* (implicit-function theorem), and
// provides aad::implicitRoot() that records exactly that as one tape node. Taping Brent's
// bisection/secant/inverse-quadratic branches would record a derivative that is only approximately
// right and costs 30× the memory. Callers under T = Real solve on value() and wrap the result.

#include <functional>
#include <limits>

namespace rl::math {

struct RootOptions {
    double xTol = 1e-12;   ///< Stop when the bracket (Brent) or step (Newton) is below this, absolute.
    double fTol = 0.0;     ///< Also stop when |f| ≤ fTol (0 = ignore).
    int maxIter = 100;
};

struct RootResult {
    double x = std::numeric_limits<double>::quiet_NaN();
    double f = std::numeric_limits<double>::quiet_NaN();
    int iterations = 0;
    bool converged = false;
};

/// Brent's method on a bracket [a, b] with f(a)·f(b) ≤ 0. Throws InputError if not bracketed,
/// NumericalError if maxIter is exhausted.
RootResult brent(const std::function<double(double)>& f, double a, double b, RootOptions opt = {});

/// Newton–Raphson with bisection safeguard: `fdf(x)` returns {f(x), f'(x)}; iterates from x0 and
/// never leaves [lo, hi], bisecting when a Newton step would exit the bracket or fails to shrink it.
/// Throws as brent().
RootResult newtonSafeguarded(const std::function<std::pair<double, double>(double)>& fdf, double x0, double lo, double hi,
                             RootOptions opt = {});

/// Expand [a, b] geometrically about its centre until f changes sign (at most maxExpansions doublings).
/// Returns false if no bracket was found. Used before brent() when a bracket is only a guess.
bool expandBracket(const std::function<double(double)>& f, double& a, double& b, int maxExpansions = 50);

} // namespace rl::math
