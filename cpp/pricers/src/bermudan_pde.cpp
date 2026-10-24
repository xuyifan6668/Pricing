#include "rl/pricers/bermudan_pde.hpp"

#include "rl/core/error.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace rl::pricers {

namespace {

/// One θ-step of V_ζ + ½V_xx = 0 backward over a clock interval `dz`:
///   (I − θ λ δ²) V_new = (I + (1−θ) λ δ²) V_old,   λ = dz / (2 h²),
/// θ = ½ is Crank–Nicolson, θ = 1 fully implicit. Edge rows hold V fixed (V_xx = 0 ⇒ V_ζ = 0).
/// Thomas algorithm with caller-owned scratch: no allocation in the time loop.
void thetaStep(std::vector<double>& v, double lambda, double theta, std::vector<double>& rhs, std::vector<double>& cp) {
    const std::size_t n = v.size();
    const double e = (1.0 - theta) * lambda;
    rhs[0] = v[0];
    rhs[n - 1] = v[n - 1];
    for (std::size_t j = 1; j + 1 < n; ++j) rhs[j] = e * v[j - 1] + (1.0 - 2.0 * e) * v[j] + e * v[j + 1];

    const double a = -theta * lambda, b = 1.0 + 2.0 * theta * lambda; // sub/super = a, diagonal = b
    cp[0] = 0.0; // first row is the identity
    v[0] = rhs[0];
    for (std::size_t j = 1; j + 1 < n; ++j) {
        const double m = b - a * cp[j - 1];
        cp[j] = a / m;
        v[j] = (rhs[j] - a * v[j - 1]) / m;
    }
    v[n - 1] = rhs[n - 1];
    for (std::size_t j = n - 1; j-- > 1;) v[j] -= cp[j] * v[j + 1];
}

} // namespace

PdeResult priceBermudanPde(const BermudanProblem& problem, const PdeOptions& opt) {
    const std::size_t n = problem.size();
    RL_REQUIRE(n > 0, "priceBermudanPde: no exercise dates");
    RL_REQUIRE(opt.spaceSteps >= 10 && opt.timeSteps >= 1 && opt.stdDevs > 0.0, "priceBermudanPde: invalid grid options");

    std::vector<double> zetas(n);
    for (std::size_t i = 0; i < n; ++i) zetas[i] = problem.zeta(i);
    const double zetaEnd = zetas.back();

    const int m = opt.spaceSteps + (opt.spaceSteps % 2); // even number of steps ⇒ x = 0 is the middle node
    const double xMax = opt.stdDevs * std::sqrt(zetaEnd);
    const double h = 2.0 * xMax / m;
    const std::size_t points = static_cast<std::size_t>(m) + 1;

    PdeResult res;
    res.spaceSteps = m;
    res.grid.resize(points);
    for (std::size_t j = 0; j < points; ++j) res.grid[j] = -xMax + h * static_cast<double>(j);
    res.boundaryState.assign(n, std::numeric_limits<double>::quiet_NaN());
    res.boundaryRate.assign(n, std::numeric_limits<double>::quiet_NaN());

    std::vector<double> v(points, 0.0), exercise(points), rhs(points), cp(points);

    for (std::size_t i = n; i-- > 0;) {
        // Exercise decision at date i.
        for (std::size_t j = 0; j < points; ++j) exercise[j] = problem.exerciseValue(i, res.grid[j]);
        // The boundary is where exercise − continuation changes sign; take the crossing nearest x = 0.
        double best = std::numeric_limits<double>::infinity();
        for (std::size_t j = 0; j + 1 < points; ++j) {
            const double d0 = exercise[j] - v[j], d1 = exercise[j + 1] - v[j + 1];
            if ((d0 > 0.0) == (d1 > 0.0)) continue;
            const double x = res.grid[j] + h * d0 / (d0 - d1);
            if (std::fabs(x) < std::fabs(best)) best = x;
        }
        if (std::isfinite(best)) {
            res.boundaryState[i] = best;
            res.boundaryRate[i] = problem.swapRate(i, best);
        }
        for (std::size_t j = 0; j < points; ++j) v[j] = std::max(v[j], exercise[j]);

        // Roll back to the previous exercise date (or to today).
        const double dz = zetas[i] - (i ? zetas[i - 1] : 0.0);
        if (dz <= 0.0) continue; // zero variance on this interval: nothing diffuses
        const int steps = std::max(4, static_cast<int>(std::ceil(opt.timeSteps * dz / zetaEnd)));
        res.timeSteps += steps;
        const double lambda = dz / steps / (2.0 * h * h);
        int k = 0;
        if (opt.rannacher) {
            for (int s = 0; s < 4; ++s) thetaStep(v, 0.5 * lambda, 1.0, rhs, cp);
            k = 2;
        }
        for (; k < steps; ++k) thetaStep(v, lambda, 0.5, rhs, cp);
    }
    res.value = v[static_cast<std::size_t>(m / 2)];
    res.values = v;
    return res;
}

} // namespace rl::pricers
