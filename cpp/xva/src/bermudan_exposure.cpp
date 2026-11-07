#include "rl/xva/bermudan_exposure.hpp"

#include "rl/core/error.hpp"
#include "rl/math/rng.hpp"

#include <algorithm>
#include <cmath>

namespace rl::xva {

namespace {

/// Continuation value at one exercise date on a uniform grid in the state.
struct Surface {
    double lo = 0.0, step = 1.0;
    std::vector<double> continuation;

    double operator()(double x) const noexcept {
        const double u = (x - lo) / step;
        const std::size_t n = continuation.size();
        if (u <= 0.0) return continuation.front();
        if (u >= static_cast<double>(n - 1)) return continuation.back();
        const std::size_t i = static_cast<std::size_t>(u);
        const double w = u - static_cast<double>(i);
        return (1.0 - w) * continuation[i] + w * continuation[i + 1];
    }
};

} // namespace

BermudanExposure bermudanExposure(const pricers::BermudanProblem& problem, const BermudanExposureOptions& opt) {
    const std::size_t n = problem.size();
    RL_REQUIRE(n > 0, "bermudanExposure: no exercise dates");
    RL_REQUIRE(opt.paths >= 100, "bermudanExposure: need at least 100 paths");
    RL_REQUIRE(opt.gridPoints >= 51 && opt.gridPoints % 2 == 1 && opt.stdDevs >= 4.0, "bermudanExposure: need an odd grid of at least 51 points over at least 4 standard deviations");
    const std::size_t m = static_cast<std::size_t>(opt.gridPoints);

    // --- Backward induction: C_i on a grid of ±stdDevs standard deviations of the state at date i ---
    std::vector<Surface> surfaces(n);
    std::vector<double> next; // V_{i+1} on date i + 1's grid
    for (std::size_t i = n; i-- > 0;) {
        const double sd = std::sqrt(problem.zeta(i));
        Surface& s = surfaces[i];
        s.lo = -opt.stdDevs * sd;
        s.step = 2.0 * opt.stdDevs * sd / static_cast<double>(m - 1);
        s.continuation.assign(m, 0.0);
        if (i + 1 < n) {
            // Trapezoid rule over date i + 1's grid against the Gaussian transition density. The step
            // is a small fraction of the transition's standard deviation, where the rule is very accurate.
            const Surface& up = surfaces[i + 1];
            const double var = problem.zeta(i + 1) - problem.zeta(i);
            RL_REQUIRE(var > 0.0, "bermudanExposure: the state has no variance between two exercise dates");
            const double norm = up.step / std::sqrt(2.0 * 3.141592653589793 * var);
            for (std::size_t a = 0; a < m; ++a) {
                const double x = s.lo + s.step * static_cast<double>(a);
                double sum = 0.0;
                for (std::size_t b = 0; b < m; ++b) {
                    const double d = up.lo + up.step * static_cast<double>(b) - x;
                    sum += next[b] * std::exp(-0.5 * d * d / var);
                }
                s.continuation[a] = sum * norm;
            }
        }
        next.resize(m);
        for (std::size_t a = 0; a < m; ++a) next[a] = std::max(problem.exerciseValue(i, s.lo + s.step * static_cast<double>(a)), s.continuation[a]);
    }

    BermudanExposure out;
    {
        // Today: one more Gaussian step, from x = 0 to the first exercise date.
        const Surface& first = surfaces[0];
        const double var = problem.zeta(0);
        double sum = 0.0;
        for (std::size_t b = 0; b < m; ++b) {
            const double x = first.lo + first.step * static_cast<double>(b);
            sum += next[b] * std::exp(-0.5 * x * x / var);
        }
        out.price = sum * first.step / std::sqrt(2.0 * 3.141592653589793 * var);
    }

    // --- Forward simulation under the optimal rule ----------------------------------------------------
    out.times.push_back(0.0);
    for (const auto& r : problem.rights()) out.times.push_back(r.time);
    const std::size_t nt = n + 1;
    out.epe.assign(nt, 0.0);
    out.ene.assign(nt, 0.0);
    out.epeStandardError.assign(nt, 0.0);
    out.optionEpe.assign(nt, 0.0);
    out.exercisedBy.assign(nt, 0.0);
    out.epe[0] = out.optionEpe[0] = out.price;

    std::vector<double> sd(n), sumSq(nt, 0.0);
    for (std::size_t i = 0; i < n; ++i) sd[i] = std::sqrt(std::max(problem.zeta(i) - (i ? problem.zeta(i - 1) : 0.0), 0.0));
    double payoff = 0.0, payoffSq = 0.0;
    math::Rng rng(opt.seed);
    for (std::size_t p = 0; p < opt.paths; ++p) {
        double x = 0.0, realised = 0.0;
        bool exercised = false;
        for (std::size_t i = 0; i < n; ++i) {
            x += sd[i] * rng.gaussian();
            const double h = problem.exerciseValue(i, x);
            double value;
            if (exercised) {
                value = h; // the remaining swap
            } else {
                const double hold = surfaces[i](x);
                if (h > 0.0 && h >= hold) {
                    exercised = true;
                    realised = h;
                    value = h;
                } else {
                    value = hold;
                    out.optionEpe[i + 1] += hold;
                }
            }
            const double pos = std::max(value, 0.0);
            out.epe[i + 1] += pos;
            sumSq[i + 1] += pos * pos;
            out.ene[i + 1] += std::min(value, 0.0);
            if (exercised) out.exercisedBy[i + 1] += 1.0;
        }
        payoff += realised;
        payoffSq += realised * realised;
    }
    const double paths = static_cast<double>(opt.paths);
    for (std::size_t i = 1; i < nt; ++i) {
        out.epe[i] /= paths;
        out.ene[i] /= paths;
        out.optionEpe[i] /= paths;
        out.exercisedBy[i] /= paths;
        out.epeStandardError[i] = std::sqrt(std::max(sumSq[i] / paths - out.epe[i] * out.epe[i], 0.0) / paths);
    }
    out.simulatedPrice = payoff / paths;
    out.simulatedPriceStandardError = std::sqrt(std::max(payoffSq / paths - out.simulatedPrice * out.simulatedPrice, 0.0) / paths);
    return out;
}

} // namespace rl::xva
