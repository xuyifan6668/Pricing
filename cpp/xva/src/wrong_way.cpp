#include "rl/xva/wrong_way.hpp"

#include "rl/core/error.hpp"
#include "rl/math/rootfind.hpp"
#include "rl/xva/xva.hpp"

#include <algorithm>
#include <cmath>

namespace rl::xva {

WrongWayResult wrongWayCva(const ExposureProfile& profile, const models::Lgm1F& model, const HazardCurve& hazard, double recovery, double beta) {
    const std::size_t nt = profile.times.size();
    RL_REQUIRE(nt >= 2 && profile.pathValue.size() == nt && profile.pathState.size() == nt, "wrongWayCva: the profile must be simulated with keepPaths");
    RL_REQUIRE(recovery >= 0.0 && recovery <= 1.0 && std::isfinite(beta), "wrongWayCva: recovery must be in [0, 1] and beta finite");
    const std::size_t n = profile.pathValue[1].size();
    const double paths = static_cast<double>(n);

    WrongWayResult out;
    out.times = profile.times;
    out.epe = profile.epe;
    out.epeGivenDefault.assign(nt, 0.0);
    out.hazardLevel.assign(nt, 0.0);
    out.survivalError.assign(nt, 0.0);
    out.independentCva = creditAdjustment(profile.times, profile.epe, hazard, recovery);

    std::vector<double> survival(n, 1.0), loss(n, 0.0), multiplier(n), weight(n), next(n);
    for (std::size_t j = 1; j < nt; ++j) {
        const double t = profile.times[j], dt = t - profile.times[j - 1];
        const double zeta = model.zeta(t), sd = std::sqrt(zeta), h = model.H(t);
        for (std::size_t p = 0; p < n; ++p) {
            const double x = profile.pathState[j][p];
            multiplier[p] = std::exp(beta * x / sd - 0.5 * beta * beta);
            weight[p] = std::exp(-h * x - 0.5 * h * h * zeta);
        }
        // The weights have expectation 1; dividing by their sample mean removes that source of sampling
        // error, so that with no dependence the calibration returns the curve's own hazard rates exactly.
        double meanWeight = 0.0;
        for (double w : weight) meanWeight += w;
        meanWeight /= paths;
        for (double& w : weight) w /= meanWeight;
        // E[w_j S_{j−1} exp(−a m_j Δt)] falls from E[w_j S_{j−1}] to 0 as a rises.
        const double target = hazard.survival(t);
        auto error = [&](double a) {
            double sum = 0.0;
            for (std::size_t p = 0; p < n; ++p) sum += weight[p] * survival[p] * std::exp(-a * multiplier[p] * dt);
            return sum / paths - target;
        };
        if (error(0.0) < 0.0)
            throw CalibrationError("wrongWayCva: survival to " + std::to_string(t) + "y cannot be matched: the earlier intervals already imply more default than the curve",
                                   {error(0.0)});
        double hi = 0.05;
        while (error(hi) > 0.0) {
            hi *= 2.0;
            RL_REQUIRE(hi < 1e4, "wrongWayCva: hazard level out of range");
        }
        math::RootOptions ro;
        ro.xTol = 1e-14;
        ro.maxIter = 200;
        const double a = math::brent(error, 0.0, hi, ro).x;
        out.hazardLevel[j] = a;
        out.survivalError[j] = error(a);

        double exposureInDefault = 0.0, defaultMass = 0.0;
        for (std::size_t p = 0; p < n; ++p) {
            next[p] = survival[p] * std::exp(-a * multiplier[p] * dt);
            const double pd = survival[p] - next[p];
            const double exposure = 0.5 * (std::max(profile.pathValue[j - 1][p], 0.0) + std::max(profile.pathValue[j][p], 0.0));
            loss[p] += exposure * pd;
            exposureInDefault += std::max(profile.pathValue[j][p], 0.0) * pd;
            defaultMass += pd;
        }
        out.epeGivenDefault[j] = defaultMass > 0.0 ? exposureInDefault / defaultMass : 0.0;
        survival.swap(next);
    }
    double sum = 0.0, sumSq = 0.0;
    for (double l : loss) {
        sum += l;
        sumSq += l * l;
    }
    const double mean = sum / paths;
    out.cva = (1.0 - recovery) * mean;
    out.standardError = (1.0 - recovery) * std::sqrt(std::max(sumSq / paths - mean * mean, 0.0) / paths);
    return out;
}

} // namespace rl::xva
