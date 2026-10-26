#include "rl/pricers/exposure.hpp"

#include "rl/core/error.hpp"
#include "rl/math/normal.hpp"
#include "rl/math/rng.hpp"
#include "rl/models/flows.hpp"
#include "rl/models/lgm_swaption.hpp"

#include <algorithm>
#include <cmath>
#include <variant>

namespace rl::pricers {

std::vector<ExposurePoint> swapExposureProfile(const models::Lgm1F& model, const instruments::Swap& swap, const market::CurveSet<double>& curves,
                                               const ExposureOptions& opt) {
    RL_REQUIRE(opt.paths >= 100 && opt.quantile > 0.5 && opt.quantile < 1.0, "swapExposureProfile: need ≥ 100 paths and a quantile in (0.5, 1)");
    const market::YieldCurve<double>& disc = curves.discount();
    std::vector<ExposurePoint> out;
    math::Rng rng(opt.seed);
    const double zq = math::normalInverseCdf<double>(opt.quantile);
    for (const auto& cf : swap.fixedLeg().cashflows()) {
        const core::Date start = std::get<instruments::FixedCoupon>(cf).accrualStart;
        const double t = disc.time(start);
        if (t <= 0.0) continue;
        const std::vector<models::BondFlow> flows = models::swapAsBonds(swap, curves, start);
        const double zeta = model.zeta(t), sd = std::sqrt(zeta);
        ExposurePoint p;
        p.time = t;
        p.expectedExposure = models::lgmOptionValue(model, zeta, flows);
        // E[V⁻] = E[V] − E[V⁺]; E[V] is today's value of the forward swap (deflated values are martingales).
        double forward = 0.0;
        for (const auto& f : flows) forward += f.amount * f.df;
        p.expectedNegativeExposure = forward - p.expectedExposure;
        // The swap's value is monotone in the state, so its q-quantile is its value at the q-quantile
        // (or (1 − q)-quantile) of the state.
        p.potentialExposure = std::max({models::deflatedValue(model, zeta, flows, zq * sd), models::deflatedValue(model, zeta, flows, -zq * sd), 0.0});
        double sum = 0.0, sumSq = 0.0;
        for (std::size_t k = 0; k < opt.paths; ++k) {
            const double v = std::max(models::deflatedValue(model, zeta, flows, sd * rng.gaussian()), 0.0);
            sum += v;
            sumSq += v * v;
        }
        const double n = static_cast<double>(opt.paths);
        p.expectedExposureMc = sum / n;
        p.mcStandardError = std::sqrt(std::max(sumSq / n - p.expectedExposureMc * p.expectedExposureMc, 0.0) / n);
        out.push_back(p);
    }
    return out;
}

} // namespace rl::pricers
