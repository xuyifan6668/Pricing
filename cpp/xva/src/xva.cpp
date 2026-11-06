#include "rl/xva/xva.hpp"

#include "rl/core/error.hpp"
#include "rl/models/flows.hpp"
#include "rl/models/lgm_swaption.hpp"
#include "rl/xva/exposure_engine.hpp"

#include <algorithm>
#include <cmath>

namespace rl::xva {

namespace {
void check(const std::vector<double>& times, const std::vector<double>& exposure, const char* who) {
    RL_REQUIRE(!times.empty() && times.size() == exposure.size(), std::string(who) + ": times and exposure must match and be non-empty");
    for (std::size_t i = 1; i < times.size(); ++i) RL_REQUIRE(times[i] > times[i - 1], std::string(who) + ": times must be strictly increasing");
}
} // namespace

double creditAdjustment(const std::vector<double>& times, const std::vector<double>& exposure, const HazardCurve& hazard, double recovery) {
    check(times, exposure, "creditAdjustment");
    RL_REQUIRE(recovery >= 0.0 && recovery <= 1.0, "creditAdjustment: recovery must be in [0, 1]");
    double sum = 0.0;
    for (std::size_t i = 1; i < times.size(); ++i) sum += 0.5 * (exposure[i - 1] + exposure[i]) * hazard.defaultProbability(times[i - 1], times[i]);
    return (1.0 - recovery) * sum;
}

double fundingAdjustment(const std::vector<double>& times, const std::vector<double>& exposure, double spread) {
    check(times, exposure, "fundingAdjustment");
    double sum = 0.0;
    for (std::size_t i = 1; i < times.size(); ++i) sum += 0.5 * (exposure[i - 1] + exposure[i]) * (times[i] - times[i - 1]);
    return spread * sum;
}

XvaResult computeXva(const std::vector<double>& times, const std::vector<double>& epe, const std::vector<double>& ene, const XvaInputs& in) {
    std::vector<double> owed(ene.size());
    for (std::size_t i = 0; i < ene.size(); ++i) owed[i] = -ene[i];
    XvaResult r;
    r.cva = creditAdjustment(times, epe, in.counterparty, in.counterpartyRecovery);
    r.dva = creditAdjustment(times, owed, in.own, in.ownRecovery);
    r.fca = fundingAdjustment(times, epe, in.fundingSpread);
    r.fba = fundingAdjustment(times, owed, in.fundingSpread);
    return r;
}

SwapExposure swapExposureClosedForm(const models::Lgm1F& model, const instruments::Swap& swap, const market::CurveSet<double>& curves) {
    SwapExposure out;
    out.times = periodStartTimes(swap, curves);
    const market::YieldCurve<double>& disc = curves.discount();
    auto forwardValue = [](const std::vector<models::BondFlow>& flows) {
        double v = 0.0;
        for (const auto& f : flows) v += f.amount * f.df;
        return v;
    };
    const double today = forwardValue(models::swapAsBonds(swap, curves, curves.referenceDate()));
    out.epe.push_back(std::max(today, 0.0));
    out.ene.push_back(std::min(today, 0.0));
    double last = 0.0;
    for (const auto& cf : swap.fixedLeg().cashflows()) {
        const core::Date start = std::get<instruments::FixedCoupon>(cf).accrualStart;
        const double t = disc.time(start);
        if (t <= last + 1e-9) continue;
        last = t;
        const std::vector<models::BondFlow> flows = models::swapAsBonds(swap, curves, start);
        // E[V⁺] is the swaption; E[V⁻] = E[V] − E[V⁺], and E[V] is today's value of the remaining swap.
        const double positive = models::lgmOptionValue(model, model.zeta(t), flows);
        out.epe.push_back(positive);
        out.ene.push_back(forwardValue(flows) - positive);
    }
    RL_REQUIRE(out.epe.size() == out.times.size(), "swapExposureClosedForm: the swap must start on or after today");
    return out;
}

} // namespace rl::xva
