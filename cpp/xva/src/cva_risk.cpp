#include "rl/xva/cva_risk.hpp"

#include "rl/aad/real.hpp"
#include "rl/aad/tape.hpp"
#include "rl/bootstrap/jacobian.hpp"
#include "rl/core/error.hpp"
#include "rl/math/normal.hpp"
#include "rl/models/lgm_swaption.hpp"
#include "rl/xva/exposure_engine.hpp"
#include "rl/xva/xva.hpp"

#include <cmath>
#include <cstdint>
#include <variant>

namespace rl::xva {

namespace {

using instruments::FixedCoupon;
using instruments::OvernightCoupon;

/// A swap as bonds with the coefficients c_k = amount × discount factor as functions of the curve:
/// one bond per fixed coupon, two per floating coupon, in leg order (the same set as
/// models::swapAsBonds, which the tests check). `first` is the index of the first fixed coupon to
/// include; floating coupons are matched by accrual start.
template <class T>
std::vector<T> coefficients(const instruments::Swap& swap, const market::CurveSet<T>& curves, core::Date from) {
    const market::YieldCurve<T>& disc = curves.discount();
    std::vector<T> c;
    for (const auto& cf : swap.fixedLeg().cashflows()) {
        const auto& k = std::get<FixedCoupon>(cf);
        if (k.accrualStart < from) continue;
        c.push_back(swap.fixedLeg().sign() * k.notional * k.rate * k.tau * disc.df(k.payment));
    }
    for (const auto& cf : swap.floatLeg().cashflows()) {
        const auto& k = std::get<OvernightCoupon>(cf);
        if (k.accrualStart < from) continue;
        const T n = swap.floatLeg().sign() * k.notional * disc.df(k.payment) / disc.df(k.accrualEnd);
        c.push_back(n * disc.df(k.accrualStart));
        c.push_back(-1.0 * n * disc.df(k.accrualEnd));
    }
    return c;
}

/// Times of the bonds, in the same order as `coefficients`.
std::vector<double> bondTimes(const instruments::Swap& swap, const market::CurveSet<double>& curves, core::Date from) {
    const market::YieldCurve<double>& disc = curves.discount();
    std::vector<double> t;
    for (const auto& cf : swap.fixedLeg().cashflows()) {
        const auto& k = std::get<FixedCoupon>(cf);
        if (k.accrualStart >= from) t.push_back(disc.time(k.payment));
    }
    for (const auto& cf : swap.floatLeg().cashflows()) {
        const auto& k = std::get<OvernightCoupon>(cf);
        if (k.accrualStart < from) continue;
        t.push_back(disc.time(k.accrualStart));
        t.push_back(disc.time(k.accrualEnd));
    }
    return t;
}

double onebp(const bootstrap::BootstrapInstrument& bi) { return bi.type == market::QuoteType::Future ? -0.01 : 1e-4; }

} // namespace

double swapCvaClosedForm(const models::Lgm1F& model, const instruments::Swap& swap, const market::CurveSet<double>& curves, const HazardCurve& hazard,
                         double recovery) {
    const SwapExposure e = swapExposureClosedForm(model, swap, curves);
    return creditAdjustment(e.times, e.epe, hazard, recovery);
}

std::vector<pricers::RiskRow> cs01Ladder(const std::vector<CdsQuote>& quotes, double recovery, const std::function<double(double)>& discount,
                                         const std::function<double(const HazardCurve&)>& adjustment, double bumpBp) {
    const double base = adjustment(HazardCurve::bootstrap(quotes, recovery, discount));
    std::vector<pricers::RiskRow> out;
    for (std::size_t i = 0; i < quotes.size(); ++i) {
        std::vector<CdsQuote> bumped = quotes;
        bumped[i].spread += bumpBp * 1e-4;
        const double years = quotes[i].maturityYears;
        const std::string id = "CDS-" + (std::fabs(years - std::round(years)) < 1e-9 ? std::to_string(static_cast<int>(std::lround(years))) : std::to_string(years)) + "Y";
        out.push_back({id, (adjustment(HazardCurve::bootstrap(bumped, recovery, discount)) - base) / bumpBp});
    }
    return out;
}

std::vector<pricers::RiskRow> swapCvaDeltaAdjoint(const models::UsdMarket& market, const models::Lgm1F& model, const instruments::Swap& swap,
                                                  const HazardCurve& hazard, double recovery) {
    const market::CurveSet<double>& curves = market.curves;
    const market::YieldCurve<double>& disc = curves.discount();

    // --- Step 1, in doubles: the weight of every bond of every swaption in the strip ------------------
    struct Date {
        core::Date from;
        double creditWeight = 0.0;       ///< (1 − R) × the trapezoid weight of this date's default probabilities
        std::vector<double> bondWeight;  ///< ∂EPE/∂c_k, in the order of `coefficients`
    };
    const std::vector<double> times = periodStartTimes(swap, curves);
    std::vector<Date> dates(times.size());
    for (std::size_t i = 0; i < times.size(); ++i) {
        const double before = i > 0 ? hazard.defaultProbability(times[i - 1], times[i]) : 0.0;
        const double after = i + 1 < times.size() ? hazard.defaultProbability(times[i], times[i + 1]) : 0.0;
        dates[i].creditWeight = (1.0 - recovery) * 0.5 * (before + after);
    }
    {
        std::size_t i = 0;
        dates[0].from = curves.referenceDate();
        for (const auto& cf : swap.fixedLeg().cashflows()) {
            const core::Date start = std::get<FixedCoupon>(cf).accrualStart;
            if (disc.time(start) > times[i] + 1e-9) dates[++i].from = start;
        }
        RL_REQUIRE(i + 1 == times.size(), "swapCvaDeltaAdjoint: the swap's period starts do not match its exposure dates");
    }
    for (std::size_t i = 0; i < dates.size(); ++i) {
        const std::vector<double> c = coefficients<double>(swap, curves, dates[i].from), t = bondTimes(swap, curves, dates[i].from);
        std::vector<models::BondFlow> flows;
        for (std::size_t k = 0; k < c.size(); ++k) flows.push_back({t[k], c[k] / disc.df(t[k]), disc.df(t[k])});
        dates[i].bondWeight.assign(c.size(), 0.0);
        if (i == 0) {
            // Today the exposure is max(V, 0): the derivative is the indicator.
            double v = 0.0;
            for (double ck : c) v += ck;
            if (v > 0.0) dates[i].bondWeight.assign(c.size(), 1.0);
            continue;
        }
        const double zeta = model.zeta(times[i]), sd = std::sqrt(zeta);
        const models::ExerciseBoundary b = models::exerciseBoundary(model, zeta, flows);
        for (std::size_t k = 0; k < c.size(); ++k) {
            if (!b.exists) dates[i].bondWeight[k] = b.exerciseAbove ? 1.0 : 0.0;
            else dates[i].bondWeight[k] = math::normalCdf((b.exerciseAbove ? -1.0 : 1.0) * (b.x + model.H(t[k]) * zeta) / sd);
        }
    }

    // --- Steps 2 and 3: one linear combination on the tape, one sweep, one transposed solve ------------
    const bootstrap::BootstrapProblem& problem = market.problem;
    const std::size_t n = problem.parameterCount();
    aad::Tape& tape = aad::Tape::current();
    const aad::Tape::Position start = tape.mark();
    math::Vector dCdz(static_cast<Eigen::Index>(n));
    {
        std::vector<aad::Real> z;
        z.reserve(n);
        for (std::size_t j = 0; j < n; ++j) z.push_back(aad::Real::input(market.z[static_cast<Eigen::Index>(j)]));
        const market::CurveSet<aad::Real> active = problem.curves(z);
        aad::Real cva = 0.0;
        for (const Date& d : dates) {
            const std::vector<aad::Real> c = coefficients<aad::Real>(swap, active, d.from);
            for (std::size_t k = 0; k < c.size(); ++k) cva += (d.creditWeight * d.bondWeight[k]) * c[k];
        }
        if (cva.isActive()) {
            tape.clearAdjointsFrom(static_cast<std::uint32_t>(start.nodes));
            tape.seedAdjoint(cva.node(), 1.0);
            tape.propagate(cva.node(), static_cast<std::uint32_t>(start.nodes));
        }
        for (std::size_t j = 0; j < n; ++j) dCdz[static_cast<Eigen::Index>(j)] = cva.isActive() ? tape.adjoint(z[j].node()) : 0.0;
    }
    tape.rewindTo(start);
    const math::Vector dCdq = bootstrap::parRisk(problem, market.z, dCdz);
    std::vector<pricers::RiskRow> out;
    for (std::size_t i = 0; i < problem.size(); ++i) out.push_back({problem.instrumentIds()[i], dCdq[static_cast<Eigen::Index>(i)] * onebp(problem.instruments()[i])});
    return out;
}

} // namespace rl::xva
