#include "rl/models/flows.hpp"

#include "rl/core/calendar.hpp"
#include "rl/core/daycount.hpp"
#include "rl/core/error.hpp"

#include <algorithm>
#include <variant>

namespace rl::models {

using instruments::FixedCoupon;
using instruments::OvernightCoupon;

double ExerciseRight::forwardValue() const noexcept {
    double v = 0.0;
    for (const BondFlow& f : flows) v += f.amount * f.df;
    return v;
}

double ExerciseRight::annuityValue() const noexcept {
    double a = 0.0;
    for (const BondFlow& f : annuity) a += f.amount * f.df;
    return a;
}

double ExerciseRight::forwardRate() const noexcept { return strike + (payer ? 1.0 : -1.0) * forwardValue() / annuityValue(); }

std::vector<BondFlow> swapAsBonds(const instruments::Swap& swap, const market::CurveSet<double>& curves, core::Date from) {
    const market::YieldCurve<double>& disc = curves.discount();
    std::vector<BondFlow> out;
    auto add = [&](core::Date d, double amount) {
        const double t = disc.time(d);
        out.push_back({t, amount, disc.df(t)});
    };

    const double fixedSign = swap.fixedLeg().sign();
    for (const auto& cf : swap.fixedLeg().cashflows()) {
        const auto* c = std::get_if<FixedCoupon>(&cf);
        RL_REQUIRE(c != nullptr, "swapAsBonds: the fixed leg of " + swap.id() + " must hold fixed coupons only");
        if (c->accrualStart < from) continue;
        add(c->payment, fixedSign * c->notional * c->rate * c->tau);
    }
    const double floatSign = swap.floatLeg().sign();
    for (const auto& cf : swap.floatLeg().cashflows()) {
        const auto* c = std::get_if<OvernightCoupon>(&cf);
        RL_REQUIRE(c != nullptr, "swapAsBonds: the floating leg of " + swap.id() + " must hold overnight coupons only");
        if (c->accrualStart < from) continue;
        RL_REQUIRE(c->spread == 0.0, "swapAsBonds: a spread on the floating leg is not supported");
        RL_REQUIRE(&curves.projection(c->index.name) == &disc, "swapAsBonds: " + c->index.name + " is not projected off the discount curve (single-curve model)");
        RL_REQUIRE(c->accrualStart >= curves.referenceDate(), "swapAsBonds: a floating coupon of " + swap.id() + " has already started accruing");
        const double delay = disc.df(c->payment) / disc.df(c->accrualEnd);
        add(c->accrualStart, floatSign * c->notional * delay);
        add(c->accrualEnd, -floatSign * c->notional * delay);
    }
    std::sort(out.begin(), out.end(), [](const BondFlow& a, const BondFlow& b) { return a.time < b.time; });
    return out;
}

ExerciseRight makeExerciseRight(const instruments::Swap& swap, const market::CurveSet<double>& curves, core::Date exerciseDate) {
    const market::YieldCurve<double>& disc = curves.discount();
    ExerciseRight r;
    r.date = exerciseDate;
    r.time = disc.time(exerciseDate);
    RL_REQUIRE(r.time > 0.0, "makeExerciseRight: the exercise date must be after the valuation date");
    r.flows = swapAsBonds(swap, curves, exerciseDate);
    r.payer = swap.payerOfFixed();
    core::Date first{}, last{};
    bool any = false;
    for (const auto& cf : swap.fixedLeg().cashflows()) {
        const auto& c = std::get<FixedCoupon>(cf);
        if (c.accrualStart < exerciseDate) continue;
        const double t = disc.time(c.payment);
        r.annuity.push_back({t, c.notional * c.tau, disc.df(t)});
        r.strike = c.rate;
        if (!any) first = c.accrualStart;
        last = c.accrualEnd;
        any = true;
    }
    RL_REQUIRE(any, "makeExerciseRight: no coupons of " + swap.id() + " start on or after the exercise date");
    r.tenorYears = core::yearFraction(core::DayCount::Act365Fixed, first, last);
    r.label = exerciseDate.toIso();
    return r;
}

std::vector<ExerciseRight> makeExerciseRights(const instruments::BermudanSwaption& bermudan, const market::CurveSet<double>& curves) {
    std::vector<ExerciseRight> out;
    for (core::Date d : bermudan.exerciseDates()) out.push_back(makeExerciseRight(bermudan.underlying(), curves, d));
    return out;
}

std::vector<core::Date> coterminalExerciseDates(const instruments::Swap& swap, double firstCallYears, int noticeDays, const std::string& calendar) {
    const core::Calendar cal = core::Calendar::get(calendar);
    const auto& cfs = swap.fixedLeg().cashflows();
    RL_REQUIRE(!cfs.empty(), "coterminalExerciseDates: empty fixed leg");
    const core::Date start = std::get<FixedCoupon>(cfs.front()).accrualStart;
    std::vector<core::Date> out;
    for (const auto& cf : cfs) {
        const core::Date s = std::get<FixedCoupon>(cf).accrualStart;
        if (core::yearFraction(core::DayCount::Act365Fixed, start, s) < firstCallYears - 0.02) continue; // 0.02y ≈ a week of date-roll slack
        out.push_back(noticeDays > 0 ? cal.advance(s, -noticeDays) : s);
    }
    return out;
}

} // namespace rl::models
