#pragma once
// Definitions for every rl::linear template. Included only by src/linear.cpp (double) and
// src/instantiate_real.cpp (Real); consumers include the .hpp files.
//
// References: Andersen & Piterbarg, Interest Rate Modeling vol. 1 (2010), §5.1–5.3 (annuity, par rate),
// §6.2 (multi-curve projection vs discounting); ISDA 2021 Definitions §7.1 (compounded overnight rate,
// observation shift, lockout); CME SR3 contract specification (compounded reference rate, price = 100 − R).

#include "rl/core/calendar.hpp"
#include "rl/core/daycount.hpp"
#include "rl/core/error.hpp"
#include "rl/linear/annuity.hpp"
#include "rl/linear/dv01.hpp"
#include "rl/linear/forward.hpp"
#include "rl/linear/leg_pv.hpp"
#include "rl/linear/par_rate.hpp"

#include <algorithm>
#include <cmath>
#include <string>
#include <variant>

namespace rl::linear {

using std::log; // ADL finds aad::log for Real; this finds std::log for double

namespace detail {

/// Last business day on or before d (the fixing that applies on a non-business day is the previous one).
inline core::Date fixingDayFor(const core::Calendar& cal, core::Date d) {
    return cal.isBusinessDay(d) ? d : cal.adjust(d, core::BusinessDayConvention::Preceding);
}

} // namespace detail

// --- coupon rates -------------------------------------------------------------------------------------

template <math::Scalar T>
T overnightRate(const OvernightCoupon& c, const CurveSet<T>& curves, const FixingStore& fixings) {
    using core::Date;
    const core::Index& idx = c.index;
    const core::Calendar cal = core::Calendar::get(idx.fixingCalendar);
    const Date asOf = curves.referenceDate();
    const market::YieldCurve<T>& proj = curves.projection(idx.name);

    // Observation period (observation shift: shift both ends back by the lookback, weights from the
    // shifted dates). SIMPLIFICATION: the ISDA "lookback without observation shift" variant keeps the
    // accrual-date weights; the two differ by a fraction of a basis point on a quarterly coupon.
    const Date obsStart = c.lookbackDays > 0 ? cal.advance(c.accrualStart, -c.lookbackDays) : c.accrualStart;
    const Date obsEnd = c.lookbackDays > 0 ? cal.advance(c.accrualEnd, -c.lookbackDays) : c.accrualEnd;
    RL_REQUIRE(obsEnd > obsStart, "overnightRate: empty observation period for coupon paying " + c.payment.toIso());
    const double tauObs = core::yearFraction(idx.dayCount, obsStart, obsEnd);
    const bool compounded = idx.compounding != core::Compounding::AveragedDaily;

    // Lockout: on the last k business days the rate is frozen at the fixing of the k-th last day.
    const Date frozenDay = c.lockoutDays > 0 ? cal.advance(obsEnd, -c.lockoutDays) : obsEnd;

    // Realised part: business days d with obsStart ≤ d < min(obsEnd, asOf). Each contributes r_d · τ(d, next).
    T growth(1.0);
    T sum(0.0);
    const Date cutoff = std::min(obsEnd, asOf);
    Date d = obsStart;
    while (d < cutoff) {
        const Date next = cal.advance(d, 1);
        const double tau = core::yearFraction(idx.dayCount, d, next);
        Date fixDay = detail::fixingDayFor(cal, d);
        if (c.lockoutDays > 0 && fixDay >= frozenDay) fixDay = frozenDay;
        const double r = withContext("coupon " + c.accrualStart.toIso() + " to " + c.accrualEnd.toIso(),
                                     [&] { return fixings.require(idx.name, fixDay); });
        if (compounded) growth *= (1.0 + r * tau);
        else sum += r * tau;
        d = next;
    }
    // Projected remainder [d, obsEnd] off the projection curve. SIMPLIFICATION: lockout is not modelled
    // on the projected part (it changes nothing at first order under a smooth curve).
    if (d < obsEnd) {
        const T pStart = proj.df(d);
        const T pEnd = proj.df(obsEnd);
        if (compounded) growth *= pStart / pEnd;
        else sum += log(pStart / pEnd); // continuous-average ≈ arithmetic average of daily rates (error ~ r²τ/2 relative)
    }
    if (compounded) return (growth - 1.0) / tauObs;
    return sum / tauObs;
}

template <math::Scalar T>
T termRate(const TermCoupon& c, const CurveSet<T>& curves, const FixingStore& fixings) {
    const core::Date asOf = curves.referenceDate();
    if (auto f = fixings.fixing(c.index.name, c.fixingDate)) return T(*f);
    if (c.fixingDate < asOf)
        throw InputError("no " + c.index.name + " fixing for " + c.fixingDate.toIso() + " (coupon " + c.accrualStart.toIso() + " to " +
                         c.accrualEnd.toIso() + ")");
    const market::YieldCurve<T>& proj = curves.projection(c.index.name);
    return (proj.df(c.accrualStart) / proj.df(c.accrualEnd) - 1.0) / c.tau;
}

// --- cashflow PV -------------------------------------------------------------------------------------

template <math::Scalar T>
T pv(const FixedCoupon& c, const CurveSet<T>& curves) {
    if (c.payment <= curves.referenceDate()) return T(0.0);
    return c.notional * c.rate * c.tau * curves.discount().df(c.payment);
}

template <math::Scalar T>
T pv(const OvernightCoupon& c, const CurveSet<T>& curves, const FixingStore& fixings) {
    if (c.payment <= curves.referenceDate()) return T(0.0);
    const T rate = overnightRate(c, curves, fixings);
    return c.notional * (rate + c.spread) * c.tau * curves.discount().df(c.payment);
}

template <math::Scalar T>
T pv(const TermCoupon& c, const CurveSet<T>& curves, const FixingStore& fixings) {
    if (c.payment <= curves.referenceDate()) return T(0.0);
    const T rate = termRate(c, curves, fixings);
    return c.notional * (rate + c.spread) * c.tau * curves.discount().df(c.payment);
}

template <math::Scalar T>
T pv(const NotionalExchange& c, const CurveSet<T>& curves) {
    if (c.payment <= curves.referenceDate()) return T(0.0);
    return c.amount * curves.discount().df(c.payment);
}

template <math::Scalar T>
T pv(const Cashflow& cf, const CurveSet<T>& curves, const FixingStore& fixings) {
    struct Visitor {
        const CurveSet<T>& curves;
        const FixingStore& fixings;
        T operator()(const FixedCoupon& c) const { return pv(c, curves); }
        T operator()(const OvernightCoupon& c) const { return pv(c, curves, fixings); }
        T operator()(const TermCoupon& c) const { return pv(c, curves, fixings); }
        T operator()(const NotionalExchange& c) const { return pv(c, curves); }
    };
    return std::visit(Visitor{curves, fixings}, cf);
}

template <math::Scalar T>
T pv(const Leg& leg, const CurveSet<T>& curves, const FixingStore& fixings) {
    T total(0.0);
    for (const Cashflow& cf : leg.cashflows()) total += pv(cf, curves, fixings);
    return leg.sign() * total;
}

// --- trades ------------------------------------------------------------------------------------------

template <math::Scalar T>
T pv(const Swap& s, const CurveSet<T>& curves, const FixingStore& fixings) {
    return pv(s.fixedLeg(), curves, fixings) + pv(s.floatLeg(), curves, fixings);
}

template <math::Scalar T>
T pv(const BasisSwap& b, const CurveSet<T>& curves, const FixingStore& fixings) {
    T total(0.0);
    for (const Leg& leg : b.legs()) total += pv(leg, curves, fixings);
    return total;
}

template <math::Scalar T>
T impliedDepositRate(const Deposit& d, const CurveSet<T>& curves) {
    const market::YieldCurve<T>& proj = curves.projection(d.index().name);
    return (proj.df(d.start()) / proj.df(d.end()) - 1.0) / d.tau();
}

template <math::Scalar T>
T pv(const Deposit& d, const CurveSet<T>& curves) {
    if (d.end() <= curves.referenceDate()) return T(0.0);
    return d.notional() * d.tau() * (impliedDepositRate(d, curves) - d.rate()) * curves.discount().df(d.end());
}

template <math::Scalar T>
T impliedFraRate(const Fra& f, const CurveSet<T>& curves) {
    const market::YieldCurve<T>& proj = curves.projection(f.index().name);
    return (proj.df(f.start()) / proj.df(f.end()) - 1.0) / f.tau();
}

template <math::Scalar T>
T pv(const Fra& f, const CurveSet<T>& curves) {
    if (f.end() <= curves.referenceDate()) return T(0.0);
    return f.notional() * f.tau() * (impliedFraRate(f, curves) - f.rate()) * curves.discount().df(f.end());
}

template <math::Scalar T>
T impliedFutureForward(const Future& f, const CurveSet<T>& curves) {
    const market::YieldCurve<T>& proj = curves.projection(f.index().name);
    return (proj.df(f.start()) / proj.df(f.end()) - 1.0) / f.tau();
}

template <math::Scalar T>
T impliedFutureRate(const Future& f, const CurveSet<T>& curves) {
    return impliedFutureForward(f, curves) + f.convexityAdjustment();
}

template <math::Scalar T>
T impliedFuturePrice(const Future& f, const CurveSet<T>& curves) {
    return 100.0 * (1.0 - impliedFutureRate(f, curves));
}

// --- annuity, par rate, par spread ------------------------------------------------------------------

template <math::Scalar T>
T annuity(const Leg& leg, const CurveSet<T>& curves) {
    const core::Date asOf = curves.referenceDate();
    const market::YieldCurve<T>& disc = curves.discount();
    struct Visitor {
        core::Date asOf;
        const market::YieldCurve<T>& disc;
        T operator()(const FixedCoupon& c) const { return c.payment <= asOf ? T(0.0) : c.tau * c.notional * disc.df(c.payment); }
        T operator()(const OvernightCoupon& c) const { return c.payment <= asOf ? T(0.0) : c.tau * c.notional * disc.df(c.payment); }
        T operator()(const TermCoupon& c) const { return c.payment <= asOf ? T(0.0) : c.tau * c.notional * disc.df(c.payment); }
        T operator()(const NotionalExchange&) const { return T(0.0); }
    };
    T a(0.0);
    for (const Cashflow& cf : leg.cashflows()) a += std::visit(Visitor{asOf, disc}, cf);
    return a;
}

template <math::Scalar T>
T parRate(const Swap& s, const CurveSet<T>& curves, const FixingStore& fixings) {
    const T a = annuity(s.fixedLeg(), curves);
    RL_REQUIRE(math::value(a) > 0.0, "parRate: fixed-leg annuity is zero (all coupons settled?) for " + s.id());
    return -pv(s.floatLeg(), curves, fixings) / (s.fixedLeg().sign() * a);
}

namespace detail {
inline double uniformSpread(const Leg& leg, const std::string& id) {
    bool any = false;
    double spread = 0.0;
    for (const Cashflow& cf : leg.cashflows()) {
        double s = 0.0;
        bool has = false;
        if (const auto* o = std::get_if<OvernightCoupon>(&cf)) { s = o->spread; has = true; }
        else if (const auto* t = std::get_if<TermCoupon>(&cf)) { s = t->spread; has = true; }
        if (!has) continue;
        if (any) RL_REQUIRE(s == spread, "parSpread: coupons on the spread leg of " + id + " carry different spreads");
        spread = s;
        any = true;
    }
    RL_REQUIRE(any, "parSpread: no floating coupons on the spread leg of " + id);
    return spread;
}
} // namespace detail

template <math::Scalar T>
T parSpread(const BasisSwap& b, const CurveSet<T>& curves, const FixingStore& fixings, std::size_t spreadLeg) {
    RL_REQUIRE(spreadLeg < 2, "parSpread: spreadLeg must be 0 or 1");
    const Leg& leg = b.legs()[spreadLeg];
    const double current = detail::uniformSpread(leg, b.id());
    const T a = annuity(leg, curves);
    RL_REQUIRE(math::value(a) > 0.0, "parSpread: annuity of the spread leg is zero for " + b.id());
    return current - pv(b, curves, fixings) / (leg.sign() * a);
}

// --- forwards -----------------------------------------------------------------------------------------

template <math::Scalar T>
T simpleForward(const market::YieldCurve<T>& curve, double t1, double t2, double tau) {
    RL_REQUIRE(tau > 0.0 && t2 > t1, "simpleForward: need t2 > t1 and tau > 0");
    return (curve.df(t1) / curve.df(t2) - 1.0) / tau;
}

template <math::Scalar T>
T simpleForward(const market::YieldCurve<T>& curve, core::Date d1, core::Date d2, core::DayCount dc) {
    return simpleForward(curve, curve.time(d1), curve.time(d2), core::yearFraction(dc, d1, d2));
}

template <math::Scalar T>
T continuousForward(const market::YieldCurve<T>& curve, double t1, double t2) {
    RL_REQUIRE(t2 > t1, "continuousForward: need t2 > t1");
    return log(curve.df(t1) / curve.df(t2)) / (t2 - t1);
}

// --- DV01 ----------------------------------------------------------------------------------------------

template <math::Scalar T>
T fixedRateDv01(const Swap& s, const CurveSet<T>& curves) {
    return s.fixedLeg().sign() * annuity(s.fixedLeg(), curves) * 1e-4;
}

template <math::Scalar T>
CurveSet<T> parallelBump(const CurveSet<T>& curves, double rateBump) {
    CurveSet<T> out = curves;
    for (const std::string& name : curves.names()) {
        const market::YieldCurve<T>& c = curves.curve(name);
        std::vector<T> ys = c.nodeValues();
        const std::vector<double>& ts = c.nodeTimes();
        for (std::size_t i = 0; i < ys.size(); ++i) ys[i] += rateBump * ts[i];
        out = out.withCurve(name, c.withNodeValues(std::move(ys)));
    }
    return out;
}

template <math::Scalar T>
T parallelDv01(const Swap& s, const CurveSet<T>& curves, const FixingStore& fixings, double bumpBp) {
    return pv(s, parallelBump(curves, bumpBp * 1e-4), fixings) - pv(s, curves, fixings);
}

template <math::Scalar T>
T parallelDv01(const BasisSwap& b, const CurveSet<T>& curves, const FixingStore& fixings, double bumpBp) {
    return pv(b, parallelBump(curves, bumpBp * 1e-4), fixings) - pv(b, curves, fixings);
}

} // namespace rl::linear

// One line per template; add here when you add a function (a missing line is an unresolved external).
#define RL_LINEAR_INSTANTIATE(T)                                                                                          \
    namespace rl::linear {                                                                                                \
    template T overnightRate<T>(const OvernightCoupon&, const CurveSet<T>&, const FixingStore&);                          \
    template T termRate<T>(const TermCoupon&, const CurveSet<T>&, const FixingStore&);                                    \
    template T pv<T>(const FixedCoupon&, const CurveSet<T>&);                                                             \
    template T pv<T>(const OvernightCoupon&, const CurveSet<T>&, const FixingStore&);                                     \
    template T pv<T>(const TermCoupon&, const CurveSet<T>&, const FixingStore&);                                          \
    template T pv<T>(const NotionalExchange&, const CurveSet<T>&);                                                        \
    template T pv<T>(const Cashflow&, const CurveSet<T>&, const FixingStore&);                                            \
    template T pv<T>(const Leg&, const CurveSet<T>&, const FixingStore&);                                                 \
    template T pv<T>(const Swap&, const CurveSet<T>&, const FixingStore&);                                                \
    template T pv<T>(const BasisSwap&, const CurveSet<T>&, const FixingStore&);                                           \
    template T impliedDepositRate<T>(const Deposit&, const CurveSet<T>&);                                                 \
    template T pv<T>(const Deposit&, const CurveSet<T>&);                                                                 \
    template T impliedFraRate<T>(const Fra&, const CurveSet<T>&);                                                         \
    template T pv<T>(const Fra&, const CurveSet<T>&);                                                                     \
    template T impliedFutureForward<T>(const Future&, const CurveSet<T>&);                                                \
    template T impliedFutureRate<T>(const Future&, const CurveSet<T>&);                                                   \
    template T impliedFuturePrice<T>(const Future&, const CurveSet<T>&);                                                  \
    template T annuity<T>(const Leg&, const CurveSet<T>&);                                                                \
    template T parRate<T>(const Swap&, const CurveSet<T>&, const FixingStore&);                                           \
    template T parSpread<T>(const BasisSwap&, const CurveSet<T>&, const FixingStore&, std::size_t);                       \
    template T simpleForward<T>(const market::YieldCurve<T>&, double, double, double);                                    \
    template T simpleForward<T>(const market::YieldCurve<T>&, core::Date, core::Date, core::DayCount);                    \
    template T continuousForward<T>(const market::YieldCurve<T>&, double, double);                                        \
    template T fixedRateDv01<T>(const Swap&, const CurveSet<T>&);                                                         \
    template CurveSet<T> parallelBump<T>(const CurveSet<T>&, double);                                                     \
    template T parallelDv01<T>(const Swap&, const CurveSet<T>&, const FixingStore&, double);                              \
    template T parallelDv01<T>(const BasisSwap&, const CurveSet<T>&, const FixingStore&, double);                         \
    }
