#pragma once
// Shared fixtures for the linear tests: flat curves (exact under linear-on-log-df), a fixing store with
// a constant rate, and coupon factories.
#include "rl/core/calendar.hpp"
#include "rl/core/date.hpp"
#include "rl/core/index.hpp"
#include "rl/instruments/cashflow.hpp"
#include "rl/instruments/leg.hpp"
#include "rl/market/curve_set.hpp"
#include "rl/market/fixings.hpp"
#include "rl/market/yield_curve.hpp"

#include <cmath>
#include <string>
#include <vector>

namespace rl::linear::test {

inline const core::Date kAsOf(2026, 9, 4); // Friday

template <class T>
market::YieldCurve<T> flatCurve(const std::string& name, double r, core::Date ref = kAsOf) {
    std::vector<double> ts{0.5, 1.0, 2.0, 5.0, 10.0, 30.0};
    std::vector<T> ys;
    for (double t : ts) ys.push_back(T(r * t));
    return market::YieldCurve<T>(ref, name, ts, ys, market::CurveInterpolation::LinearLogDf);
}

/// SOFR curve (discount + USD-SOFR projection) at rSofr, FF projection at rFf.
template <class T>
market::CurveSet<T> flatSet(double rSofr, double rFf, core::Date ref = kAsOf) {
    return market::CurveSet<T>({flatCurve<T>("USD-SOFR-OIS", rSofr, ref), flatCurve<T>("USD-FF", rFf, ref)}, "USD-SOFR-OIS",
                               {{"USD-SOFR", "USD-SOFR-OIS"}, {"USD-FEDFUNDS", "USD-FF"}});
}

/// Constant fixings for every business day in [from, to] on the index's fixing calendar.
inline market::FixingStore constantFixings(const core::Index& idx, core::Date from, core::Date to, double r) {
    market::FixingStore fs;
    const core::Calendar cal = core::Calendar::get(idx.fixingCalendar);
    for (core::Date d = from; d <= to; ++d)
        if (cal.isBusinessDay(d)) fs.add(idx.name, d, r);
    return fs;
}

inline instruments::FixedCoupon fixedCoupon(core::Date s, core::Date e, core::Date pay, double tau, double notional, double rate) {
    instruments::FixedCoupon c;
    c.accrualStart = s;
    c.accrualEnd = e;
    c.payment = pay;
    c.tau = tau;
    c.notional = notional;
    c.rate = rate;
    return c;
}

inline instruments::OvernightCoupon oisCoupon(core::Date s, core::Date e, core::Date pay, double notional, const core::Index& idx = core::Index::sofr(),
                                              double spread = 0.0, int lookback = 0, int lockout = 0) {
    instruments::OvernightCoupon c;
    c.accrualStart = s;
    c.accrualEnd = e;
    c.payment = pay;
    c.tau = core::yearFraction(idx.dayCount, s, e);
    c.notional = notional;
    c.index = idx;
    c.spread = spread;
    c.lookbackDays = lookback;
    c.lockoutDays = lockout;
    return c;
}

} // namespace rl::linear::test
