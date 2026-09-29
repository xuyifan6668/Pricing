#include "rl/instruments/builder.hpp"

#include "rl/core/daycount.hpp"
#include "rl/core/error.hpp"

#include <cmath>

namespace rl::instruments {

using core::Calendar;
using core::Date;
using core::Schedule;
using core::ScheduleSpec;
using core::Tenor;

// --- leg factories ---------------------------------------------------------------------------------

Leg makeFixedLeg(const Schedule& schedule, core::DayCount dayCount, double rate, double notional, Side side, core::Currency currency) {
    RL_REQUIRE(!schedule.empty(), "makeFixedLeg: empty schedule");
    RL_REQUIRE(std::isfinite(rate), "makeFixedLeg: non-finite rate");
    RL_REQUIRE(notional > 0.0, "makeFixedLeg: notional must be positive (the leg's Side carries the sign)");
    std::vector<Cashflow> cfs;
    cfs.reserve(schedule.size());
    for (const auto& p : schedule) {
        FixedCoupon c;
        c.accrualStart = p.start;
        c.accrualEnd = p.end;
        c.payment = p.payment;
        c.tau = core::yearFraction(dayCount, p.start, p.end);
        c.notional = notional;
        c.rate = rate;
        cfs.emplace_back(std::move(c));
    }
    return Leg(std::move(cfs), side, currency);
}

Leg makeOvernightLeg(const Schedule& schedule, const core::Index& index, double spread, double notional, Side side, core::Currency currency,
                     int lookbackDays, int lockoutDays) {
    RL_REQUIRE(!schedule.empty(), "makeOvernightLeg: empty schedule");
    RL_REQUIRE(index.isOvernight(), "makeOvernightLeg: " + index.name + " is not an overnight index");
    RL_REQUIRE(notional > 0.0, "makeOvernightLeg: notional must be positive");
    RL_REQUIRE(lookbackDays >= 0 && lockoutDays >= 0, "makeOvernightLeg: lookback/lockout must be non-negative");
    std::vector<Cashflow> cfs;
    cfs.reserve(schedule.size());
    for (const auto& p : schedule) {
        OvernightCoupon c;
        c.accrualStart = p.start;
        c.accrualEnd = p.end;
        c.payment = p.payment;
        c.tau = core::yearFraction(index.dayCount, p.start, p.end);
        c.notional = notional;
        c.index = index;
        c.spread = spread;
        c.lookbackDays = lookbackDays;
        c.lockoutDays = lockoutDays;
        cfs.emplace_back(std::move(c));
    }
    return Leg(std::move(cfs), side, currency);
}

Leg makeTermLeg(const Schedule& schedule, const core::Index& index, double spread, double notional, Side side, core::Currency currency) {
    RL_REQUIRE(!schedule.empty(), "makeTermLeg: empty schedule");
    RL_REQUIRE(!index.isOvernight(), "makeTermLeg: " + index.name + " is an overnight index; use makeOvernightLeg");
    RL_REQUIRE(notional > 0.0, "makeTermLeg: notional must be positive");
    const Calendar fixingCal = Calendar::get(index.fixingCalendar);
    std::vector<Cashflow> cfs;
    cfs.reserve(schedule.size());
    for (const auto& p : schedule) {
        TermCoupon c;
        c.fixingDate = fixingCal.advance(p.start, -index.fixingLag);
        c.accrualStart = p.start;
        c.accrualEnd = p.end;
        c.payment = p.payment;
        c.tau = core::yearFraction(index.dayCount, p.start, p.end);
        c.notional = notional;
        c.index = index;
        c.spread = spread;
        cfs.emplace_back(std::move(c));
    }
    return Leg(std::move(cfs), side, currency);
}

double tenorDaysApprox(Tenor t) noexcept {
    switch (t.unit) {
    case core::TenorUnit::Days: return t.n;
    case core::TenorUnit::Weeks: return 7.0 * t.n;
    case core::TenorUnit::Months: return 30.4375 * t.n;
    case core::TenorUnit::Years: return 365.25 * t.n;
    }
    return 0.0;
}

namespace {
/// A schedule from effective to termination: one period if the whole span is ≤ `singlePeriodUpTo`,
/// otherwise regular periods of `frequency` with the given stub rule.
Schedule makeSchedule(Date effective, Date termination, Tenor frequency, Tenor singlePeriodUpTo, const std::string& calendarName,
                      core::BusinessDayConvention bdc, core::StubRule stub, bool eom, int paymentLag) {
    RL_REQUIRE(termination > effective, "schedule: termination must follow effective");
    ScheduleSpec spec;
    spec.effective = effective;
    spec.termination = termination;
    spec.calendar = Calendar::get(calendarName);
    spec.convention = bdc;
    spec.terminationConvention = bdc;
    spec.stub = stub;
    spec.endOfMonth = eom;
    spec.paymentLag = paymentLag;
    const double spanDays = static_cast<double>(termination - effective);
    if (spanDays <= tenorDaysApprox(singlePeriodUpTo) + 0.5) {
        spec.frequency = Tenor::days(static_cast<int>(spanDays));
        return Schedule::build(spec, {effective, termination});
    }
    spec.frequency = frequency;
    return Schedule::generate(spec);
}
} // namespace

// --- OisSwapBuilder --------------------------------------------------------------------------------

OisSwapBuilder::OisSwapBuilder(OisSwapConventions conventions, Date tradeDate) : conv_(std::move(conventions)), tradeDate_(tradeDate) {}

Date OisSwapBuilder::effectiveDate() const {
    if (effective_) return *effective_;
    return spotDate(tradeDate_, conv_.spotLag, Calendar::get(conv_.spotCalendar));
}

Date OisSwapBuilder::terminationDate() const {
    if (maturity_) return *maturity_;
    RL_REQUIRE(tenor_.has_value(), "OisSwapBuilder: set tenor() or maturity()");
    return addTenor(effectiveDate(), *tenor_);
}

Schedule OisSwapBuilder::schedule(Tenor frequency) const {
    return makeSchedule(effectiveDate(), terminationDate(), frequency, conv_.singlePeriodUpTo, conv_.scheduleCalendar, conv_.convention,
                        conv_.stub, conv_.endOfMonth, conv_.paymentLag);
}

Schedule OisSwapBuilder::fixedSchedule() const { return schedule(conv_.fixedFrequency); }
Schedule OisSwapBuilder::floatSchedule() const { return schedule(conv_.floatFrequency); }

Swap OisSwapBuilder::build() const {
    RL_REQUIRE(!id_.empty(), "OisSwapBuilder: id is required");
    RL_REQUIRE(fixedRate_.has_value(), "OisSwapBuilder: fixedRate is required");
    return withContext("build " + id_, [&] {
        const core::Currency ccy = conv_.index.currency;
        Leg fixed = makeFixedLeg(fixedSchedule(), conv_.fixedDayCount, *fixedRate_, notional_, payFixed_ ? Side::Pay : Side::Receive, ccy);
        Leg floating = makeOvernightLeg(floatSchedule(), conv_.index, spread_, notional_, payFixed_ ? Side::Receive : Side::Pay, ccy, lookback_, lockout_);
        return Swap(id_, std::move(fixed), std::move(floating));
    });
}

// --- BasisSwapBuilder ------------------------------------------------------------------------------

BasisSwapBuilder::BasisSwapBuilder(BasisSwapConventions conventions, Date tradeDate) : conv_(std::move(conventions)), tradeDate_(tradeDate) {}

Date BasisSwapBuilder::effectiveDate() const {
    if (effective_) return *effective_;
    return spotDate(tradeDate_, conv_.spotLag, Calendar::get(conv_.spotCalendar));
}

Date BasisSwapBuilder::terminationDate() const {
    if (maturity_) return *maturity_;
    RL_REQUIRE(tenor_.has_value(), "BasisSwapBuilder: set tenor() or maturity()");
    return addTenor(effectiveDate(), *tenor_);
}

Schedule BasisSwapBuilder::schedule() const {
    return makeSchedule(effectiveDate(), terminationDate(), conv_.frequency, conv_.singlePeriodUpTo, conv_.scheduleCalendar, conv_.convention,
                        conv_.stub, conv_.endOfMonth, conv_.paymentLag);
}

BasisSwap BasisSwapBuilder::build() const {
    RL_REQUIRE(!id_.empty(), "BasisSwapBuilder: id is required");
    return withContext("build " + id_, [&] {
        const core::Currency ccy = conv_.spreadIndex.currency;
        const Schedule s = schedule();
        Leg spreadLeg = makeOvernightLeg(s, conv_.spreadIndex, spread_, notional_, paySpreadLeg_ ? Side::Pay : Side::Receive, ccy);
        Leg flatLeg = makeOvernightLeg(s, conv_.flatIndex, 0.0, notional_, paySpreadLeg_ ? Side::Receive : Side::Pay, ccy);
        return BasisSwap(id_, std::move(spreadLeg), std::move(flatLeg));
    });
}

// --- DepositBuilder --------------------------------------------------------------------------------

DepositBuilder::DepositBuilder(DepositConventions conventions, Date tradeDate) : conv_(std::move(conventions)), tradeDate_(tradeDate) {}

Deposit DepositBuilder::build() const {
    RL_REQUIRE(!id_.empty(), "DepositBuilder: id is required");
    RL_REQUIRE(rate_.has_value(), "DepositBuilder: rate is required");
    return withContext("build " + id_, [&] {
        const Calendar cal = Calendar::get(conv_.spotCalendar);
        const Date start = start_ ? *start_ : spotDate(tradeDate_, conv_.spotLag, cal);
        Date end;
        if (end_) end = *end_;
        else {
            RL_REQUIRE(tenor_.has_value(), "DepositBuilder: set tenor() or end()");
            end = cal.advance(start, *tenor_, conv_.convention, conv_.endOfMonth);
        }
        RL_REQUIRE(end > start, "DepositBuilder: end must follow start");
        const double tau = core::yearFraction(conv_.index.dayCount, start, end);
        return Deposit(id_, conv_.index, start, end, tau, notional_, *rate_);
    });
}

// --- FraBuilder ------------------------------------------------------------------------------------

FraBuilder::FraBuilder(core::Index index) : index_(std::move(index)) {}

Fra FraBuilder::build() const {
    RL_REQUIRE(!id_.empty(), "FraBuilder: id is required");
    RL_REQUIRE(start_ && end_, "FraBuilder: start and end are required");
    RL_REQUIRE(rate_.has_value(), "FraBuilder: rate is required");
    RL_REQUIRE(*end_ > *start_, "FraBuilder: end must follow start");
    return withContext("build " + id_, [&] {
        const Calendar fixingCal = Calendar::get(index_.fixingCalendar);
        const Date fixing = fixingCal.advance(*start_, -index_.fixingLag);
        const double tau = core::yearFraction(index_.dayCount, *start_, *end_);
        return Fra(id_, index_, fixing, *start_, *end_, tau, notional_, *rate_);
    });
}

// --- FutureBuilder ---------------------------------------------------------------------------------

FutureBuilder::FutureBuilder(FutureConventions conventions) : conv_(std::move(conventions)) {}

FutureBuilder& FutureBuilder::contract(int year, unsigned month) {
    RL_REQUIRE(month % 3 == 0, "FutureBuilder::contract: month must be a quarterly IMM month (3, 6, 9, 12)");
    start_ = thirdWednesday(year, month);
    end_ = nextQuarterlyImm(*start_);
    return *this;
}

Future FutureBuilder::build() const {
    RL_REQUIRE(!id_.empty(), "FutureBuilder: id is required");
    RL_REQUIRE(start_.has_value(), "FutureBuilder: set contract() or start()");
    RL_REQUIRE(price_.has_value(), "FutureBuilder: price is required");
    return withContext("build " + id_, [&] {
        const Date start = *start_;
        const Date end = end_ ? *end_ : nextQuarterlyImm(start);
        RL_REQUIRE(end > start, "FutureBuilder: end must follow start");
        const double tau = core::yearFraction(conv_.dayCount, start, end);
        return Future(id_, conv_.index, start, end, tau, *price_, convexityBp_ * 1e-4);
    });
}

} // namespace rl::instruments
