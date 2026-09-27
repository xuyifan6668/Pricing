#include "rl/instruments/conventions.hpp"

#include "rl/core/error.hpp"

namespace rl::instruments {

using core::BusinessDayConvention;
using core::DayCount;
using core::Index;
using core::StubRule;
using core::Tenor;

OisSwapConventions Conventions::usdSofrOis() {
    OisSwapConventions c;
    c.index = Index::sofr();
    c.spotLag = 2;
    c.spotCalendar = "USD-FEDWIRE";
    c.scheduleCalendar = "USD-SIFMA";
    c.convention = BusinessDayConvention::ModifiedFollowing;
    c.fixedFrequency = Tenor::years(1);
    c.fixedDayCount = DayCount::Act360;
    c.floatFrequency = Tenor::years(1);
    c.floatDayCount = DayCount::Act360;
    c.stub = StubRule::ShortFront;
    c.endOfMonth = false;
    c.paymentLag = 2;
    c.singlePeriodUpTo = Tenor::years(1);
    return c;
}

OisSwapConventions Conventions::usdFedFundsOis() {
    OisSwapConventions c = usdSofrOis();
    c.index = Index::fedFunds();
    return c;
}

BasisSwapConventions Conventions::usdFedFundsSofrBasis() {
    // SIMPLIFICATION: FF–SOFR basis swaps trade with the spread on the Fed Funds leg and both legs paying
    // on the same quarterly schedule. Dealers also quote an annual-pay version; the frequency is a data
    // field here so either can be built. The FF leg averages arithmetically (Index::fedFunds() carries
    // Compounding::AveragedDaily); the SOFR leg compounds.
    BasisSwapConventions c;
    c.spreadIndex = Index::fedFunds();
    c.flatIndex = Index::sofr();
    c.spotLag = 2;
    c.spotCalendar = "USD-FEDWIRE";
    c.scheduleCalendar = "USD-SIFMA";
    c.convention = BusinessDayConvention::ModifiedFollowing;
    c.frequency = Tenor::months(3);
    c.stub = StubRule::ShortFront;
    c.endOfMonth = false;
    c.paymentLag = 2;
    c.singlePeriodUpTo = Tenor::months(3);
    return c;
}

DepositConventions Conventions::usdSofrDeposit() {
    DepositConventions c;
    c.index = Index::sofr();
    c.spotLag = 2;
    c.spotCalendar = "USD-FEDWIRE";
    c.convention = BusinessDayConvention::ModifiedFollowing;
    c.endOfMonth = true;
    return c;
}

FutureConventions Conventions::usdSofr3mFuture() {
    FutureConventions c;
    c.index = Index::sofr();
    c.referenceTenor = Tenor::months(3);
    c.dayCount = DayCount::Act360;
    return c;
}

OisSwapConventions Conventions::oisSwap(std::string_view indexName) {
    if (indexName == "USD-SOFR") return usdSofrOis();
    if (indexName == "USD-FEDFUNDS") return usdFedFundsOis();
    throw InputError("no OIS swap conventions for index " + std::string(indexName));
}

BasisSwapConventions Conventions::basisSwap(std::string_view spreadIndex, std::string_view flatIndex) {
    if (spreadIndex == "USD-FEDFUNDS" && flatIndex == "USD-SOFR") return usdFedFundsSofrBasis();
    throw InputError("no basis swap conventions for " + std::string(spreadIndex) + " vs " + std::string(flatIndex));
}

DepositConventions Conventions::deposit(std::string_view indexName) {
    if (indexName == "USD-SOFR") return usdSofrDeposit();
    if (indexName == "USD-FEDFUNDS") {
        DepositConventions c = usdSofrDeposit();
        c.index = Index::fedFunds();
        return c;
    }
    throw InputError("no deposit conventions for index " + std::string(indexName));
}

FutureConventions Conventions::future(std::string_view indexName) {
    if (indexName == "USD-SOFR") return usdSofr3mFuture();
    throw InputError("no futures conventions for index " + std::string(indexName));
}

// --- IMM -----------------------------------------------------------------------------------------

core::Date thirdWednesday(int year, unsigned month) {
    const core::Date first(year, month, 1);
    const int wd = static_cast<int>(first.weekday()); // Sunday = 0 ... Wednesday = 3
    const int offset = (3 - wd + 7) % 7;
    return first + (offset + 14);
}

bool isQuarterlyImm(core::Date d) {
    const auto ymd = d.ymd();
    if (ymd.month % 3 != 0) return false;
    return d == thirdWednesday(ymd.year, ymd.month);
}

core::Date nextQuarterlyImm(core::Date d) {
    auto ymd = d.ymd();
    int y = ymd.year;
    unsigned m = ymd.month;
    // Candidate: the IMM date of the current quarter-month (round m up to a multiple of 3).
    m = ((m + 2) / 3) * 3;
    for (int k = 0; k < 8; ++k) {
        core::Date imm = thirdWednesday(y, m);
        if (imm > d) return imm;
        m += 3;
        if (m > 12) {
            m = 3;
            ++y;
        }
    }
    throw NumericalError("nextQuarterlyImm: no IMM date found"); // unreachable
}

core::Date nthQuarterlyImm(core::Date d, int n) {
    RL_REQUIRE(n >= 1, "nthQuarterlyImm: n must be ≥ 1");
    core::Date x = d;
    for (int i = 0; i < n; ++i) x = nextQuarterlyImm(x);
    return x;
}

char immMonthCode(unsigned month) {
    switch (month) {
    case 3: return 'H';
    case 6: return 'M';
    case 9: return 'U';
    case 12: return 'Z';
    default: throw InputError("not a quarterly IMM month: " + std::to_string(month));
    }
}

unsigned immMonthFromCode(char code) {
    switch (code) {
    case 'H': case 'h': return 3;
    case 'M': case 'm': return 6;
    case 'U': case 'u': return 9;
    case 'Z': case 'z': return 12;
    default: throw InputError(std::string("not a quarterly IMM month code: ") + code);
    }
}

core::Date addTenor(core::Date d, Tenor t) {
    switch (t.unit) {
    case core::TenorUnit::Days: return d.addDays(t.n);
    case core::TenorUnit::Weeks: return d.addDays(7 * t.n);
    case core::TenorUnit::Months: return d.addMonths(t.n);
    case core::TenorUnit::Years: return d.addYears(t.n);
    }
    return d;
}

core::Date spotDate(core::Date tradeDate, int lag, const core::Calendar& calendar) { return calendar.advance(tradeDate, lag); }

} // namespace rl::instruments
