#pragma once
// Market-standard trade conventions per (currency, index): the data a builder needs to turn "10y USD SOFR
// OIS at 3.45%" into dated cashflows. Pure data. Where a convention is a desk choice rather than
// a market fact it is marked and lives here so that changing it is a data edit, not a code edit.
//
// Sources: ISDA 2021 Definitions (overnight compounding, payment delay), CME SR3 contract specs (IMM
// reference periods), and the SIFMA/Bloomberg OIS conventions everyone uses. See the per-field comments.

#include "rl/core/calendar.hpp"
#include "rl/core/daycount.hpp"
#include "rl/core/date.hpp"
#include "rl/core/index.hpp"
#include "rl/core/schedule.hpp"
#include "rl/core/tenor.hpp"

#include <string>
#include <string_view>

namespace rl::instruments {

/// Fixed-vs-overnight swap (OIS).
struct OisSwapConventions {
    core::Index index;                                   ///< the floating index (SOFR, FF, ESTR, ...)
    int spotLag = 2;                                     ///< business days from trade date to effective date
    std::string spotCalendar;                            ///< calendar for the spot lag
    std::string scheduleCalendar;                        ///< calendar for period-date adjustment
    core::BusinessDayConvention convention = core::BusinessDayConvention::ModifiedFollowing;
    core::Tenor fixedFrequency = core::Tenor::years(1);
    core::DayCount fixedDayCount = core::DayCount::Act360;
    core::Tenor floatFrequency = core::Tenor::years(1);
    core::DayCount floatDayCount = core::DayCount::Act360;
    core::StubRule stub = core::StubRule::ShortFront;
    bool endOfMonth = false;
    int paymentLag = 2;                                  ///< business days from accrual end to payment (both legs)
    /// Swaps with maturity ≤ this tenor are one period, paid at maturity (market: SOFR OIS ≤ 1Y).
    core::Tenor singlePeriodUpTo = core::Tenor::years(1);
};

/// Overnight-vs-overnight basis swap (FF vs SOFR). The spread is on the `spreadIndex` leg.
struct BasisSwapConventions {
    core::Index spreadIndex;                             ///< leg that carries the quoted spread (FF)
    core::Index flatIndex;                               ///< leg paid flat (SOFR)
    int spotLag = 2;
    std::string spotCalendar;
    std::string scheduleCalendar;
    core::BusinessDayConvention convention = core::BusinessDayConvention::ModifiedFollowing;
    core::Tenor frequency = core::Tenor::months(3);      ///< both legs; SIMPLIFICATION: desk convention, see .cpp
    core::StubRule stub = core::StubRule::ShortFront;
    bool endOfMonth = false;
    int paymentLag = 2;
    core::Tenor singlePeriodUpTo = core::Tenor::months(3);
};

/// Cash deposit on an overnight or term index.
struct DepositConventions {
    core::Index index;
    int spotLag = 2;
    std::string spotCalendar;
    core::BusinessDayConvention convention = core::BusinessDayConvention::ModifiedFollowing;
    bool endOfMonth = true;                              ///< money-market EOM rule
};

/// Three-month compounded-overnight future (CME SR3): reference period IMM date to next IMM date,
/// price = 100 − 100·R where R is the daily-compounded rate over the period, ACT/360.
struct FutureConventions {
    core::Index index;
    core::Tenor referenceTenor = core::Tenor::months(3);
    core::DayCount dayCount = core::DayCount::Act360;
};

struct Conventions {
    /// USD SOFR OIS: annual/annual ACT/360 both legs, MF, USD-SIFMA, T+2 (Fedwire), 2-day payment delay,
    /// single period up to 1Y. (Bloomberg USOSFR; ISDA 2021 "USD-SOFR-OIS Compound".)
    static OisSwapConventions usdSofrOis();
    /// USD Fed Funds OIS: same as SOFR with the FF index (annual, ACT/360, 2-day delay).
    static OisSwapConventions usdFedFundsOis();
    /// FF vs SOFR basis: FF (arithmetic average) + spread vs SOFR (compounded) flat, quarterly, ACT/360,
    /// 2-day payment delay. Spread on the FF leg.
    static BasisSwapConventions usdFedFundsSofrBasis();
    static DepositConventions usdSofrDeposit();
    static FutureConventions usdSofr3mFuture();

    /// Lookup by index name for each family; throws InputError for an index without conventions.
    static OisSwapConventions oisSwap(std::string_view indexName);
    static BasisSwapConventions basisSwap(std::string_view spreadIndex, std::string_view flatIndex);
    static DepositConventions deposit(std::string_view indexName);
    static FutureConventions future(std::string_view indexName);
};

// --- IMM dates (CME quarterly cycle): third Wednesday of March, June, September, December -------------
/// Third Wednesday of the given month (any month, not only the quarterly ones).
core::Date thirdWednesday(int year, unsigned month);
/// True if d is the third Wednesday of March/June/September/December.
bool isQuarterlyImm(core::Date d);
/// First quarterly IMM date strictly after d.
core::Date nextQuarterlyImm(core::Date d);
/// The n-th quarterly IMM date strictly after d (n ≥ 1); n = 1 is nextQuarterlyImm.
core::Date nthQuarterlyImm(core::Date d, int n);
/// Contract-month code letter for the quarterly cycle: H (Mar), M (Jun), U (Sep), Z (Dec). Throws otherwise.
char immMonthCode(unsigned month);
unsigned immMonthFromCode(char code);

/// d + tenor without any adjustment (days/weeks add days; months/years add months with clamping).
core::Date addTenor(core::Date d, core::Tenor t);
/// Spot date: trade date advanced by `lag` business days on `calendar`.
core::Date spotDate(core::Date tradeDate, int lag, const core::Calendar& calendar);

} // namespace rl::instruments
