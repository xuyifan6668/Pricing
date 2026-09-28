#pragma once
// Trade builders: conventions + a few numbers → a fully dated trade.
//   OisSwapBuilder::usdSofr(tradeDate).id("SOFR-10Y").tenor(Tenor::years(10)).fixedRate(0.0345).build()
// produces the market-standard swap. Everything tedious lives here: spot lag, schedule generation with
// the single-period rule for short tenors, payment delays, year fractions, fixing dates, IMM periods.
//
// Builders are single-use value objects; the trade types they produce are in cashflow.hpp, leg.hpp,
// swap.hpp, basis_swap.hpp, deposit.hpp, fra.hpp, future.hpp.

#include "rl/core/calendar.hpp"
#include "rl/core/date.hpp"
#include "rl/core/schedule.hpp"
#include "rl/core/tenor.hpp"
#include "rl/instruments/basis_swap.hpp"
#include "rl/instruments/cashflow.hpp"
#include "rl/instruments/conventions.hpp"
#include "rl/instruments/deposit.hpp"
#include "rl/instruments/fra.hpp"
#include "rl/instruments/future.hpp"
#include "rl/instruments/leg.hpp"
#include "rl/instruments/swap.hpp"

#include <optional>
#include <string>
#include <vector>

namespace rl::instruments {

// --- leg factories (also used directly by tests) --------------------------------------------------

/// One FixedCoupon per schedule period; tau under `dayCount` on the adjusted accrual dates.
Leg makeFixedLeg(const core::Schedule& schedule, core::DayCount dayCount, double rate, double notional, Side side,
                 core::Currency currency);

/// One OvernightCoupon per period; tau under index.dayCount; the compounding/averaging rule travels
/// with the index (core::Index::compounding).
Leg makeOvernightLeg(const core::Schedule& schedule, const core::Index& index, double spread, double notional, Side side,
                     core::Currency currency, int lookbackDays = 0, int lockoutDays = 0);

/// One TermCoupon per period; fixing date = accrual start moved back index.fixingLag business days on
/// the index's fixing calendar.
Leg makeTermLeg(const core::Schedule& schedule, const core::Index& index, double spread, double notional, Side side,
                core::Currency currency);

/// Approximate calendar days in a tenor, for ordering tenors of different units (1Y ≈ 365.25 d).
double tenorDaysApprox(core::Tenor t) noexcept;

// --- OIS swap --------------------------------------------------------------------------------------

class OisSwapBuilder {
public:
    OisSwapBuilder(OisSwapConventions conventions, core::Date tradeDate);
    static OisSwapBuilder usdSofr(core::Date tradeDate) { return OisSwapBuilder(Conventions::usdSofrOis(), tradeDate); }
    static OisSwapBuilder usdFedFunds(core::Date tradeDate) { return OisSwapBuilder(Conventions::usdFedFundsOis(), tradeDate); }

    OisSwapBuilder& id(std::string id) { id_ = std::move(id); return *this; }
    OisSwapBuilder& tenor(core::Tenor t) { tenor_ = t; maturity_.reset(); return *this; }
    OisSwapBuilder& effective(core::Date d) { effective_ = d; return *this; }   ///< overrides the spot rule
    OisSwapBuilder& maturity(core::Date d) { maturity_ = d; tenor_.reset(); return *this; } ///< unadjusted termination
    OisSwapBuilder& fixedRate(double r) { fixedRate_ = r; return *this; }
    OisSwapBuilder& notional(double n) { notional_ = n; return *this; }
    OisSwapBuilder& payFixed(bool p) { payFixed_ = p; return *this; }           ///< true: payer swap (default)
    OisSwapBuilder& floatSpread(double s) { spread_ = s; return *this; }
    OisSwapBuilder& lookback(int days) { lookback_ = days; return *this; }
    OisSwapBuilder& lockout(int days) { lockout_ = days; return *this; }

    const OisSwapConventions& conventions() const noexcept { return conv_; }
    core::Date tradeDate() const noexcept { return tradeDate_; }
    core::Date effectiveDate() const;
    core::Date terminationDate() const;        ///< unadjusted
    core::Schedule fixedSchedule() const;
    core::Schedule floatSchedule() const;
    /// Throws InputError if id, tenor/maturity or fixed rate are missing.
    Swap build() const;

private:
    core::Schedule schedule(core::Tenor frequency) const;
    OisSwapConventions conv_;
    core::Date tradeDate_;
    std::string id_;
    std::optional<core::Tenor> tenor_;
    std::optional<core::Date> effective_;
    std::optional<core::Date> maturity_;
    std::optional<double> fixedRate_;
    double notional_ = 1.0e6;
    bool payFixed_ = true;
    double spread_ = 0.0;
    int lookback_ = 0;
    int lockout_ = 0;
};

// --- overnight/overnight basis swap ------------------------------------------------------------------

class BasisSwapBuilder {
public:
    BasisSwapBuilder(BasisSwapConventions conventions, core::Date tradeDate);
    static BasisSwapBuilder usdFedFundsSofr(core::Date tradeDate) { return BasisSwapBuilder(Conventions::usdFedFundsSofrBasis(), tradeDate); }

    BasisSwapBuilder& id(std::string id) { id_ = std::move(id); return *this; }
    BasisSwapBuilder& tenor(core::Tenor t) { tenor_ = t; maturity_.reset(); return *this; }
    BasisSwapBuilder& effective(core::Date d) { effective_ = d; return *this; }
    BasisSwapBuilder& maturity(core::Date d) { maturity_ = d; tenor_.reset(); return *this; }
    BasisSwapBuilder& spread(double s) { spread_ = s; return *this; }            ///< on the spread leg, decimal
    BasisSwapBuilder& notional(double n) { notional_ = n; return *this; }
    BasisSwapBuilder& paySpreadLeg(bool p) { paySpreadLeg_ = p; return *this; } ///< default true: pay FF+s, receive SOFR

    const BasisSwapConventions& conventions() const noexcept { return conv_; }
    core::Date effectiveDate() const;
    core::Date terminationDate() const;
    core::Schedule schedule() const;
    /// legs()[0] is the spread leg, legs()[1] the flat leg. Throws if id or tenor/maturity are missing.
    BasisSwap build() const;

private:
    BasisSwapConventions conv_;
    core::Date tradeDate_;
    std::string id_;
    std::optional<core::Tenor> tenor_;
    std::optional<core::Date> effective_;
    std::optional<core::Date> maturity_;
    double spread_ = 0.0;
    double notional_ = 1.0e6;
    bool paySpreadLeg_ = true;
};

// --- deposit ---------------------------------------------------------------------------------------

class DepositBuilder {
public:
    DepositBuilder(DepositConventions conventions, core::Date tradeDate);
    static DepositBuilder usdSofr(core::Date tradeDate) { return DepositBuilder(Conventions::usdSofrDeposit(), tradeDate); }

    DepositBuilder& id(std::string id) { id_ = std::move(id); return *this; }
    DepositBuilder& tenor(core::Tenor t) { tenor_ = t; return *this; }
    DepositBuilder& start(core::Date d) { start_ = d; return *this; }
    DepositBuilder& end(core::Date d) { end_ = d; return *this; }   ///< explicit (adjusted) end; overrides tenor
    DepositBuilder& rate(double r) { rate_ = r; return *this; }
    DepositBuilder& notional(double n) { notional_ = n; return *this; }
    Deposit build() const;

private:
    DepositConventions conv_;
    core::Date tradeDate_;
    std::string id_;
    std::optional<core::Tenor> tenor_;
    std::optional<core::Date> start_, end_;
    std::optional<double> rate_;
    double notional_ = 1.0e6;
};

// --- FRA -------------------------------------------------------------------------------------------

class FraBuilder {
public:
    explicit FraBuilder(core::Index index);
    FraBuilder& id(std::string id) { id_ = std::move(id); return *this; }
    FraBuilder& start(core::Date d) { start_ = d; return *this; }   ///< adjusted accrual start
    FraBuilder& end(core::Date d) { end_ = d; return *this; }       ///< adjusted accrual end
    FraBuilder& rate(double r) { rate_ = r; return *this; }
    FraBuilder& notional(double n) { notional_ = n; return *this; }
    Fra build() const;

private:
    core::Index index_;
    std::string id_;
    std::optional<core::Date> start_, end_;
    std::optional<double> rate_;
    double notional_ = 1.0e6;
};

// --- 3M compounded-overnight future -------------------------------------------------------------------

class FutureBuilder {
public:
    explicit FutureBuilder(FutureConventions conventions);
    static FutureBuilder usdSofr3m() { return FutureBuilder(Conventions::usdSofr3mFuture()); }

    FutureBuilder& id(std::string id) { id_ = std::move(id); return *this; }
    /// Quarterly contract month: reference period = IMM(year, month) → next quarterly IMM.
    FutureBuilder& contract(int year, unsigned month);
    FutureBuilder& start(core::Date d) { start_ = d; return *this; } ///< explicit reference start
    FutureBuilder& end(core::Date d) { end_ = d; return *this; }     ///< explicit reference end
    FutureBuilder& price(double p) { price_ = p; return *this; }
    FutureBuilder& convexityAdjustmentBp(double bp) { convexityBp_ = bp; return *this; }
    Future build() const;

private:
    FutureConventions conv_;
    std::string id_;
    std::optional<core::Date> start_, end_;
    std::optional<double> price_;
    double convexityBp_ = 0.0;
};

} // namespace rl::instruments
