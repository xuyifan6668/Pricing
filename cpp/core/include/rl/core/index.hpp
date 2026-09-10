#pragma once
// Interest-rate index descriptions. Pure data: an Index knows its conventions, never its fixings
// or its curve (those live in market::). Compare by name.
//
// Two families:
//   Overnight (SOFR, FF, ESTR, SONIA, TONA): one fixing per business day, published the next
//     business day. A coupon on such an index compounds or averages the daily fixings over its
//     accrual period, with lookback / lockout / payment-delay conventions.
//   Term (EURIBOR-3M/6M, Term SOFR-3M): one fixing per coupon, set `fixingLag` business days before
//     the accrual start, accruing simply.

#include "rl/core/currency.hpp"
#include "rl/core/daycount.hpp"
#include "rl/core/tenor.hpp"

#include <string>
#include <string_view>

namespace rl::core {

enum class IndexKind { Overnight, Term };

/// How daily fixings aggregate into a coupon rate (overnight indices only).
enum class Compounding {
    Simple,          ///< Term indices: rate × τ.
    CompoundedDaily, ///< ∏(1 + rᵢ τᵢ) − 1 over the period. SOFR/ESTR/SONIA/TONA OIS.
    AveragedDaily,   ///< Σ rᵢ τᵢ / τ. Fed Funds legs of FF–SOFR basis swaps (arithmetic average).
};

std::string_view toString(IndexKind k);
std::string_view toString(Compounding c);

struct Index {
    std::string name;            ///< "USD-SOFR", "USD-FEDFUNDS", "EUR-ESTR", "GBP-SONIA", "JPY-TONA", "EUR-EURIBOR-3M", ...
    Currency currency = Currency::USD;
    IndexKind kind = IndexKind::Overnight;
    Tenor tenor = Tenor::days(1);        ///< 1D for overnight; 3M/6M for term indices.
    DayCount dayCount = DayCount::Act360;
    std::string fixingCalendar;          ///< Registry name; the days on which the index publishes.
    int fixingLag = 0;                   ///< Term: business days from fixing to accrual start (EURIBOR 2).
    Compounding compounding = Compounding::CompoundedDaily;
    // Overnight coupon conventions (market standard for the swap, not the index, but they travel with it):
    int lookbackDays = 0;      ///< Observation shift/lookback in business days (SONIA swaps: 0 since 2021; some loans 5).
    int lockoutDays = 0;       ///< Last k fixings frozen at the k-th last rate (SOFR futures/FRNs use 2; swaps 0).
    int paymentDelayDays = 0;  ///< Business days from accrual end to payment (SOFR/FF/ESTR OIS: 2; SONIA/TONA: 0/1... see below).

    bool isOvernight() const noexcept { return kind == IndexKind::Overnight; }
    friend bool operator==(const Index& a, const Index& b) noexcept { return a.name == b.name; }
    friend bool operator!=(const Index& a, const Index& b) noexcept { return !(a == b); }

    // --- market-standard definitions --------------------------------------------------------
    static Index sofr();       ///< ACT/360, USD-SIFMA fixing days, compounded, 2-day payment delay.
    static Index fedFunds();   ///< ACT/360, USD-FEDWIRE, arithmetic average in basis swaps, 2-day delay.
    static Index estr();       ///< ACT/360, TARGET2, compounded, 1-day payment delay.
    static Index sonia();      ///< ACT/365F, GBP, compounded, 0-day delay.
    static Index tona();       ///< ACT/365F, JPY, compounded, 2-day delay.
    static Index euribor3m();  ///< ACT/360, TARGET2, 2-day fixing lag, simple.
    static Index euribor6m();
    static Index termSofr3m(); ///< ACT/360, USD-SIFMA, 2-day fixing lag, simple.
    /// Lookup by name (the strings in the `name` field above). Throws InputError.
    static Index get(std::string_view name);
};

} // namespace rl::core
