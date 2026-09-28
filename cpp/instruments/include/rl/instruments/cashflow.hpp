#pragma once
// Resolved cashflows: the atoms every trade is made of after trade resolution (dates adjusted, year
// fractions computed, notionals per coupon). Passive, market-free, model-free; a variant over a closed
// set traversed with std::visit.
//
// DESIGN: public-member aggregates, not classes with accessors. These structs have no invariant of
// their own (ordering is the Leg's invariant), they serialise field by field, and the visitor in
// linear/ reads them a million times per bootstrap. Accessors would add nothing but keystrokes.
//
// DESIGN: `tau` is stored, not recomputed from the dates. The cashflow is the *resolved* trade; linear/
// must be able to price a fixed coupon without a DayCount, and the day count is a property of the leg
// the coupon was generated from, not of the coupon.
//
// DESIGN: core::Index is carried by value on floating coupons (~100 bytes). The OIS kernel reads the
// fixing calendar, day count and compounding rule per coupon inside the hot loop; a by-name lookup
// there would be a map search per cashflow.
//
// The variant is closed at four alternatives in v1. v2 coupon types (CMS, range accrual, capped
// floater) are new alternatives, and every std::visit in linear/ then fails to compile until it
// handles them — the point of a closed variant over a virtual hierarchy.

#include "rl/core/date.hpp"
#include "rl/core/index.hpp"

#include <variant>

namespace rl::instruments {

struct FixedCoupon {
    core::Date accrualStart, accrualEnd, payment;
    double tau = 0.0;      ///< year fraction of [accrualStart, accrualEnd] under the leg's day count
    double notional = 0.0; ///< positive; the Leg's Side carries the sign
    double rate = 0.0;     ///< decimal
};

struct OvernightCoupon {
    core::Date accrualStart, accrualEnd, payment;
    double tau = 0.0;
    double notional = 0.0;
    core::Index index;         ///< overnight index; its compounding rule (compounded / averaged) applies
    double spread = 0.0;       ///< decimal, added to the compounded rate
    int lookbackDays = 0;      ///< observation shift in business days
    int lockoutDays = 0;       ///< last k fixings frozen at the k-th last
};

struct TermCoupon {
    core::Date fixingDate, accrualStart, accrualEnd, payment;
    double tau = 0.0;
    double notional = 0.0;
    core::Index index;         ///< term index (EURIBOR, Term SOFR); fixes once on fixingDate
    double spread = 0.0;
};

struct NotionalExchange {
    core::Date payment;
    double amount = 0.0;       ///< signed before the leg's sign: +receive, −pay
};

using Cashflow = std::variant<FixedCoupon, OvernightCoupon, TermCoupon, NotionalExchange>;

enum class Side { Pay, Receive };
enum class Settlement { Physical, Cash };

inline core::Date paymentDate(const Cashflow& cf) {
    return std::visit([](const auto& c) { return c.payment; }, cf);
}
/// Accrual end for coupons; the payment date for a notional exchange.
inline core::Date accrualEndOrPayment(const Cashflow& cf) {
    return std::visit(
        [](const auto& c) {
            if constexpr (std::is_same_v<std::decay_t<decltype(c)>, NotionalExchange>) return c.payment;
            else return c.accrualEnd;
        },
        cf);
}

} // namespace rl::instruments
