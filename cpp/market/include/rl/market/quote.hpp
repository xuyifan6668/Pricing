#pragma once
// Market quotes as they arrive from the screen: an id, an instrument type, the identifying dates or
// tenor, a number and its unit. No pricing, no conventions beyond "what kind of thing is this". The
// bootstrap turns a QuoteSet into instruments (bootstrap/instrument_set.hpp) using the conventions table
// in instruments/conventions.hpp; the snapshot stores the QuoteSet it was built from so that ∂F/∂q has
// a q to differentiate with respect to.
//
// Value type; validated on construction; immutable afterwards except through with*() which return copies.
//
// JSON (data/quotes/usd_2026-09-04.json):
//   {"asOf": "2026-09-04", "currency": "USD", "quotes": [
//     {"id": "SOFR-OIS-1Y", "type": "ois_swap",  "index": "USD-SOFR", "tenor": "1Y", "value": 3.62, "unit": "pct"},
//     {"id": "SR3-Z6",      "type": "future",    "index": "USD-SOFR", "start": "2026-12-16", "end": "2027-03-17",
//                            "value": 96.35, "unit": "price", "convexity_bp": 0.4},
//     {"id": "FF-SOFR-1Y",  "type": "basis_swap","index": "USD-FEDFUNDS", "reference_index": "USD-SOFR",
//                            "tenor": "1Y", "value": -2.5, "unit": "bp"},
//     {"id": "FF-SOFR-30Y", "type": "basis_swap","index": "USD-FEDFUNDS", "reference_index": "USD-SOFR",
//                            "tenor": "30Y", "value": -4.0, "unit": "bp", "curve": "USD-SOFR"},
//     {"id": "TURN-2026",   "type": "turn", "index": "USD-SOFR", "start": "2026-12-31", "end": "2027-01-04",
//                            "value": 8.0, "unit": "bp"} ]}

#include "rl/core/currency.hpp"
#include "rl/core/date.hpp"
#include "rl/core/tenor.hpp"

#include <cstddef>
#include <filesystem>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace rl::market {

class Hasher;

enum class QuoteType {
    Deposit,   ///< cash deposit from spot: rate, tenor
    Fra,       ///< forward rate agreement on a term index: rate, start/end
    Future,    ///< 3M compounded-overnight future: price = 100 − 100·R, start/end are the IMM reference dates
    OisSwap,   ///< fixed vs compounded overnight: par fixed rate, tenor
    BasisSwap, ///< two overnight legs: spread on `index` leg vs `referenceIndex` flat, tenor
    Turn,      ///< turn-of-year jump in the overnight forward over [start, end], in rate units
};
enum class QuoteUnit { Decimal, Percent, BasisPoints, Price };

std::string_view toString(QuoteType t);
std::string_view toString(QuoteUnit u);
QuoteType quoteTypeFromString(std::string_view s); // "deposit","fra","future","ois_swap","basis_swap","turn"
QuoteUnit quoteUnitFromString(std::string_view s); // "decimal","pct","percent","bp","bps","price"

struct Quote {
    std::string id;
    QuoteType type = QuoteType::OisSwap;
    std::string index;                    ///< index the instrument is on ("USD-SOFR"); basis: the spread leg's index
    std::string referenceIndex;           ///< basis swaps only: the flat leg's index
    std::string curve;                    ///< which curve this instrument pins; empty = its index (basis: the spread index).
                                          ///< Set it when a basis swap pins the *other* curve, e.g. 25Y/30Y FF–SOFR basis
                                          ///< pinning USD-SOFR where no SOFR OIS is quoted (that is what couples the curves).
    std::optional<core::Tenor> tenor;     ///< deposits, swaps, basis swaps: maturity from spot
    std::optional<core::Date> start;      ///< FRA, future, turn
    std::optional<core::Date> end;        ///< FRA, future, turn
    double value = 0.0;                   ///< as quoted, in `unit`
    QuoteUnit unit = QuoteUnit::Percent;
    double convexityAdjustmentBp = 0.0;   ///< futures: forward = futures rate − adjustment

    /// `value` as a decimal rate (Percent/100, BasisPoints/1e4, Decimal as is). Throws for Price.
    double asRate() const;
    /// `value` for Price quotes; throws otherwise.
    double asPrice() const;
    /// The quoted number in its own unit, from a decimal rate (inverse of asRate; identity for Price).
    double fromRate(double rate) const;
    /// Multiplicative factor such that (value in unit) = factor × (value in decimal rate); 1 for Price.
    double unitScale() const noexcept;
    /// Throws InputError with the id in the message if required fields for `type` are missing.
    void validate() const;
};

class QuoteSet {
public:
    QuoteSet() = default;
    /// Validates every quote and rejects duplicate ids.
    QuoteSet(core::Date asOf, core::Currency currency, std::vector<Quote> quotes);

    static QuoteSet fromJson(std::string_view text);
    static QuoteSet fromFile(const std::filesystem::path& file);
    std::string toJson() const;

    core::Date asOf() const noexcept { return asOf_; }
    core::Currency currency() const noexcept { return currency_; }
    const std::vector<Quote>& quotes() const noexcept { return quotes_; }
    std::size_t size() const noexcept { return quotes_.size(); }
    bool empty() const noexcept { return quotes_.empty(); }
    const Quote& operator[](std::size_t i) const noexcept { return quotes_[i]; }
    std::vector<Quote>::const_iterator begin() const noexcept { return quotes_.begin(); }
    std::vector<Quote>::const_iterator end() const noexcept { return quotes_.end(); }

    bool contains(std::string_view id) const noexcept;
    std::size_t indexOf(std::string_view id) const; // throws InputError if absent
    const Quote& find(std::string_view id) const;   // throws InputError if absent
    std::vector<Quote> ofType(QuoteType t) const;
    std::vector<std::string> ids() const;

    /// Copies with new `value`s (same order and count). What the bump-and-rebuild driver uses.
    QuoteSet withValues(const std::vector<double>& values) const;
    QuoteSet withValue(std::string_view id, double value) const;
    /// Copy with one quote's value bumped by `bump` in the quote's own unit.
    QuoteSet bumped(std::string_view id, double bump) const;

    /// Deterministic contribution to a snapshot hash: asOf, currency, every id and value.
    void mixInto(Hasher& h) const noexcept;

private:
    core::Date asOf_;
    core::Currency currency_ = core::Currency::USD;
    std::vector<Quote> quotes_;
};

} // namespace rl::market
