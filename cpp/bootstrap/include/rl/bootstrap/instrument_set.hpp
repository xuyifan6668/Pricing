#pragma once
// Quotes → dated bootstrap instruments. Each quote becomes one trade built with the market
// conventions for its index, tagged with the curve it pins and the pillar date it pins it at.
//
// Ordering: instruments are sorted by (curve, in order of first appearance in the quote set) and then by
// maturity, so that residual i and parameter i refer to the same pillar once makeCurveSpecs() has run.
// Two instruments on one curve with the same maturity make the Jacobian singular by construction and are
// rejected here, not discovered in the solver.
//
// Curve naming: a curve is named after the index it projects ("USD-SOFR", "USD-FEDFUNDS"); a basis swap
// pins the curve of its spread index; a turn pins the curve of its index.

#include "rl/core/date.hpp"
#include "rl/instruments/basis_swap.hpp"
#include "rl/instruments/deposit.hpp"
#include "rl/instruments/fra.hpp"
#include "rl/instruments/future.hpp"
#include "rl/instruments/swap.hpp"
#include "rl/market/quote.hpp"

#include <string>
#include <variant>
#include <vector>

namespace rl::bootstrap {

/// Turn-of-year: the overnight forward over [start, end] is `jump` above the smooth curve.
struct TurnOfYear {
    std::string id;
    std::string index;
    core::Date start;
    core::Date end;
    double jump; ///< rate units (8 bp → 0.0008)
};

using BootstrapTrade = std::variant<instruments::Deposit, instruments::Fra, instruments::Future, instruments::Swap, instruments::BasisSwap, TurnOfYear>;

struct BootstrapInstrument {
    std::string id;              ///< the quote id
    market::QuoteType type;
    std::string curveName;       ///< the curve this instrument pins
    core::Date maturity;         ///< pillar date on that curve (last accrual end; futures: reference end)
    double quote;                ///< decimal rate / spread / jump, or price for futures
    market::Quote source;        ///< the quote it came from (unit, convexity, ...)
    BootstrapTrade trade;
};

using InstrumentSet = std::vector<BootstrapInstrument>;

struct InstrumentSetOptions {
    double notional = 1.0;        ///< per-unit notional keeps PV residuals O(quote error)
    bool includeTurns = false;    ///< include turn-of-year quotes as jump parameters (off by default: notebook 04 shows the curve without them first)
};

/// Builds every quote in the set (trade date = quotes.asOf()). Throws InputError (with the quote id as
/// context) for unsupported index/type combinations or duplicate (curve, maturity) pillars.
InstrumentSet buildInstrumentSet(const market::QuoteSet& quotes, const InstrumentSetOptions& options = {});

/// Curve names in order of first appearance.
std::vector<std::string> curveNames(const InstrumentSet& instruments);
/// Instruments pinning one curve, in maturity order.
std::vector<const BootstrapInstrument*> instrumentsOn(const InstrumentSet& instruments, std::string_view curveName);
std::string_view tradeKind(const BootstrapTrade& t) noexcept; // "deposit", "fra", "future", "swap", "basis_swap", "turn"

} // namespace rl::bootstrap
