#pragma once
// Shared fixture for the bootstrap tests: the shipped USD quote and fixing files, the instrument set,
// curve specs and a BootstrapProblem built from them.
#include "rl/bootstrap/bootstrap_problem.hpp"
#include "rl/bootstrap/instrument_set.hpp"
#include "rl/bootstrap/pillars.hpp"
#include "rl/core/calendar.hpp"
#include "rl/market/fixings.hpp"
#include "rl/market/quote.hpp"

namespace rl::bootstrap::test {

inline market::QuoteSet usdQuotes() { return market::QuoteSet::fromFile(core::dataDir() / "quotes" / "usd_2026-09-04.json"); }
inline market::FixingStore usdFixings() { return market::FixingStore::fromFile(core::dataDir() / "fixings" / "usd_2026-09-04.json"); }

inline BootstrapProblem usdProblem(market::CurveInterpolation interp = market::CurveInterpolation::LinearLogDf, bool includeTurns = false) {
    market::QuoteSet quotes = usdQuotes();
    InstrumentSetOptions io;
    io.includeTurns = includeTurns;
    InstrumentSet instruments = buildInstrumentSet(quotes, io);
    PillarOptions po;
    po.interpolation = interp;
    std::vector<CurveSpec> specs = makeCurveSpecs(instruments, quotes.asOf(), po);
    return BootstrapProblem(quotes.asOf(), std::move(specs), std::move(instruments), usdFixings(), std::move(quotes));
}

} // namespace rl::bootstrap::test
