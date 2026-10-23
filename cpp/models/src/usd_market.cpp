#include "rl/models/usd_market.hpp"

#include "rl/bootstrap/instrument_set.hpp"
#include "rl/bootstrap/pillars.hpp"
#include "rl/core/calendar.hpp"
#include "rl/core/error.hpp"
#include "rl/market/quote.hpp"

namespace rl::models {

namespace {
UsdMarket assemble(bootstrap::BootstrapProblem problem, market::SwaptionVolQuotes volQuotes) {
    const bootstrap::BootstrapResult solved = bootstrap::solve(problem);
    if (!solved.converged) throw NumericalError("loadUsdMarket: the curve bootstrap did not converge");
    vol::CubeFit fit = vol::fitCube(volQuotes, vol::SabrFitOptions{});
    market::CurveSet<double> curves = problem.curves(solved.z);
    market::FixingStore fixings = problem.fixings();
    market::SwaptionCube<double> cube = fit.cube;
    const core::Date asOf = problem.asOf();
    return UsdMarket{asOf, std::move(problem), solved.z, std::move(curves), std::move(fixings), std::move(volQuotes), std::move(cube), std::move(fit)};
}
} // namespace

UsdMarket loadUsdMarket(market::CurveInterpolation interpolation) {
    market::QuoteSet quotes = market::QuoteSet::fromFile(core::dataDir() / "quotes" / "usd_2026-09-04.json");
    market::FixingStore fixings = market::FixingStore::fromFile(core::dataDir() / "fixings" / "usd_2026-09-04.json");
    bootstrap::InstrumentSet instruments = bootstrap::buildInstrumentSet(quotes);
    bootstrap::PillarOptions po;
    po.interpolation = interpolation;
    std::vector<bootstrap::CurveSpec> specs = bootstrap::makeCurveSpecs(instruments, quotes.asOf(), po);
    bootstrap::BootstrapProblem problem(quotes.asOf(), std::move(specs), std::move(instruments), std::move(fixings), std::move(quotes));
    return assemble(std::move(problem), market::SwaptionVolQuotes::fromFile(core::dataDir() / "vol" / "usd_swaption_cube_2026-09-04.json"));
}

UsdMarket UsdMarket::withCurveQuotes(const math::Vector& quotes) const { return assemble(problem.withQuotes(quotes), volQuotes); }

UsdMarket UsdMarket::withVolShift(double bump) const {
    std::vector<market::SwaptionVolNodeQuote> nodes = volQuotes.nodes();
    for (auto& n : nodes) {
        n.atmNormalVol += bump;
        for (auto& p : n.smile) p.normalVol += bump;
    }
    return assemble(problem, market::SwaptionVolQuotes(volQuotes.asOf(), volQuotes.currency(), volQuotes.index(), std::move(nodes)));
}

} // namespace rl::models
