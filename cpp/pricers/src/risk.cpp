#include "rl/pricers/risk.hpp"

#include "rl/aad/real.hpp"
#include "rl/aad/tape.hpp"
#include "rl/bootstrap/jacobian.hpp"
#include "rl/core/error.hpp"
#include "rl/linear/leg_pv.hpp"
#include "rl/models/lgm_swaption.hpp"

#include <cstdint>

namespace rl::pricers {

namespace {
/// +1 bp of a quote in the units the bootstrap holds it: a rate goes up 1e-4; a futures price goes
/// DOWN 0.01, so that every row means "this point of the curve 1 bp higher".
double onebp(const bootstrap::BootstrapInstrument& bi) { return bi.type == market::QuoteType::Future ? -0.01 : 1e-4; }
} // namespace

std::vector<RiskRow> swapParDeltaAdjoint(const models::UsdMarket& market, const std::vector<instruments::Swap>& swaps) {
    const bootstrap::BootstrapProblem& problem = market.problem;
    const std::size_t n = problem.parameterCount();
    aad::Tape& tape = aad::Tape::current();
    const aad::Tape::Position start = tape.mark();
    math::Vector dPdz(static_cast<Eigen::Index>(n));
    {
        // Record: nodes → curves → PV of every swap. One tape, one output.
        std::vector<aad::Real> z;
        z.reserve(n);
        for (std::size_t j = 0; j < n; ++j) z.push_back(aad::Real::input(market.z[static_cast<Eigen::Index>(j)]));
        const market::CurveSet<aad::Real> curves = problem.curves(z);
        aad::Real pv = 0.0;
        for (const auto& s : swaps) pv += linear::pv(s, curves, market.fixings);
        // Sweep: ∂PV/∂(every node) in one backward pass.
        if (pv.isActive()) {
            tape.clearAdjointsFrom(static_cast<std::uint32_t>(start.nodes));
            tape.seedAdjoint(pv.node(), 1.0);
            tape.propagate(pv.node(), static_cast<std::uint32_t>(start.nodes));
        }
        for (std::size_t j = 0; j < n; ++j) dPdz[static_cast<Eigen::Index>(j)] = pv.isActive() ? tape.adjoint(z[j].node()) : 0.0;
    }
    tape.rewindTo(start);
    // Node risk → par-quote risk by the implicit-function theorem (one transposed solve).
    const math::Vector dPdq = bootstrap::parRisk(problem, market.z, dPdz);
    std::vector<RiskRow> out;
    for (std::size_t i = 0; i < problem.size(); ++i)
        out.push_back({problem.instrumentIds()[i], dPdq[static_cast<Eigen::Index>(i)] * onebp(problem.instruments()[i])});
    return out;
}

std::vector<RiskRow> parDeltaBump(const models::UsdMarket& market, const std::function<double(const models::UsdMarket&)>& price, double bumpBp) {
    const double base = price(market);
    const math::Vector q = market.problem.quoteValues();
    std::vector<RiskRow> out;
    for (std::size_t i = 0; i < market.problem.size(); ++i) {
        math::Vector qb = q;
        qb[static_cast<Eigen::Index>(i)] += bumpBp * onebp(market.problem.instruments()[i]);
        out.push_back({market.problem.instrumentIds()[i], (price(market.withCurveQuotes(qb)) - base) / bumpBp});
    }
    return out;
}

std::vector<RiskRow> swapParDeltaBump(const models::UsdMarket& market, const std::vector<instruments::Swap>& swaps, double bumpBp) {
    // Curves only: re-bootstrap per quote, without re-fitting the vol cube a swap does not use.
    auto pv = [&](const market::CurveSet<double>& curves) {
        double v = 0.0;
        for (const auto& s : swaps) v += linear::pv(s, curves, market.fixings);
        return v;
    };
    const double base = pv(market.curves);
    const math::Vector q = market.problem.quoteValues();
    std::vector<RiskRow> out;
    for (std::size_t i = 0; i < market.problem.size(); ++i) {
        math::Vector qb = q;
        qb[static_cast<Eigen::Index>(i)] += bumpBp * onebp(market.problem.instruments()[i]);
        const bootstrap::BootstrapProblem bumped = market.problem.withQuotes(qb);
        const bootstrap::BootstrapResult r = bootstrap::solve(bumped);
        if (!r.converged) throw NumericalError("swapParDeltaBump: bumped bootstrap did not converge for " + market.problem.instrumentIds()[i]);
        out.push_back({market.problem.instrumentIds()[i], (pv(bumped.curves(r.z)) - base) / bumpBp});
    }
    return out;
}

std::vector<RiskRow> bermudanVegaByExpiry(const models::UsdMarket& market, const instruments::Swap& swap, double firstCallYears, int noticeDays, double kappa,
                                          const std::function<double(const BermudanProblem&)>& price) {
    const CalibratedBermudan base = calibrateBermudan(market, swap, firstCallYears, noticeDays, kappa);
    const double v0 = price(base.problem);
    const std::vector<models::CalibrationInstrument> strip = models::makeStrip(base.problem.rights(), market.cube);
    std::vector<RiskRow> out;
    for (std::size_t i = 0; i < strip.size(); ++i) {
        std::vector<models::CalibrationInstrument> bumped = strip;
        bumped[i].marketVol += 1e-4;
        bumped[i].marketPrice = models::marketSwaptionPrice(bumped[i].right, bumped[i].marketVol);
        const models::LgmCalibration cal = models::calibrateLgm(kappa, bumped);
        out.push_back({strip[i].right.label, price(BermudanProblem(cal.model, base.problem.rights())) - v0});
    }
    return out;
}

} // namespace rl::pricers
