#pragma once
// Shared fixtures: the USD market, a calibrated one-factor model, standard trades and credit curves.
#include "rl/core/tenor.hpp"
#include "rl/instruments/builder.hpp"
#include "rl/models/usd_market.hpp"
#include "rl/pricers/bermudan.hpp"
#include "rl/xva/credit.hpp"

#include <string>
#include <vector>

namespace rl::xva::test {

inline const models::UsdMarket& usdMarket() {
    static const models::UsdMarket m = models::loadUsdMarket();
    return m;
}

inline instruments::Swap spotSwap(int years, double rate, bool payFixed = true, double notional = 1.0, const std::string& id = "swap") {
    return instruments::OisSwapBuilder::usdSofr(usdMarket().asOf).id(id).tenor(core::Tenor::years(years)).fixedRate(rate).notional(notional).payFixed(payFixed).build();
}

/// LGM at 3% mean reversion, calibrated to the co-terminal strip of the 10-year swap at 3.50%.
inline const models::Lgm1F& model() {
    static const models::Lgm1F m = pricers::calibrateBermudan(usdMarket(), spotSwap(10, 0.035), 1.0, 0, 0.03).problem.model();
    return m;
}

inline double discount(double t) { return usdMarket().curves.discount().df(t); }

/// A BBB-like CDS curve: 80 bp at one year rising to 150 bp at ten.
inline std::vector<CdsQuote> cdsQuotes() { return {{1.0, 0.0080}, {3.0, 0.0105}, {5.0, 0.0125}, {7.0, 0.0140}, {10.0, 0.0150}}; }

inline HazardCurve counterpartyCurve() { return HazardCurve::bootstrap(cdsQuotes(), 0.4, discount); }

} // namespace rl::xva::test
