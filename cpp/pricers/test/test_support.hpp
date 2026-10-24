#pragma once
#include "rl/core/tenor.hpp"
#include "rl/instruments/builder.hpp"
#include "rl/models/usd_market.hpp"
#include "rl/pricers/bermudan.hpp"

namespace rl::pricers::test {

inline const models::UsdMarket& usdMarket() {
    static const models::UsdMarket m = models::loadUsdMarket();
    return m;
}

inline instruments::Swap spotSwap(int years, double rate, bool payFixed = true, double notional = 1.0) {
    return instruments::OisSwapBuilder::usdSofr(usdMarket().asOf).id("swap").tenor(core::Tenor::years(years)).fixedRate(rate).notional(notional).payFixed(payFixed).build();
}

/// The reference trade: 10-year swap, callable annually from year 1, at `rate`, mean reversion `kappa`.
inline CalibratedBermudan bermudan10nc1(double rate = 0.035, bool payFixed = true, double kappa = 0.03) {
    return calibrateBermudan(usdMarket(), spotSwap(10, rate, payFixed), 1.0, 0, kappa);
}

/// The same problem restricted to a subset of its exercise dates.
inline BermudanProblem subset(const BermudanProblem& p, const std::vector<std::size_t>& keep) {
    std::vector<models::ExerciseRight> rights;
    for (std::size_t i : keep) rights.push_back(p.rights()[i]);
    return BermudanProblem(p.model(), std::move(rights));
}

} // namespace rl::pricers::test
