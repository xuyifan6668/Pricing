#pragma once
// Shared fixtures: the USD market (built once per test binary) and a few standard trades.
#include "rl/core/calendar.hpp"
#include "rl/core/tenor.hpp"
#include "rl/instruments/builder.hpp"
#include "rl/models/flows.hpp"
#include "rl/models/usd_market.hpp"

namespace rl::models::test {

inline const UsdMarket& usdMarket() {
    static const UsdMarket m = loadUsdMarket();
    return m;
}

/// A spot-starting USD SOFR swap, notional 1.
inline instruments::Swap spotSwap(int years, double rate, bool payFixed = true) {
    return instruments::OisSwapBuilder::usdSofr(usdMarket().asOf).id("swap").tenor(core::Tenor::years(years)).fixedRate(rate).notional(1.0).payFixed(payFixed).build();
}

/// The annual exercise rights of a 10-year swap callable after one year ("10nc1"): 9 rights.
inline std::vector<ExerciseRight> rights10nc1(double rate, bool payFixed = true, int noticeDays = 0) {
    const instruments::Swap swap = spotSwap(10, rate, payFixed);
    std::vector<ExerciseRight> out;
    for (core::Date d : coterminalExerciseDates(swap, 1.0, noticeDays, "USD-SIFMA")) out.push_back(makeExerciseRight(swap, usdMarket().curves, d));
    return out;
}

} // namespace rl::models::test
