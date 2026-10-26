// Risk: the adjoint and bump routes agree, and the ladders say what the product's economics say.
#include "rl/pricers/risk.hpp"

#include "rl/linear/dv01.hpp"
#include "rl/linear/leg_pv.hpp"
#include "rl/pricers/bermudan_pde.hpp"
#include "test_support.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <algorithm>
#include <cmath>

using namespace rl::pricers;
using namespace rl::pricers::test;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {
double find(const std::vector<RiskRow>& rows, const std::string& id) {
    for (const auto& r : rows)
        if (r.id == id) return r.value;
    FAIL("no row " << id);
    return 0.0;
}
double total(const std::vector<RiskRow>& rows) {
    double s = 0.0;
    for (const auto& r : rows) s += r.value;
    return s;
}
} // namespace

TEST_CASE("Swap par delta: adjoint route equals bump-and-rebuild on every quote", "[pricers][risk][aad]") {
    const auto& m = usdMarket();
    const std::vector<rl::instruments::Swap> book{spotSwap(10, 0.030, true, 1e6), spotSwap(7, 0.041, false, 2e6), spotSwap(30, 0.036, true, 5e5)};
    const auto aad = swapParDeltaAdjoint(m, book), bump = swapParDeltaBump(m, book, 0.01), bump1bp = swapParDeltaBump(m, book);
    double largest = 0.0;
    for (const auto& r : aad) largest = std::max(largest, std::fabs(r.value));
    REQUIRE(aad.size() == bump.size());
    REQUIRE(aad.size() == m.problem.size());
    for (std::size_t i = 0; i < aad.size(); ++i) {
        INFO(aad[i].id);
        // The adjoint is the exact derivative; a small bump converges to it.
        CHECK_THAT(aad[i].value, WithinAbs(bump[i].value, 1e-4 * largest));
        // A one-sided 1 bp bump also carries the (cross-)gamma: up to ~1e-3 of the largest delta.
        CHECK_THAT(aad[i].value, WithinAbs(bump1bp[i].value, 2e-3 * largest));
    }
    // The 30y SOFR swap's risk lands on the 30y Fed Funds OIS and the 30y FF-SOFR basis quotes: there is
    // no 30y SOFR OIS in the curve, so those two instruments ARE the 30y SOFR point.
    CHECK(std::fabs(find(aad, "FF-OIS-30Y")) > 0.5 * largest);
    CHECK(std::fabs(find(aad, "FF-SOFR-30Y")) > 0.5 * largest);
}

TEST_CASE("A 10y payer swap's delta sits on the 10y quote and sums to its DV01", "[pricers][risk]") {
    const auto& m = usdMarket();
    const auto swap = spotSwap(10, 0.030, true, 1e6);
    const auto rows = swapParDeltaAdjoint(m, {swap});
    const double dv01 = rl::linear::parallelDv01(swap, m.curves, m.fixings);
    CHECK(dv01 > 0.0); // a payer gains when rates rise
    CHECK(find(rows, "SOFR-OIS-10Y") > 0.9 * dv01);
    CHECK(std::fabs(find(rows, "SOFR-OIS-20Y")) < 1e-6 * dv01);
    CHECK(std::fabs(find(rows, "FF-SOFR-10Y")) < 1e-6 * dv01); // a SOFR swap has no Fed Funds basis risk
    CHECK_THAT(total(rows), WithinRel(dv01, 0.05));
}

TEST_CASE("Bermudan delta: a payer is long rates, less than the swap it can enter, more than nothing", "[pricers][risk][bermudan]") {
    const auto& m = usdMarket();
    const auto swap = spotSwap(10, 0.035, true, 1e6);
    auto price = [&](const rl::models::UsdMarket& mk) { return priceBermudanPde(calibrateBermudan(mk, swap, 1.0, 0, 0.03).problem).value; };
    const auto rows = parDeltaBump(m, price);
    const double delta = total(rows);
    // The forward swap from the first call date: its delta is the cap on the option's delta.
    const auto first = calibrateBermudan(m, swap, 1.0, 0, 0.03).problem.rights().front();
    const double annuity = first.annuityValue() * 1e-4; // value of 1 bp on the swap entered at the first call
    CHECK(delta > 0.2 * annuity);
    CHECK(delta < annuity);
    // No exposure beyond the swap's maturity, and none to the Fed Funds basis.
    CHECK(std::fabs(find(rows, "SOFR-OIS-20Y")) < 1e-3 * delta);
    CHECK(std::fabs(find(rows, "FF-SOFR-5Y")) < 1e-3 * delta);
    // The receiver has the opposite sign.
    const auto rec = spotSwap(10, 0.035, false, 1e6);
    auto priceRec = [&](const rl::models::UsdMarket& mk) { return priceBermudanPde(calibrateBermudan(mk, rec, 1.0, 0, 0.03).problem).value; };
    CHECK(total(parDeltaBump(m, priceRec)) < 0.0);
}

TEST_CASE("Bermudan vega: positive in total, and bounded by the Europeans' own vegas", "[pricers][risk][vega]") {
    const auto& m = usdMarket();
    const auto swap = spotSwap(10, 0.035, true, 1e6);
    auto pde = [](const BermudanProblem& p) { return priceBermudanPde(p).value; };
    const auto rows = bermudanVegaByExpiry(m, swap, 1.0, 0, 0.03, pde);
    REQUIRE(rows.size() == 9);
    CHECK(total(rows) > 0.0);
    // Vega to one co-terminal European cannot exceed that European's own vega (you hold at most one of it),
    // though individual buckets may be negative: more vol early makes early exercise likelier and
    // takes value from the later dates.
    const CalibratedBermudan b = calibrateBermudan(m, swap, 1.0, 0, 0.03);
    for (std::size_t i = 0; i < rows.size(); ++i) {
        const auto& r = b.problem.rights()[i];
        const double sd = 0.0090 * std::sqrt(r.time);
        const double europeanVega = r.annuityValue() * std::sqrt(r.time) * std::exp(-0.5 * std::pow((r.forwardRate() - r.strike) / sd, 2)) / std::sqrt(2 * M_PI) * 1e-4;
        INFO(rows[i].id << " vega " << rows[i].value << " european " << europeanVega);
        CHECK(rows[i].value < 1.1 * europeanVega);
    }
    // A parallel 1 bp vol shift moves the price by about the sum of the buckets.
    auto full = [&](const rl::models::UsdMarket& mk) { return priceBermudanPde(calibrateBermudan(mk, swap, 1.0, 0, 0.03).problem).value; };
    const double parallel = full(m.withVolShift(1e-4)) - full(m);
    CHECK_THAT(parallel, WithinRel(total(rows), 0.05));
}
