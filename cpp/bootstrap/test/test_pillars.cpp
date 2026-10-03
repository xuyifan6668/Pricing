#include "rl/bootstrap/pillars.hpp"

#include "rl/core/daycount.hpp"
#include "rl/core/error.hpp"
#include "test_support.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace rl::bootstrap;
using namespace rl::bootstrap::test;
using rl::core::Date;
using Catch::Matchers::WithinAbs;

TEST_CASE("Curve specs: one pillar per instrument, ACT/365F times, discount role, interpolation overrides", "[bootstrap][pillars]") {
    auto quotes = usdQuotes();
    auto is = buildInstrumentSet(quotes);
    PillarOptions po;
    po.interpolation = rl::market::CurveInterpolation::MonotoneConvex;
    po.perCurve["USD-FEDFUNDS"] = rl::market::CurveInterpolation::LinearLogDf;
    auto specs = makeCurveSpecs(is, quotes.asOf(), po);
    REQUIRE(specs.size() == 2);
    CHECK(specs[0].name == "USD-SOFR");
    CHECK(specs[0].isDiscount);
    CHECK(specs[0].projects == std::vector<std::string>{"USD-SOFR"});
    CHECK(specs[1].projects == std::vector<std::string>{"USD-FEDFUNDS"});
    CHECK(specs[0].interpolation == rl::market::CurveInterpolation::MonotoneConvex);
    CHECK(specs[1].name == "USD-FEDFUNDS");
    CHECK_FALSE(specs[1].isDiscount);
    CHECK(specs[1].interpolation == rl::market::CurveInterpolation::LinearLogDf);
    CHECK(parameterCount(specs) == is.size());
    for (const auto& s : specs) {
        REQUIRE(s.pillarTimes.size() == s.pillarDates.size());
        REQUIRE(s.pillarIds.size() == s.pillarDates.size());
        CHECK(s.turns.empty());
        for (std::size_t i = 0; i < s.pillarTimes.size(); ++i) {
            CHECK_THAT(s.pillarTimes[i], WithinAbs(rl::core::yearFraction(rl::core::DayCount::Act365Fixed, quotes.asOf(), s.pillarDates[i]), 1e-15));
            if (i > 0) CHECK(s.pillarTimes[i] > s.pillarTimes[i - 1]);
        }
        CHECK(s.pillarTimes.front() > 0.0);
    }
    CHECK(specs[0].pillarIds.front() == "SOFR-OIS-1W");
    CHECK(specs[0].pillarIds.back() == "FF-SOFR-30Y");
    CHECK(specs[0].pillarTimes.back() > 29.5);
    CHECK(specs[0].pillarTimes.back() < 30.5);
}

TEST_CASE("Curve specs: turns become jump parameters, not nodes; missing discount curve is an error", "[bootstrap][pillars]") {
    auto quotes = usdQuotes();
    InstrumentSetOptions io;
    io.includeTurns = true;
    auto is = buildInstrumentSet(quotes, io);
    auto specs = makeCurveSpecs(is, quotes.asOf());
    CHECK(specs[0].turns.size() == 1);
    CHECK(specs[0].pillarTimes.size() + 1 == specs[0].parameterCount());
    CHECK(parameterCount(specs) == is.size());
    PillarOptions bad;
    bad.discountCurve = "EUR-ESTR";
    CHECK_THROWS_AS(makeCurveSpecs(is, quotes.asOf(), bad), rl::InputError);
}
