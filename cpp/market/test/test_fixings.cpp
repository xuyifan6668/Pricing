#include "rl/market/fixings.hpp"

#include "rl/core/calendar.hpp"
#include "rl/core/error.hpp"
#include "rl/market/hash.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <filesystem>

using namespace rl::market;
using rl::core::Date;
using Catch::Matchers::WithinAbs;

TEST_CASE("FixingStore: add, lookup, range, overwrite", "[market][fixings]") {
    FixingStore fs;
    fs.add("USD-SOFR", Date(2026, 9, 1), 0.0431);
    fs.add("USD-SOFR", Date(2026, 9, 2), 0.0432);
    fs.add("USD-SOFR", Date(2026, 9, 3), 0.0430);
    fs.add("USD-FEDFUNDS", Date(2026, 9, 3), 0.0433);
    CHECK(fs.size() == 4);
    CHECK(fs.has("USD-SOFR"));
    CHECK_FALSE(fs.has("EUR-ESTR"));
    CHECK(fs.fixing("USD-SOFR", Date(2026, 9, 2)) == 0.0432);
    CHECK_FALSE(fs.fixing("USD-SOFR", Date(2026, 9, 4)).has_value());
    CHECK_FALSE(fs.fixing("GBP-SONIA", Date(2026, 9, 2)).has_value());
    CHECK(fs.require("USD-FEDFUNDS", Date(2026, 9, 3)) == 0.0433);
    CHECK_THROWS_AS(fs.require("USD-SOFR", Date(2026, 9, 4)), rl::InputError);
    CHECK(fs.firstDate("USD-SOFR") == Date(2026, 9, 1));
    CHECK(fs.lastDate("USD-SOFR") == Date(2026, 9, 3));
    CHECK_FALSE(fs.lastDate("JPY-TONA").has_value());
    auto r = fs.range("USD-SOFR", Date(2026, 9, 2), Date(2026, 9, 30));
    REQUIRE(r.size() == 2);
    CHECK(r[0].first == Date(2026, 9, 2));
    CHECK(r[1].second == 0.0430);
    fs.add("USD-SOFR", Date(2026, 9, 2), 0.0435); // correction overwrites
    CHECK(fs.fixing("USD-SOFR", Date(2026, 9, 2)) == 0.0435);
    CHECK(fs.size() == 4);
    CHECK_THROWS_AS(fs.add("", Date(2026, 9, 2), 0.01), rl::InputError);
    CHECK_THROWS_AS(fs.add("USD-SOFR", Date(2026, 9, 2), std::numeric_limits<double>::quiet_NaN()), rl::InputError);
    CHECK(fs.indices() == std::vector<std::string>{"USD-FEDFUNDS", "USD-SOFR"});
}

TEST_CASE("FixingStore JSON round trip (percent in the file, decimal in memory)", "[market][fixings]") {
    const char* text = R"({"USD-SOFR": {"2026-09-01": 4.31, "2026-09-02": 4.32}, "USD-FEDFUNDS": {"2026-09-02": 4.33}})";
    FixingStore fs = FixingStore::fromJson(text);
    CHECK(fs.size() == 3);
    CHECK_THAT(fs.require("USD-SOFR", Date(2026, 9, 1)), WithinAbs(0.0431, 1e-16));
    FixingStore again = FixingStore::fromJson(fs.toJson());
    CHECK(again == fs);
    CHECK_THROWS_AS(FixingStore::fromJson("[1,2]"), rl::InputError);
    CHECK_THROWS_AS(FixingStore::fromJson(R"({"USD-SOFR": {"2026-09-01": "x"}})"), rl::InputError);
    CHECK_THROWS_AS(FixingStore::fromJson("{not json"), rl::InputError);
}

TEST_CASE("FixingStore hash contribution is deterministic and value-sensitive", "[market][fixings][hash]") {
    FixingStore a, b;
    a.add("USD-SOFR", Date(2026, 9, 1), 0.0431);
    b.add("USD-SOFR", Date(2026, 9, 1), 0.0431);
    Hasher ha, hb;
    a.mixInto(ha);
    b.mixInto(hb);
    CHECK(ha.value() == hb.value());
    b.add("USD-SOFR", Date(2026, 9, 1), 0.0431 + 1e-12);
    Hasher hc;
    b.mixInto(hc);
    CHECK(hc.value() != ha.value());
}

TEST_CASE("The shipped USD fixings file loads and covers the ten business days before as-of", "[market][fixings][data]") {
    const auto file = rl::core::dataDir() / "fixings" / "usd_2026-09-04.json";
    REQUIRE(std::filesystem::exists(file));
    FixingStore fs = FixingStore::fromFile(file);
    CHECK(fs.has("USD-SOFR"));
    CHECK(fs.has("USD-FEDFUNDS"));
    // Every SIFMA business day in [2026-08-21, 2026-09-03] has a SOFR print.
    auto cal = rl::core::Calendar::get("USD-SIFMA");
    for (Date d(2026, 8, 21); d <= Date(2026, 9, 3); ++d)
        if (cal.isBusinessDay(d)) CHECK(fs.fixing("USD-SOFR", d).has_value());
    CHECK_FALSE(fs.fixing("USD-SOFR", Date(2026, 9, 4)).has_value()); // the as-of day's fixing is not known yet
}
