// SwaptionVolQuotes: JSON units, grid validation, hashing, and the shipped USD file.
#include "rl/market/vol_quotes.hpp"

#include "rl/core/calendar.hpp" // core::dataDir
#include "rl/core/error.hpp"
#include "rl/market/hash.hpp"
#include "rl/market/swaption_cube.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <string>

using namespace rl::market;
using rl::core::Date;
using rl::core::Tenor;
using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {
std::string nodeJson(const char* e, const char* t, double atmBp = 88.0, const char* smile = R"([{"offset_bp": -50, "vol_bp": 90.5}, {"offset_bp": 50, "vol_bp": 87.1}])") {
    return std::string(R"({"expiry": ")") + e + R"(", "tenor": ")" + t + R"(", "forward_pct": 3.6, "atm_bp": )" + std::to_string(atmBp) +
           R"(, "smile": )" + smile + "}";
}
std::string doc(const std::string& nodes, const char* unit = "normal_bp") {
    return std::string(R"({"asOf": "2026-09-04", "currency": "USD", "index": "USD-SOFR", "unit": ")") + unit + R"(", "nodes": [)" + nodes + "]}";
}
std::string grid2x2() {
    // Deliberately listed out of axis order: 10Y before 2Y, 5Y expiry before 1Y.
    return nodeJson("5Y", "10Y") + "," + nodeJson("1Y", "10Y") + "," + nodeJson("5Y", "2Y") + "," + nodeJson("1Y", "2Y");
}
} // namespace

TEST_CASE("Vol quotes: JSON is in human units, memory is decimal; axes are sorted", "[market][volquotes]") {
    const auto q = SwaptionVolQuotes::fromJson(doc(grid2x2()));
    CHECK(q.asOf() == Date(2026, 9, 4));
    CHECK(q.index() == "USD-SOFR");
    CHECK(q.size() == 4);
    REQUIRE(q.expiries().size() == 2);
    CHECK(q.expiries()[0] == Tenor::years(1));
    CHECK(q.tenors()[0] == Tenor::years(2));
    const auto& n = q.find(Tenor::years(1), Tenor::years(10));
    CHECK(n.id() == "1Yx10Y");
    CHECK_THAT(n.forward, WithinAbs(0.036, 1e-16));
    CHECK_THAT(n.atmNormalVol, WithinAbs(0.0088, 1e-16));
    REQUIRE(n.smile.size() == 2);
    CHECK_THAT(n.smile[0].offset, WithinAbs(-0.005, 1e-16));
    CHECK_THAT(n.smile[1].normalVol, WithinAbs(0.00871, 1e-16));
    CHECK_THROWS_AS(q.find(Tenor::years(2), Tenor::years(10)), rl::InputError);

    const auto back = SwaptionVolQuotes::fromJson(q.toJson());
    CHECK(back.size() == q.size());
    CHECK_THAT(back.find(Tenor::years(5), Tenor::years(2)).smile[1].normalVol, WithinRel(0.00871, 1e-15));
}

TEST_CASE("Vol quotes: validation names the node and rejects what a cube cannot use", "[market][volquotes]") {
    CHECK_THROWS_AS(SwaptionVolQuotes::fromJson(doc(grid2x2() + "," + nodeJson("1Y", "2Y"))), rl::InputError); // duplicate
    CHECK_THROWS_AS(SwaptionVolQuotes::fromJson(doc(nodeJson("5Y", "10Y") + "," + nodeJson("1Y", "10Y") + "," + nodeJson("5Y", "2Y"))),
                    rl::InputError); // hole at 1Yx2Y
    CHECK_THROWS_AS(SwaptionVolQuotes::fromJson(doc(grid2x2(), "lognormal_pct")), rl::InputError);
    CHECK_THROWS_AS(SwaptionVolQuotes::fromJson(doc(nodeJson("1Y", "2Y", 0.0))), rl::InputError);
    CHECK_THROWS_AS(SwaptionVolQuotes::fromJson(doc(nodeJson("1Y", "2Y", 1500.0))), rl::InputError); // a unit error, not a market
    CHECK_THROWS_AS(SwaptionVolQuotes::fromJson(doc(nodeJson("1Y", "2Y", 88.0, R"([{"offset_bp": 0, "vol_bp": 88}])"))), rl::InputError);
    CHECK_THROWS_AS(SwaptionVolQuotes::fromJson(doc(nodeJson("1Y", "2Y", 88.0, R"([{"offset_bp": 50, "vol_bp": 88}, {"offset_bp": -50, "vol_bp": 90}])"))),
                    rl::InputError);
    CHECK_THROWS_AS(SwaptionVolQuotes::fromJson("{\"asOf\": \"2026-09-04\"}"), rl::InputError);
    try {
        SwaptionVolQuotes::fromJson(doc(nodeJson("1Y", "2Y", 88.0, R"([{"offset_bp": -25, "vol_bp": -1}])")));
        FAIL("expected InputError");
    } catch (const rl::InputError& e) {
        CHECK_THAT(std::string(e.what()), ContainsSubstring("1Yx2Y"));
    }
    // A node without a smile is a valid ATM-only node.
    CHECK_NOTHROW(SwaptionVolQuotes::fromJson(doc(nodeJson("1Y", "2Y", 88.0, "[]"))));
}

TEST_CASE("Vol quotes: hash is content-determined and sensitive to one smile vol", "[market][volquotes][hash]") {
    auto hashOf = [](const SwaptionVolQuotes& q) {
        Hasher h;
        q.mixInto(h);
        return h.value();
    };
    const auto a = SwaptionVolQuotes::fromJson(doc(grid2x2()));
    CHECK(hashOf(a) == hashOf(SwaptionVolQuotes::fromJson(doc(grid2x2()))));
    std::string bumped = grid2x2();
    bumped.replace(bumped.rfind("87.1"), 4, "87.2");
    CHECK(hashOf(a) != hashOf(SwaptionVolQuotes::fromJson(doc(bumped))));
}

TEST_CASE("Vol quotes: the shipped USD cube file is a complete 9 x 9 grid in plausible units", "[market][volquotes][data]") {
    const auto q = SwaptionVolQuotes::fromFile(rl::core::dataDir() / "vol" / "usd_swaption_cube_2026-09-04.json");
    CHECK(q.asOf() == Date(2026, 9, 4));
    CHECK(q.index() == "USD-SOFR");
    CHECK(q.expiries().size() == 9);
    CHECK(q.tenors().size() == 9);
    CHECK(q.size() == 81);
    CHECK(q.expiries().front() == Tenor::months(1));
    CHECK(q.tenors().back() == Tenor::years(30));
    for (const auto& n : q.nodes()) {
        INFO(n.id());
        CHECK(n.forward > 0.02);
        CHECK(n.forward < 0.06);
        CHECK(n.atmNormalVol > 0.0040);
        CHECK(n.atmNormalVol < 0.0200);
        REQUIRE(n.smile.size() == 6);
        CHECK_THAT(n.smile.front().offset, WithinAbs(-0.01, 1e-15));
        CHECK_THAT(n.smile.back().offset, WithinAbs(0.01, 1e-15));
        for (const auto& p : n.smile) CHECK(n.forward + p.offset > 0.0); // every strike positive: β = 0.5 with no shift needs it
    }
    CHECK_NOTHROW(q.find(Tenor::years(10), Tenor::years(30))); // the stressed node
}
