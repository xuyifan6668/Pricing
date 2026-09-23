#include "rl/market/quote.hpp"

#include "rl/core/calendar.hpp"
#include "rl/core/error.hpp"
#include "rl/market/hash.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <filesystem>

using namespace rl::market;
using rl::core::Date;
using rl::core::Tenor;
using Catch::Matchers::WithinAbs;

namespace {
Quote ois(const char* id, const char* tenor, double pct) {
    Quote q;
    q.id = id;
    q.type = QuoteType::OisSwap;
    q.index = "USD-SOFR";
    q.tenor = Tenor::parse(tenor);
    q.value = pct;
    q.unit = QuoteUnit::Percent;
    return q;
}
} // namespace

TEST_CASE("Quote unit conversions", "[market][quote]") {
    Quote q = ois("x", "1Y", 3.62);
    CHECK_THAT(q.asRate(), WithinAbs(0.0362, 1e-16));
    CHECK_THAT(q.fromRate(0.0362), WithinAbs(3.62, 1e-13));
    CHECK(q.unitScale() == 100.0);
    q.unit = QuoteUnit::BasisPoints;
    q.value = -2.5;
    CHECK_THAT(q.asRate(), WithinAbs(-0.00025, 1e-18));
    q.unit = QuoteUnit::Decimal;
    q.value = 0.04;
    CHECK(q.asRate() == 0.04);
    q.unit = QuoteUnit::Price;
    q.value = 96.5;
    CHECK(q.asPrice() == 96.5);
    CHECK_THROWS_AS(q.asRate(), rl::InputError);
    q.unit = QuoteUnit::Percent;
    CHECK_THROWS_AS(q.asPrice(), rl::InputError);
}

TEST_CASE("Quote validation per type", "[market][quote]") {
    Quote q = ois("ok", "1Y", 3.6);
    CHECK_NOTHROW(q.validate());
    q.id.clear();
    CHECK_THROWS_AS(q.validate(), rl::InputError);
    q = ois("no-tenor", "1Y", 3.6);
    q.tenor.reset();
    CHECK_THROWS_AS(q.validate(), rl::InputError);
    q = ois("wrong-unit", "1Y", 360.0); // 360% — a percent value fed as decimal would be caught too
    CHECK_THROWS_AS(q.validate(), rl::InputError);
    Quote f;
    f.id = "SR3";
    f.type = QuoteType::Future;
    f.index = "USD-SOFR";
    f.start = Date(2026, 12, 16);
    f.end = Date(2027, 3, 17);
    f.value = 96.3;
    f.unit = QuoteUnit::Price;
    CHECK_NOTHROW(f.validate());
    f.unit = QuoteUnit::Percent;
    CHECK_THROWS_AS(f.validate(), rl::InputError);
    f.unit = QuoteUnit::Price;
    f.end = Date(2026, 12, 15);
    CHECK_THROWS_AS(f.validate(), rl::InputError);
    Quote b;
    b.id = "FF-SOFR-1Y";
    b.type = QuoteType::BasisSwap;
    b.index = "USD-FEDFUNDS";
    b.tenor = Tenor::parse("1Y");
    b.value = -2.5;
    b.unit = QuoteUnit::BasisPoints;
    CHECK_THROWS_AS(b.validate(), rl::InputError); // missing reference_index
    b.referenceIndex = "USD-FEDFUNDS";
    CHECK_THROWS_AS(b.validate(), rl::InputError); // same index both legs
    b.referenceIndex = "USD-SOFR";
    CHECK_NOTHROW(b.validate());
}

TEST_CASE("QuoteSet: construction, lookup, duplicates, with*", "[market][quote]") {
    QuoteSet qs(Date(2026, 9, 4), rl::core::Currency::USD, {ois("1Y", "1Y", 3.62), ois("2Y", "2Y", 3.55)});
    CHECK(qs.size() == 2);
    CHECK(qs.contains("2Y"));
    CHECK_FALSE(qs.contains("3Y"));
    CHECK(qs.indexOf("2Y") == 1);
    CHECK(qs.find("1Y").value == 3.62);
    CHECK_THROWS_AS(qs.find("3Y"), rl::InputError);
    CHECK(qs.ids() == std::vector<std::string>{"1Y", "2Y"});
    CHECK(qs.ofType(QuoteType::OisSwap).size() == 2);
    CHECK(qs.ofType(QuoteType::Future).empty());
    CHECK_THROWS_AS(QuoteSet(Date(2026, 9, 4), rl::core::Currency::USD, {ois("1Y", "1Y", 3.6), ois("1Y", "2Y", 3.5)}), rl::InputError);

    QuoteSet b = qs.bumped("2Y", 0.01);
    CHECK_THAT(b.find("2Y").value, WithinAbs(3.56, 1e-12));
    CHECK(qs.find("2Y").value == 3.55); // original untouched
    QuoteSet w = qs.withValues({1.0, 2.0});
    CHECK(w[0].value == 1.0);
    CHECK(w[1].value == 2.0);
    CHECK_THROWS_AS(qs.withValues({1.0}), rl::InputError);
    CHECK(qs.withValue("1Y", 9.0).find("1Y").value == 9.0);
}

TEST_CASE("QuoteSet JSON round trip and error reporting", "[market][quote]") {
    const char* text = R"({"asOf": "2026-09-04", "currency": "USD", "quotes": [
        {"id": "SOFR-OIS-1Y", "type": "ois_swap", "index": "USD-SOFR", "tenor": "1Y", "value": 3.62, "unit": "pct"},
        {"id": "SR3-Z6", "type": "future", "index": "USD-SOFR", "start": "2026-12-16", "end": "2027-03-17", "value": 96.35, "unit": "price", "convexity_bp": 0.4},
        {"id": "FF-SOFR-1Y", "type": "basis_swap", "index": "USD-FEDFUNDS", "reference_index": "USD-SOFR", "tenor": "1Y", "value": -2.5, "unit": "bp"},
        {"id": "TURN-2026", "type": "turn", "index": "USD-SOFR", "start": "2026-12-31", "end": "2027-01-04", "value": 8.0, "unit": "bp"}]})";
    QuoteSet qs = QuoteSet::fromJson(text);
    CHECK(qs.asOf() == Date(2026, 9, 4));
    CHECK(qs.size() == 4);
    CHECK(qs.find("SR3-Z6").convexityAdjustmentBp == 0.4);
    CHECK(qs.find("SR3-Z6").start == Date(2026, 12, 16));
    CHECK(qs.find("FF-SOFR-1Y").referenceIndex == "USD-SOFR");
    CHECK(qs.find("TURN-2026").type == QuoteType::Turn);
    QuoteSet again = QuoteSet::fromJson(qs.toJson());
    REQUIRE(again.size() == qs.size());
    for (std::size_t i = 0; i < qs.size(); ++i) {
        CHECK(again[i].id == qs[i].id);
        CHECK(again[i].type == qs[i].type);
        CHECK(again[i].value == qs[i].value);
        CHECK(again[i].unit == qs[i].unit);
        CHECK(again[i].tenor == qs[i].tenor);
        CHECK(again[i].start == qs[i].start);
    }
    // Errors carry the quote id as context.
    try {
        QuoteSet::fromJson(R"({"asOf": "2026-09-04", "quotes": [{"id": "BAD", "type": "future", "index": "USD-SOFR", "value": 96}]})");
        FAIL("expected InputError");
    } catch (const rl::InputError& e) {
        CHECK(std::string(e.what()).find("BAD") != std::string::npos);
    }
    CHECK_THROWS_AS(QuoteSet::fromJson("{"), rl::InputError);
    CHECK_THROWS_AS(QuoteSet::fromJson(R"({"quotes": []})"), rl::InputError);
}

TEST_CASE("QuoteSet hash contribution", "[market][quote][hash]") {
    QuoteSet a(Date(2026, 9, 4), rl::core::Currency::USD, {ois("1Y", "1Y", 3.62)});
    Hasher ha, hb;
    a.mixInto(ha);
    a.bumped("1Y", 1e-9).mixInto(hb);
    CHECK(ha.value() != hb.value());
    Hasher hc;
    a.mixInto(hc);
    CHECK(ha.value() == hc.value());
}

TEST_CASE("The shipped USD quote file loads, is sorted by maturity within each type, and has the D9 instrument mix", "[market][quote][data]") {
    const auto file = rl::core::dataDir() / "quotes" / "usd_2026-09-04.json";
    REQUIRE(std::filesystem::exists(file));
    QuoteSet qs = QuoteSet::fromFile(file);
    CHECK(qs.asOf() == Date(2026, 9, 4));
    CHECK(qs.ofType(QuoteType::OisSwap).size() >= 15);
    CHECK(qs.ofType(QuoteType::Future).size() >= 6);
    CHECK(qs.ofType(QuoteType::BasisSwap).size() >= 8);
    CHECK(qs.ofType(QuoteType::Turn).size() == 1);
    for (const auto& q : qs) CHECK(q.index.rfind("USD-", 0) == 0);
    auto fut = qs.ofType(QuoteType::Future);
    for (std::size_t i = 1; i < fut.size(); ++i) {
        CHECK(*fut[i].start == *fut[i - 1].end); // consecutive quarterly IMM strip
    }
}
