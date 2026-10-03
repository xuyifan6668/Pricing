#include "rl/bootstrap/instrument_set.hpp"

#include "rl/core/error.hpp"
#include "test_support.hpp"

#include <catch2/catch_test_macros.hpp>

#include <variant>

using namespace rl::bootstrap;
using namespace rl::bootstrap::test;
using rl::core::Date;
using rl::market::QuoteType;

TEST_CASE("Instrument set from the USD file: sorted by curve then maturity, turns dropped by default", "[bootstrap][instruments]") {
    auto quotes = usdQuotes();
    InstrumentSet is = buildInstrumentSet(quotes);
    CHECK(is.size() == quotes.size() - 1); // the turn
    CHECK(curveNames(is) == std::vector<std::string>{"USD-SOFR", "USD-FEDFUNDS"});
    std::string curve;
    Date last = Date::minDate();
    for (const auto& bi : is) {
        if (bi.curveName != curve) {
            curve = bi.curveName;
            last = Date::minDate();
        }
        CHECK(bi.maturity > last);
        last = bi.maturity;
        CHECK(bi.maturity > quotes.asOf());
    }
    // 25Y/30Y basis swaps are assigned to the SOFR curve; 25Y/30Y FF OIS to the FF curve.
    auto sofr = instrumentsOn(is, "USD-SOFR");
    auto ff = instrumentsOn(is, "USD-FEDFUNDS");
    CHECK(sofr.back()->id == "FF-SOFR-30Y");
    CHECK(sofr.back()->type == QuoteType::BasisSwap);
    CHECK(sofr[sofr.size() - 2]->id == "FF-SOFR-25Y");
    CHECK(ff.back()->id == "FF-OIS-30Y");
    CHECK(ff.back()->type == QuoteType::OisSwap);
    CHECK(sofr.size() + ff.size() == is.size());
    // Front of the SOFR curve: OIS 1W first, then the SR3 strip after the 2M OIS.
    CHECK(sofr.front()->id == "SOFR-OIS-1W");
    CHECK(tradeKind(sofr.front()->trade) == "swap");
    std::size_t nFut = 0;
    for (auto* bi : sofr)
        if (bi->type == QuoteType::Future) {
            ++nFut;
            CHECK(std::get<rl::instruments::Future>(bi->trade).end() == bi->maturity);
        }
    CHECK(nFut == 8);
    // Quotes are carried in decimal / price units.
    CHECK(is.front().quote == quotes.find("SOFR-OIS-1W").asRate());
    CHECK(std::get<rl::instruments::Swap>(sofr[3]->trade).id() == sofr[3]->id); // the 2M OIS, before the futures strip
}

TEST_CASE("Instrument set: turns included on request; duplicate pillars rejected", "[bootstrap][instruments]") {
    auto quotes = usdQuotes();
    InstrumentSetOptions io;
    io.includeTurns = true;
    InstrumentSet is = buildInstrumentSet(quotes, io);
    CHECK(is.size() == quotes.size());
    bool turn = false;
    for (const auto& bi : is)
        if (auto* t = std::get_if<TurnOfYear>(&bi.trade)) {
            turn = true;
            CHECK(t->start == Date(2026, 12, 31));
            CHECK(t->jump == 8e-4);
            CHECK(bi.curveName == "USD-SOFR");
        }
    CHECK(turn);
    // Two 10Y SOFR OIS quotes pin the same pillar.
    rl::market::Quote dup = quotes.find("SOFR-OIS-10Y");
    dup.id = "SOFR-OIS-10Y-bis";
    std::vector<rl::market::Quote> qs = quotes.quotes();
    qs.push_back(dup);
    rl::market::QuoteSet bad(quotes.asOf(), quotes.currency(), qs);
    CHECK_THROWS_AS(buildInstrumentSet(bad), rl::InputError);
    // Unsupported index.
    rl::market::Quote estr = quotes.find("SOFR-OIS-10Y");
    estr.id = "ESTR";
    estr.index = "EUR-ESTR";
    CHECK_THROWS_AS(buildInstrumentSet(rl::market::QuoteSet(quotes.asOf(), quotes.currency(), {estr})), rl::InputError);
}

TEST_CASE("Notional option scales the trades", "[bootstrap][instruments]") {
    auto quotes = usdQuotes();
    InstrumentSetOptions io;
    io.notional = 1e6;
    InstrumentSet is = buildInstrumentSet(quotes, io);
    const auto& s = std::get<rl::instruments::Swap>(instrumentsOn(is, "USD-SOFR")[0]->trade);
    CHECK(std::get<rl::instruments::FixedCoupon>(s.fixedLeg().cashflows()[0]).notional == 1e6);
}
