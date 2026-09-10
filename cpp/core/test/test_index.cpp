#include "rl/core/index.hpp"

#include "rl/core/error.hpp"

#include <catch2/catch_test_macros.hpp>

using namespace rl::core;

TEST_CASE("Currency strings and settlement conventions", "[core][currency]") {
    CHECK(currencyFromString("usd") == Currency::USD);
    CHECK(toString(Currency::JPY) == "JPY");
    CHECK_THROWS_AS(currencyFromString("XAU"), rl::InputError);
    CHECK(settlementDays(Currency::USD) == 2);
    CHECK(settlementDays(Currency::GBP) == 0);
    CHECK(settlementCalendarName(Currency::EUR) == "EUR-TARGET2");
}

TEST_CASE("Standard indices carry market conventions", "[core][index]") {
    Index sofr = Index::sofr();
    CHECK(sofr.name == "USD-SOFR");
    CHECK(sofr.isOvernight());
    CHECK(sofr.dayCount == DayCount::Act360);
    CHECK(sofr.fixingCalendar == "USD-SIFMA");
    CHECK(sofr.compounding == Compounding::CompoundedDaily);
    CHECK(sofr.paymentDelayDays == 2);
    CHECK(sofr.fixingLag == 0);

    Index ff = Index::fedFunds();
    CHECK(ff.compounding == Compounding::AveragedDaily);
    CHECK(ff.fixingCalendar == "USD-FEDWIRE");

    CHECK(Index::sonia().dayCount == DayCount::Act365Fixed);
    CHECK(Index::sonia().paymentDelayDays == 0);
    CHECK(Index::estr().paymentDelayDays == 1);
    CHECK(Index::tona().currency == Currency::JPY);

    Index e6 = Index::euribor6m();
    CHECK_FALSE(e6.isOvernight());
    CHECK(e6.tenor == Tenor::months(6));
    CHECK(e6.fixingLag == 2);
    CHECK(e6.compounding == Compounding::Simple);
    CHECK(Index::termSofr3m().tenor == Tenor::months(3));
}

TEST_CASE("Index registry and equality by name", "[core][index]") {
    for (auto n : {"USD-SOFR", "USD-FEDFUNDS", "EUR-ESTR", "GBP-SONIA", "JPY-TONA", "EUR-EURIBOR-3M", "EUR-EURIBOR-6M", "USD-TERMSOFR-3M"})
        CHECK(Index::get(n).name == n);
    CHECK(Index::get("USD-SOFR") == Index::sofr());
    CHECK(Index::sofr() != Index::fedFunds());
    CHECK_THROWS_AS(Index::get("USD-LIBOR-3M"), rl::InputError);
}

TEST_CASE("Enum names", "[core][index]") {
    CHECK(toString(Compounding::AveragedDaily) == "AveragedDaily");
    CHECK(toString(IndexKind::Term) == "Term");
}
