#include "rl/core/tenor.hpp"

#include "rl/core/error.hpp"

#include <catch2/catch_test_macros.hpp>

using namespace rl::core;

TEST_CASE("Tenor parsing", "[core][tenor]") {
    CHECK(Tenor::parse("3M") == Tenor::months(3));
    CHECK(Tenor::parse("10y") == Tenor::years(10));
    CHECK(Tenor::parse("1W") == Tenor::weeks(1));
    CHECK(Tenor::parse("7D") == Tenor::weeks(1)); // normalised equality
    CHECK(Tenor::parse(" 6 m ") == Tenor::months(6));
    CHECK(Tenor::parse("ON") == Tenor::days(1));
    CHECK(Tenor::parse("TN") == Tenor::days(2));
    CHECK(Tenor::parse("12M") == Tenor::years(1));
    CHECK(Tenor::parse("-3M") == Tenor::months(-3));
    CHECK_THROWS_AS(Tenor::parse(""), rl::InputError);
    CHECK_THROWS_AS(Tenor::parse("3X"), rl::InputError);
    CHECK_THROWS_AS(Tenor::parse("M"), rl::InputError);
    CHECK_THROWS_AS(Tenor::parse("3.5Y"), rl::InputError);
}

TEST_CASE("Tenor formatting", "[core][tenor]") {
    CHECK(Tenor::months(3).str() == "3M");
    CHECK(Tenor::years(10).str() == "10Y");
    CHECK(Tenor::days(1).str() == "1D");
    CHECK(Tenor::weeks(2).str() == "2W");
}

TEST_CASE("Tenor conversions", "[core][tenor]") {
    CHECK(Tenor::years(2).inMonths() == 24);
    CHECK(Tenor::months(6).inMonths() == 6);
    CHECK(Tenor::weeks(2).inDays() == 14);
    CHECK(Tenor::days(3).inDays() == 3);
    CHECK_THROWS_AS(Tenor::days(3).inMonths(), rl::InputError);
    CHECK_THROWS_AS(Tenor::months(3).inDays(), rl::InputError);
    CHECK(Tenor::months(3).isMonthly());
    CHECK_FALSE(Tenor::weeks(1).isMonthly());
}

TEST_CASE("Tenor arithmetic", "[core][tenor]") {
    CHECK(Tenor::months(3) * 4 == Tenor::years(1));
    CHECK(3 * Tenor::months(2) == Tenor::months(6));
    CHECK(-Tenor::months(3) == Tenor::months(-3));
}
