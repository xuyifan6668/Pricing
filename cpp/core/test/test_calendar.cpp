#include "rl/core/calendar.hpp"

#include "rl/core/error.hpp"

#include <catch2/catch_test_macros.hpp>

#include <filesystem>
#include <fstream>
#include <string>

using namespace rl::core;
using BDC = BusinessDayConvention;

// Holiday facts below are independently known (Federal Reserve K.8, SIFMA holiday schedule,
// gov.uk/bank-holidays, ECB TARGET2 calendar, Japanese Cabinet Office) and also pinned in
// tools/test_gen_holidays.py, so a failure here separates "loader bug" from "data bug".

TEST_CASE("Registry loads every calendar", "[core][calendar]") {
    for (auto name : {"USD-SIFMA", "USD-FEDWIRE", "GBP", "EUR-TARGET2", "JPY"}) {
        Calendar c = Calendar::get(name);
        CHECK(c.name() == name);
        CHECK(c.holidays().size() > 400); // TARGET2 has 6 closures a year; the files never list weekends
        CHECK(c.firstHoliday().year() == 2000);
        CHECK(c.lastHoliday().year() == 2075);
    }
    CHECK(Calendar::get("NONE").holidays().empty());
    CHECK_THROWS_AS(Calendar::get("MARS"), rl::InputError);
}

TEST_CASE("Holiday table is sorted and unique", "[core][calendar]") {
    const Calendar jpy = Calendar::get("JPY");
    const auto& h = jpy.holidays();
    for (std::size_t i = 1; i < h.size(); ++i) CHECK(h[i] > h[i - 1]);
}

TEST_CASE("Known US holidays 2026: SIFMA vs Fedwire", "[core][calendar]") {
    Calendar sifma = Calendar::get("USD-SIFMA"), fed = Calendar::get("USD-FEDWIRE");
    for (auto d : {Date(2026, 1, 1), Date(2026, 1, 19), Date(2026, 2, 16), Date(2026, 5, 25), Date(2026, 6, 19),
                   Date(2026, 9, 7), Date(2026, 10, 12), Date(2026, 11, 11), Date(2026, 11, 26), Date(2026, 12, 25)}) {
        CHECK(sifma.isHoliday(d));
        CHECK(fed.isHoliday(d));
    }
    // 4 July 2026 is a Saturday: SIFMA closes Friday 3 July; Fedwire is open.
    CHECK(sifma.isHoliday(Date(2026, 7, 3)));
    CHECK_FALSE(fed.isHoliday(Date(2026, 7, 3)));
    CHECK(fed.isBusinessDay(Date(2026, 7, 3)));
    // Good Friday is a bond-market close, not a Fed holiday.
    CHECK(sifma.isHoliday(Date(2026, 4, 3)));
    CHECK_FALSE(fed.isHoliday(Date(2026, 4, 3)));
}

TEST_CASE("Known GBP, TARGET2 and JPY holidays", "[core][calendar]") {
    Calendar gbp = Calendar::get("GBP");
    CHECK(gbp.isHoliday(Date(2026, 4, 6)));   // Easter Monday
    CHECK(gbp.isHoliday(Date(2026, 8, 31)));  // Summer bank holiday
    CHECK(gbp.isHoliday(Date(2026, 12, 28))); // Boxing Day substitute (26th is a Saturday)
    CHECK(gbp.isHoliday(Date(2023, 5, 8)));   // Coronation
    CHECK_FALSE(gbp.isHoliday(Date(2026, 5, 1)));
    Calendar eur = Calendar::get("EUR-TARGET2");
    CHECK(eur.isHoliday(Date(2026, 5, 1)));
    CHECK(eur.isHoliday(Date(2026, 4, 3)));
    CHECK_FALSE(eur.isHoliday(Date(2026, 8, 31)));
    Calendar jpy = Calendar::get("JPY");
    CHECK(jpy.isHoliday(Date(2026, 1, 2)));   // bank holiday
    CHECK(jpy.isHoliday(Date(2026, 3, 20)));  // Vernal Equinox
    CHECK(jpy.isHoliday(Date(2026, 5, 6)));   // substitute for Constitution Day on Sunday
    CHECK(jpy.isHoliday(Date(2026, 9, 22)));  // Citizens' holiday between Respect for the Aged and Autumnal Equinox
    CHECK(jpy.isHoliday(Date(2026, 12, 31)));
}

TEST_CASE("Weekends are not holidays but are not business days", "[core][calendar]") {
    Calendar c = Calendar::get("USD-SIFMA");
    CHECK_FALSE(c.isHoliday(Date(2026, 9, 5)));
    CHECK(c.isWeekend(Date(2026, 9, 5)));
    CHECK_FALSE(c.isBusinessDay(Date(2026, 9, 5)));
    CHECK(c.isBusinessDay(Date(2026, 9, 4)));
}

TEST_CASE("fromHolidays, fromFile parse, join", "[core][calendar]") {
    Calendar a = Calendar::fromHolidays("A", {Date(2026, 9, 9), Date(2026, 9, 9), Date(2026, 9, 1)});
    CHECK(a.holidays().size() == 2);
    CHECK(a.holidays()[0] == Date(2026, 9, 1));
    Calendar b = Calendar::fromHolidays("B", {Date(2026, 9, 10)});
    Calendar j = Calendar::join(a, b);
    CHECK(j.name() == "A+B");
    CHECK(j.holidays().size() == 3);
    CHECK(j.isHoliday(Date(2026, 9, 10)));
    Calendar reg = Calendar::get("USD-SIFMA+GBP");
    CHECK(reg.isHoliday(Date(2026, 8, 31)));
    CHECK(reg.isHoliday(Date(2026, 9, 7)));
}

TEST_CASE("Bad holiday file lines carry file:line context", "[core][calendar]") {
    auto p = std::filesystem::temp_directory_path() / "rl_bad_holidays.txt";
    {
        std::ofstream(p) << "# comment\n2026-01-01\n\n2026-02-30  # nope\n";
    }
    try {
        Calendar::fromFile("BAD", p);
        FAIL("should throw");
    } catch (const rl::Error& e) {
        REQUIRE(e.context().size() == 1);
        CHECK(e.context()[0].find(":4") != std::string::npos);
    }
    std::filesystem::remove(p);
}

TEST_CASE("Business-day convention strings", "[core][calendar]") {
    CHECK(bdcFromString("MF") == BDC::ModifiedFollowing);
    CHECK(bdcFromString("Modified Following") == BDC::ModifiedFollowing);
    CHECK(bdcFromString("F") == BDC::Following);
    CHECK(bdcFromString("P") == BDC::Preceding);
    CHECK(bdcFromString("Unadjusted") == BDC::Unadjusted);
    CHECK_THROWS_AS(bdcFromString("Sideways"), rl::InputError);
}

// --- the maths -------------------------------------------------------------------------------

TEST_CASE("adjust", "[core][calendar]") {
    Calendar sifma = Calendar::get("USD-SIFMA"), fed = Calendar::get("USD-FEDWIRE");
    Date sun(2027, 1, 31); // Sunday
    CHECK(sifma.adjust(sun, BDC::Unadjusted) == sun);
    CHECK(sifma.adjust(sun, BDC::Following) == Date(2027, 2, 1));
    CHECK(sifma.adjust(sun, BDC::ModifiedFollowing) == Date(2027, 1, 29));
    CHECK(sifma.adjust(sun, BDC::Preceding) == Date(2027, 1, 29));
    CHECK(sifma.adjust(sun, BDC::ModifiedPreceding) == Date(2027, 1, 29));
    Date sat(2026, 8, 1); // Saturday, start of month
    CHECK(sifma.adjust(sat, BDC::Preceding) == Date(2026, 7, 31));
    CHECK(sifma.adjust(sat, BDC::ModifiedPreceding) == Date(2026, 8, 3));
    CHECK(sifma.adjust(sat, BDC::ModifiedFollowing) == Date(2026, 8, 3));
    // Business days are fixed points.
    CHECK(sifma.adjust(Date(2026, 9, 4), BDC::ModifiedFollowing) == Date(2026, 9, 4));
    // Holiday plus weekend run: Fri 3 Jul 2026 (SIFMA) → Mon 6 Jul; Fedwire keeps the Friday.
    CHECK(sifma.adjust(Date(2026, 7, 3), BDC::Following) == Date(2026, 7, 6));
    CHECK(fed.adjust(Date(2026, 7, 3), BDC::Following) == Date(2026, 7, 3));
    // Golden Week 2026 (Tokyo): 2026-05-03 (Sun) → next business day is Thu 7 May.
    CHECK(Calendar::get("JPY").adjust(Date(2026, 5, 3), BDC::Following) == Date(2026, 5, 7));
}

TEST_CASE("advance by business days", "[core][calendar]") {
    Calendar c = Calendar::get("USD-SIFMA");
    CHECK(c.advance(Date(2026, 9, 4), 1) == Date(2026, 9, 8)); // over the weekend and Labor Day
    CHECK(c.advance(Date(2026, 9, 4), 2) == Date(2026, 9, 9));
    CHECK(c.advance(Date(2026, 9, 8), -1) == Date(2026, 9, 4));
    CHECK(c.advance(Date(2026, 9, 8), 0) == Date(2026, 9, 8));
    CHECK(c.advance(Date(2026, 9, 5), 0) == Date(2026, 9, 5)); // zero leaves a Saturday alone
    CHECK(c.advance(Date(2026, 9, 5), 1) == Date(2026, 9, 8));
    CHECK(c.advance(Date(2026, 9, 2), 2) == Date(2026, 9, 4)); // USD spot from Wed 2 Sep
}

TEST_CASE("advance by tenor with and without EOM", "[core][calendar]") {
    Calendar c = Calendar::get("USD-SIFMA");
    // 27 Feb 2026 is the last business day of February.
    CHECK(c.advance(Date(2026, 2, 27), Tenor::months(6), BDC::ModifiedFollowing, true) == Date(2026, 8, 31));
    CHECK(c.advance(Date(2026, 2, 27), Tenor::months(6), BDC::ModifiedFollowing, false) == Date(2026, 8, 27));
    // Not the last business day: EOM flag has no effect.
    CHECK(c.advance(Date(2026, 2, 26), Tenor::months(6), BDC::ModifiedFollowing, true) == Date(2026, 8, 26));
    // 1W deposit is 7 calendar days then adjusted; 2026-09-01 (Tue) + 1W = 2026-09-08 (Tue).
    CHECK(c.advance(Date(2026, 9, 1), Tenor::weeks(1), BDC::Following) == Date(2026, 9, 8));
    // 2026-09-04 + 1D = Sat 5 Sep → Following → Tue 8 Sep (Labor Day).
    CHECK(c.advance(Date(2026, 9, 4), Tenor::days(1), BDC::Following) == Date(2026, 9, 8));
    // MF month crossing: 2027-07-31 (Sat) → Fri 30 Jul.
    CHECK(c.advance(Date(2027, 1, 31), Tenor::months(6), BDC::ModifiedFollowing) == Date(2027, 7, 30));
    // Negative tenor.
    CHECK(c.advance(Date(2026, 9, 8), Tenor::months(-1), BDC::Following) == Date(2026, 8, 10)); // 8 Aug is Saturday
}

TEST_CASE("businessDaysBetween", "[core][calendar]") {
    Calendar c = Calendar::get("USD-SIFMA");
    CHECK(c.businessDaysBetween(Date(2026, 9, 4), Date(2026, 9, 8)) == 1);
    CHECK(c.businessDaysBetween(Date(2026, 9, 4), Date(2026, 9, 8), true, true) == 2);
    CHECK(c.businessDaysBetween(Date(2026, 9, 4), Date(2026, 9, 8), false, false) == 0);
    CHECK(c.businessDaysBetween(Date(2026, 9, 8), Date(2026, 9, 4)) == -1);
    CHECK(c.businessDaysBetween(Date(2026, 9, 8), Date(2026, 9, 8)) == 0);
    // A calendar year of SIFMA business days: 2026 has 261 weekdays and 12 SIFMA closes on weekdays
    // (1 Jan, 19 Jan, 16 Feb, 3 Apr, 25 May, 19 Jun, 3 Jul, 7 Sep, 12 Oct, 11 Nov, 26 Nov, 25 Dec).
    CHECK(c.businessDaysBetween(Date(2026, 1, 1), Date(2027, 1, 1)) == 261 - 12);
}
