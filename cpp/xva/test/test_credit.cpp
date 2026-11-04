// The credit curve, tested against the par spreads it was built from and the credit triangle.
#include "rl/xva/credit.hpp"

#include "rl/core/error.hpp"
#include "test_support.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>

using namespace rl::xva;
using namespace rl::xva::test;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

TEST_CASE("Hazard curve: survival is the exponential of the integrated hazard rate", "[xva][credit]") {
    const HazardCurve c({1.0, 3.0}, {0.01, 0.03});
    CHECK(c.survival(0.0) == 1.0);
    CHECK_THAT(c.survival(0.5), WithinRel(std::exp(-0.005), 1e-14));
    CHECK_THAT(c.survival(2.0), WithinRel(std::exp(-0.01 - 0.03), 1e-14));
    CHECK_THAT(c.survival(5.0), WithinRel(std::exp(-0.01 - 0.06 - 0.06), 1e-14)); // the last rate continues
    CHECK_THAT(c.defaultProbability(1.0, 3.0), WithinRel(std::exp(-0.01) - std::exp(-0.07), 1e-13));
    CHECK_THAT(HazardCurve::flat(0.02).survival(10.0), WithinRel(std::exp(-0.2), 1e-14));
    CHECK_THROWS_AS(HazardCurve({2.0, 1.0}, {0.01, 0.01}), rl::InputError);
    CHECK_THROWS_AS(HazardCurve({1.0}, {-0.01}), rl::InputError);
}

TEST_CASE("Bootstrap reprices every CDS it was built from", "[xva][credit]") {
    const HazardCurve c = counterpartyCurve();
    REQUIRE(c.times().size() == 5);
    for (const CdsQuote& q : cdsQuotes()) CHECK_THAT(cdsParSpread(c, q.maturityYears, 0.4, discount), WithinAbs(q.spread, 1e-13));
    // Survival falls with time; a rising spread curve means default risk is back-loaded.
    for (double t = 0.5; t < 12.0; t += 0.5) CHECK(c.survival(t) < c.survival(t - 0.5));
    CHECK(c.hazards().back() > 1.5 * c.hazards().front());
}

TEST_CASE("Credit triangle: a flat spread gives a hazard rate close to spread / (1 - recovery)", "[xva][credit]") {
    for (double recovery : {0.0, 0.4}) {
        const double spread = 0.0120;
        const HazardCurve c = HazardCurve::bootstrap({{1.0, spread}, {5.0, spread}, {10.0, spread}}, recovery, discount);
        for (double h : c.hazards()) CHECK_THAT(h, WithinRel(spread / (1.0 - recovery), 0.01));
        CHECK_THAT(c.hazards()[0], WithinRel(c.hazards()[2], 1e-3)); // flat in, flat out
    }
}

TEST_CASE("Wider spreads and lower recovery both mean more default risk priced in", "[xva][credit]") {
    auto quotes = cdsQuotes();
    const double base = counterpartyCurve().survival(10.0);
    for (CdsQuote& q : quotes) q.spread += 0.0010;
    CHECK(HazardCurve::bootstrap(quotes, 0.4, discount).survival(10.0) < base);
    // For the same spreads a higher recovery needs MORE defaults to justify the premium.
    CHECK(HazardCurve::bootstrap(cdsQuotes(), 0.6, discount).survival(10.0) < base);
}

TEST_CASE("Bootstrap rejects curves that cannot be built", "[xva][credit]") {
    CHECK_THROWS_AS(HazardCurve::bootstrap({}, 0.4, discount), rl::InputError);
    CHECK_THROWS_AS(HazardCurve::bootstrap({{5.0, 0.01}, {3.0, 0.01}}, 0.4, discount), rl::InputError);
    CHECK_THROWS_AS(HazardCurve::bootstrap({{5.0, 0.01}}, 1.0, discount), rl::InputError);
    // A 10-year spread far below the 5-year one needs negative default probabilities in between.
    CHECK_THROWS_AS(HazardCurve::bootstrap({{5.0, 0.0300}, {10.0, 0.0100}}, 0.4, discount), rl::CalibrationError);
}
