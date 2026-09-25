// CurveSet<T>: roles, lookup, replacement, validation.
#include "rl/market/curve_set.hpp"

#include "rl/aad/real.hpp"
#include "rl/core/error.hpp"

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>

#include <vector>

using namespace rl::market;
using rl::aad::Real;
using rl::core::Date;
using rl::math::value;

namespace {
const Date kRef(2026, 9, 4);
template <class T>
YieldCurve<T> flat(const char* name, double r, Date ref = kRef) {
    return YieldCurve<T>(ref, name, {1.0, 5.0, 30.0}, {T(r * 1.0), T(r * 5.0), T(r * 30.0)}, CurveInterpolation::LinearLogDf);
}
} // namespace

TEMPLATE_TEST_CASE("CurveSet: roles, lookup, names", "[market][curveset]", double, Real) {
    using T = TestType;
    CurveSet<T> cs({flat<T>("USD-SOFR-OIS", 0.035), flat<T>("USD-FF", 0.0353)}, "USD-SOFR-OIS",
                   {{"USD-SOFR", "USD-SOFR-OIS"}, {"USD-FEDFUNDS", "USD-FF"}});
    CHECK(cs.referenceDate() == kRef);
    CHECK(cs.discount().name() == "USD-SOFR-OIS");
    CHECK(cs.projection("USD-SOFR").name() == "USD-SOFR-OIS");
    CHECK(cs.projection("USD-FEDFUNDS").name() == "USD-FF");
    CHECK(&cs.projection("USD-SOFR") == &cs.discount()); // stored once
    CHECK(cs.hasProjection("USD-FEDFUNDS"));
    CHECK_FALSE(cs.hasProjection("EUR-ESTR"));
    CHECK_THROWS_AS(cs.projection("EUR-ESTR"), rl::InputError);
    CHECK(cs.curve("USD-FF").name() == "USD-FF");
    CHECK_THROWS_AS(cs.curve("nope"), rl::InputError);
    CHECK(cs.names() == std::vector<std::string>{"USD-SOFR-OIS", "USD-FF"});
    CHECK(value(cs.projection("USD-FEDFUNDS").df(1.0)) < value(cs.discount().df(1.0)));
}

TEST_CASE("CurveSet: withCurve replaces one curve and keeps roles", "[market][curveset]") {
    CurveSet<double> cs({flat<double>("D", 0.03), flat<double>("P", 0.04)}, "D", {{"IDX", "P"}, {"DISC-IDX", "D"}});
    CurveSet<double> cs2 = cs.withCurve("D", flat<double>("D", 0.031));
    CHECK(cs2.discount().df(1.0) < cs.discount().df(1.0));
    CHECK(cs2.projection("DISC-IDX").df(1.0) == cs2.discount().df(1.0));
    CHECK(cs2.projection("IDX").df(1.0) == cs.projection("IDX").df(1.0));
    CHECK(cs.discount().df(1.0) == flat<double>("D", 0.03).df(1.0)); // original untouched
    CHECK_THROWS_AS(cs.withCurve("X", flat<double>("X", 0.01)), rl::InputError);
}

TEST_CASE("CurveSet validation", "[market][curveset]") {
    using CS = CurveSet<double>;
    CHECK_THROWS_AS(CS({flat<double>("A", 0.03), flat<double>("A", 0.04)}, "A", {}), rl::InputError);        // duplicate name
    CHECK_THROWS_AS(CS({flat<double>("A", 0.03)}, "B", {}), rl::InputError);                                 // unknown discount
    CHECK_THROWS_AS(CS({flat<double>("A", 0.03)}, "A", {{"IDX", "B"}}), rl::InputError);                     // projection → unknown curve
    CHECK_THROWS_AS(CS({flat<double>("A", 0.03), flat<double>("B", 0.03, Date(2026, 9, 5))}, "A", {}), rl::InputError); // ref dates differ
    CHECK_NOTHROW(CS({flat<double>("A", 0.03)}, "A", {}));                                                   // a discount-only set is fine
}
