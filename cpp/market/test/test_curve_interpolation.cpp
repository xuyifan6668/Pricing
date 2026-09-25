#include "rl/market/curve_interpolation.hpp"

#include "rl/core/error.hpp"

#include <catch2/catch_test_macros.hpp>

using namespace rl::market;

TEST_CASE("CurveInterpolation strings round-trip and accept aliases", "[market][interp]") {
    for (auto k : {CurveInterpolation::LinearLogDf, CurveInterpolation::MonotoneConvex})
        CHECK(curveInterpolationFromString(toString(k)) == k);
    CHECK(curveInterpolationFromString("Linear") == CurveInterpolation::LinearLogDf);
    CHECK(curveInterpolationFromString("log-linear") == CurveInterpolation::LinearLogDf);
    CHECK(curveInterpolationFromString("MC") == CurveInterpolation::MonotoneConvex);
    CHECK(curveInterpolationFromString("Monotone Convex") == CurveInterpolation::MonotoneConvex);
    CHECK_THROWS_AS(curveInterpolationFromString("cubic"), rl::InputError);
}
