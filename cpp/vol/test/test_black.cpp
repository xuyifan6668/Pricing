// Shifted Black-76: parity, a published reference value, shift identity, degenerate strikes, greeks, implied vol.
#include "rl/vol/black.hpp"

#include "rl/aad/real.hpp"
#include "rl/core/error.hpp"

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>

using namespace rl::vol;
using rl::aad::Real;
using rl::math::value;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {
constexpr OptionType kTypes[] = {OptionType::Call, OptionType::Put};
}

TEMPLATE_TEST_CASE("Black: put-call parity with and without a shift", "[vol][black]", double, Real) {
    using T = TestType;
    for (double shift : {0.0, 0.01, 0.03}) {
        for (double f : {0.035, -0.004}) {
            if (f + shift <= 0.0) continue;
            for (double k : {0.02, 0.035, 0.001, 0.08}) {
                for (double tt : {0.25, 5.0}) {
                    const double c = value(blackPrice(T(f), k, tt, T(0.3), shift, OptionType::Call));
                    const double p = value(blackPrice(T(f), k, tt, T(0.3), shift, OptionType::Put));
                    CHECK_THAT(c - p, WithinAbs(f - k, 1e-15));
                }
            }
        }
    }
    rl::aad::Tape::current().reset();
}

TEST_CASE("Black: Haug (2007) Black-76 reference value", "[vol][black]") {
    // F = 19, K = 19, T = 0.75, r = 10%, σ = 28%: discounted call = put = 1.7011 (Haug 2e, §1.1.4).
    const double disc = std::exp(-0.10 * 0.75);
    CHECK_THAT(disc * blackPrice(19.0, 19.0, 0.75, 0.28, 0.0, OptionType::Call), WithinAbs(1.7011, 5e-5));
    CHECK_THAT(disc * blackPrice(19.0, 19.0, 0.75, 0.28, 0.0, OptionType::Put), WithinAbs(1.7011, 5e-5));
}

TEST_CASE("Black: a shift is exactly a displaced forward and strike", "[vol][black]") {
    for (OptionType type : kTypes) {
        CHECK(blackPrice(0.002, -0.003, 3.0, 0.4, 0.02, type) == blackPrice(0.002 + 0.02, -0.003 + 0.02, 3.0, 0.4, 0.0, type));
        CHECK(blackVega(0.002, -0.003, 3.0, 0.4, 0.02) == blackVega(0.022, 0.017, 3.0, 0.4, 0.0));
    }
}

TEST_CASE("Black: strike below minus the shift, forward below minus the shift, zero vol", "[vol][black]") {
    // K + s ≤ 0 < F + s: the call is exercised almost surely, the put never.
    CHECK_THAT(blackPrice(0.01, -0.02, 1.0, 0.3, 0.01, OptionType::Call), WithinAbs(0.03, 1e-17));
    CHECK(blackPrice(0.01, -0.02, 1.0, 0.3, 0.01, OptionType::Put) == 0.0);
    CHECK(blackDelta(0.01, -0.02, 1.0, 0.3, 0.01, OptionType::Call) == 1.0);
    CHECK(blackVega(0.01, -0.02, 1.0, 0.3, 0.01) == 0.0);
    CHECK_THROWS_AS(blackPrice(-0.02, 0.01, 1.0, 0.3, 0.01, OptionType::Call), rl::InputError);
    CHECK_THAT(blackPrice(0.04, 0.03, 1.0, 0.0, 0.0, OptionType::Call), WithinAbs(0.01, 1e-17));
    CHECK(blackPrice(0.04, 0.03, 1.0, 0.0, 0.0, OptionType::Put) == 0.0);
}

TEST_CASE("Black: greeks match central finite differences", "[vol][black][greeks]") {
    for (double shift : {0.0, 0.02}) {
        for (double k : {0.02, 0.035, 0.06}) {
            for (double tt : {0.25, 2.0, 10.0}) {
                for (OptionType type : kTypes) {
                    const double f = 0.035, s = 0.25, h = 1e-7;
                    const double vegaFd = (blackPrice(f, k, tt, s + h, shift, type) - blackPrice(f, k, tt, s - h, shift, type)) / (2.0 * h);
                    const double deltaFd = (blackPrice(f + h, k, tt, s, shift, type) - blackPrice(f - h, k, tt, s, shift, type)) / (2.0 * h);
                    const double gammaFd = (blackDelta(f + 1e-6, k, tt, s, shift, type) - blackDelta(f - 1e-6, k, tt, s, shift, type)) / 2e-6;
                    CHECK(std::fabs(blackVega(f, k, tt, s, shift) - vegaFd) <= 1e-7 * (1.0 + std::fabs(vegaFd)));
                    CHECK(std::fabs(blackDelta(f, k, tt, s, shift, type) - deltaFd) <= 1e-7 * (1.0 + std::fabs(deltaFd)));
                    CHECK(std::fabs(blackGamma(f, k, tt, s, shift) - gammaFd) <= 1e-5 * (1.0 + std::fabs(gammaFd)));
                }
            }
        }
    }
}

TEST_CASE("Black implied vol round-trips on OTM options to 1e-12", "[vol][black][implied]") {
    const double fwd = 0.03;
    for (double shift : {0.0, 0.02}) {
        for (double s : {0.10, 0.25, 0.60}) {
            for (double tt : {0.25, 1.0, 10.0}) {
                const double f = fwd + shift, sd = s * std::sqrt(tt);
                for (int x = -3; x <= 3; ++x) {
                    const double k = f * std::exp(x * sd);
                    const double strike = k - shift;
                    const OptionType otm = k >= f ? OptionType::Call : OptionType::Put;
                    INFO("shift " << shift << " sigma " << s << " T " << tt << " x " << x);
                    CHECK_THAT(blackImpliedVol(blackPrice(fwd, strike, tt, s, shift, otm), fwd, strike, tt, shift, otm), WithinRel(s, 1e-12));
                    if (std::abs(x) <= 1) {
                        const OptionType itm = otm == OptionType::Call ? OptionType::Put : OptionType::Call;
                        CHECK_THAT(blackImpliedVol(blackPrice(fwd, strike, tt, s, shift, itm), fwd, strike, tt, shift, itm), WithinRel(s, 1e-11));
                    }
                }
            }
        }
    }
}

TEST_CASE("Black implied vol: domain errors and bounds", "[vol][black][implied]") {
    CHECK(blackImpliedVol(0.01, 0.04, 0.03, 1.0, 0.0, OptionType::Call) == 0.0);
    CHECK_THROWS_AS(blackImpliedVol(0.005, 0.04, 0.03, 1.0, 0.0, OptionType::Call), rl::InputError);    // below intrinsic
    CHECK_THROWS_AS(blackImpliedVol(0.031, 0.03, 0.04, 1.0, 0.0, OptionType::Call), rl::InputError);    // OTM call ≥ F
    CHECK_THROWS_AS(blackImpliedVol(0.01, 0.03, -0.01, 1.0, 0.0, OptionType::Call), rl::InputError);    // K + s ≤ 0
    CHECK_THROWS_AS(blackImpliedVol(0.01, -0.01, 0.01, 1.0, 0.0, OptionType::Call), rl::InputError);    // F + s ≤ 0
    CHECK_THROWS_AS(blackImpliedVol(0.01, 0.03, 0.03, 0.0, 0.0, OptionType::Call), rl::InputError);     // T = 0
}

TEST_CASE("Black under Real: bit-identical to double and active", "[vol][black][aad]") {
    const Real f = Real::input(0.035), s = Real::input(0.3);
    for (OptionType type : kTypes) {
        const Real p = blackPrice(f, 0.03, 2.0, s, 0.01, type);
        CHECK(p.value() == blackPrice(0.035, 0.03, 2.0, 0.3, 0.01, type));
        CHECK(p.isActive());
        CHECK(blackDelta(f, 0.03, 2.0, s, 0.01, type).value() == blackDelta(0.035, 0.03, 2.0, 0.3, 0.01, type));
    }
    CHECK(blackVega(f, 0.03, 2.0, s, 0.01).value() == blackVega(0.035, 0.03, 2.0, 0.3, 0.01));
    CHECK(blackGamma(f, 0.03, 2.0, s, 0.01).value() == blackGamma(0.035, 0.03, 2.0, 0.3, 0.01));
    rl::aad::Tape::current().reset();
}
