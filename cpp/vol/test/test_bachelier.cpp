// Bachelier: parity, ATM closed form, greeks against finite differences, degenerate limits, implied vol, Real.
#include "rl/vol/bachelier.hpp"

#include "rl/aad/real.hpp"
#include "rl/core/error.hpp"

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <algorithm>
#include <cmath>

using namespace rl::vol;
using rl::aad::Real;
using rl::math::value;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {
constexpr double kInvSqrt2Pi = 0.39894228040143267794;
constexpr OptionType kTypes[] = {OptionType::Call, OptionType::Put};
} // namespace

TEMPLATE_TEST_CASE("Bachelier: put-call parity, no-arbitrage bounds, ATM closed form", "[vol][bachelier]", double, Real) {
    using T = TestType;
    for (double f : {0.035, -0.004, 0.0}) {
        for (double k : {0.02, 0.035, -0.01, 0.06}) {
            for (double tt : {0.25, 1.0, 10.0}) {
                for (double s : {0.0050, 0.0120}) {
                    const double c = value(bachelierPrice(T(f), k, tt, T(s), OptionType::Call));
                    const double p = value(bachelierPrice(T(f), k, tt, T(s), OptionType::Put));
                    CHECK_THAT(c - p, WithinAbs(f - k, 1e-15));
                    CHECK(c >= std::max(f - k, 0.0) - 1e-18);
                    CHECK(p >= std::max(k - f, 0.0) - 1e-18);
                }
            }
        }
    }
    // ATM: C = P = σ √T / √(2π).
    CHECK_THAT(value(bachelierPrice(T(0.03), 0.03, 2.0, T(0.009), OptionType::Call)), WithinRel(0.009 * std::sqrt(2.0) * kInvSqrt2Pi, 1e-15));
    CHECK_THAT(value(bachelierPrice(T(0.03), 0.03, 2.0, T(0.009), OptionType::Put)), WithinRel(0.009 * std::sqrt(2.0) * kInvSqrt2Pi, 1e-15));
    rl::aad::Tape::current().reset();
}

TEST_CASE("Bachelier: greeks match central finite differences", "[vol][bachelier][greeks]") {
    for (double k : {0.015, 0.030, 0.035, 0.050}) {
        for (double tt : {1.0 / 12.0, 1.0, 10.0}) {
            for (OptionType type : kTypes) {
                const double f = 0.035, s = 0.0085;
                const double hv = 1e-7, hf = 1e-7, hg = 1e-6;
                const double vegaFd = (bachelierPrice(f, k, tt, s + hv, type) - bachelierPrice(f, k, tt, s - hv, type)) / (2.0 * hv);
                const double deltaFd = (bachelierPrice(f + hf, k, tt, s, type) - bachelierPrice(f - hf, k, tt, s, type)) / (2.0 * hf);
                const double gammaFd = (bachelierDelta(f + hg, k, tt, s, type) - bachelierDelta(f - hg, k, tt, s, type)) / (2.0 * hg);
                CHECK(std::fabs(bachelierVega(f, k, tt, s) - vegaFd) <= 1e-7 * (1.0 + std::fabs(vegaFd)));
                CHECK(std::fabs(bachelierDelta(f, k, tt, s, type) - deltaFd) <= 1e-7 * (1.0 + std::fabs(deltaFd)));
                CHECK(std::fabs(bachelierGamma(f, k, tt, s) - gammaFd) <= 1e-6 * (1.0 + std::fabs(gammaFd)));
            }
            const double dc = bachelierDelta(0.035, k, tt, 0.0085, OptionType::Call);
            CHECK(dc >= 0.0);
            CHECK(dc <= 1.0);
            CHECK_THAT(bachelierDelta(0.035, k, tt, 0.0085, OptionType::Put), WithinAbs(dc - 1.0, 1e-15));
        }
    }
}

TEST_CASE("Bachelier: zero vol or zero expiry is intrinsic value", "[vol][bachelier]") {
    CHECK(bachelierPrice(0.04, 0.03, 1.0, 0.0, OptionType::Call) == 0.04 - 0.03);
    CHECK(bachelierPrice(0.04, 0.03, 1.0, 0.0, OptionType::Put) == 0.0);
    CHECK(bachelierPrice(0.02, 0.03, 0.0, 0.01, OptionType::Put) == 0.03 - 0.02);
    CHECK(bachelierDelta(0.04, 0.03, 1.0, 0.0, OptionType::Call) == 1.0);
    CHECK(bachelierDelta(0.04, 0.03, 1.0, 0.0, OptionType::Put) == 0.0);
    CHECK(bachelierVega(0.04, 0.03, 1.0, 0.0) == 0.0);
    CHECK_THAT(bachelierVega(0.03, 0.03, 4.0, 0.0), WithinRel(2.0 * kInvSqrt2Pi, 1e-15)); // the ATM limit is not zero
    CHECK(bachelierGamma(0.04, 0.03, 1.0, 0.0) == 0.0);
    CHECK_THROWS_AS(bachelierPrice(0.04, 0.03, -1.0, 0.01, OptionType::Call), rl::InputError);
    CHECK_THROWS_AS(bachelierPrice(0.04, 0.03, 1.0, -0.01, OptionType::Call), rl::InputError);
}

TEST_CASE("Bachelier implied vol round-trips to 1e-12 on OTM options out to six standard deviations", "[vol][bachelier][implied]") {
    const double f = 0.035;
    for (double s : {0.0030, 0.0085, 0.0150}) {
        for (double tt : {1.0 / 12.0, 1.0, 10.0}) {
            for (int x = -6; x <= 6; ++x) {
                const double k = f + x * s * std::sqrt(tt);
                const OptionType otm = k >= f ? OptionType::Call : OptionType::Put;
                const double price = bachelierPrice(f, k, tt, s, otm);
                INFO("sigma " << s << " T " << tt << " x " << x);
                CHECK_THAT(bachelierImpliedVol(price, f, k, tt, otm), WithinRel(s, 1e-12));
                // In-the-money prices carry intrinsic value; the time value loses digits in the subtraction, so the
                // inversion is only as good as the price allows (this is why the conversion spec says: invert OTM).
                if (std::abs(x) <= 2) {
                    const OptionType itm = otm == OptionType::Call ? OptionType::Put : OptionType::Call;
                    CHECK_THAT(bachelierImpliedVol(bachelierPrice(f, k, tt, s, itm), f, k, tt, itm), WithinRel(s, 1e-12));
                }
            }
        }
    }
    // Negative forward and strike: normal vols do not care.
    CHECK_THAT(bachelierImpliedVol(bachelierPrice(-0.004, -0.001, 2.0, 0.006, OptionType::Call), -0.004, -0.001, 2.0, OptionType::Call), WithinRel(0.006, 1e-12));
}

TEST_CASE("Bachelier implied vol: intrinsic, below intrinsic, bad expiry", "[vol][bachelier][implied]") {
    CHECK(bachelierImpliedVol(0.01, 0.04, 0.03, 1.0, OptionType::Call) == 0.0);
    CHECK(bachelierImpliedVol(0.0, 0.02, 0.03, 1.0, OptionType::Call) == 0.0);
    CHECK_THROWS_AS(bachelierImpliedVol(0.009, 0.04, 0.03, 1.0, OptionType::Call), rl::InputError);
    CHECK_THROWS_AS(bachelierImpliedVol(0.01, 0.04, 0.03, 0.0, OptionType::Call), rl::InputError);
}

TEST_CASE("Bachelier under Real: bit-identical to double and active", "[vol][bachelier][aad]") {
    const Real f = Real::input(0.035), s = Real::input(0.0085);
    for (OptionType type : kTypes) {
        const Real p = bachelierPrice(f, 0.03, 2.0, s, type);
        CHECK(p.value() == bachelierPrice(0.035, 0.03, 2.0, 0.0085, type));
        CHECK(p.isActive());
        CHECK(bachelierDelta(f, 0.03, 2.0, s, type).value() == bachelierDelta(0.035, 0.03, 2.0, 0.0085, type));
    }
    CHECK(bachelierVega(f, 0.03, 2.0, s).value() == bachelierVega(0.035, 0.03, 2.0, 0.0085));
    CHECK(bachelierGamma(f, 0.03, 2.0, s).value() == bachelierGamma(0.035, 0.03, 2.0, 0.0085));
    rl::aad::Tape::current().reset();
}
