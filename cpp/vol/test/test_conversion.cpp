// Normal <-> shifted-lognormal conversion: round trips, price agreement, shift handling, domain errors.
#include "rl/vol/conversion.hpp"

#include "rl/core/error.hpp"
#include "rl/vol/bachelier.hpp"
#include "rl/vol/black.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <algorithm>
#include <cmath>

using namespace rl::vol;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {
OptionType otmType(double f, double k) { return k >= f ? OptionType::Call : OptionType::Put; }
} // namespace

TEST_CASE("ATM series: sigma_N = sigma_LN (F + s)(1 - sigma^2 T/24) to O(sigma^4 T^2)", "[vol][conversion][atm]") {
    // Exact ATM conversion from pieces: price under Black, invert under Bachelier.
    for (double shift : {0.0, 0.02}) {
        for (auto [s, tt] : {std::pair{0.2, 1.0}, std::pair{0.3, 10.0}, std::pair{0.5, 5.0}}) {
            const double f = 0.03;
            const double exact = bachelierImpliedVol(blackPrice(f, f, tt, s, shift, OptionType::Call), f, f, tt, OptionType::Call);
            const double approx = atmNormalFromLognormalApprox(f, tt, s, shift);
            const double u2 = s * s * tt;
            // The next term of the series is +u⁴/640: the approximation is low by that much, relatively.
            CHECK_THAT(1.0 - approx / exact, WithinRel(u2 * u2 / 640.0, 0.1));
            // The inverse series undoes the forward series to the same order.
            CHECK_THAT(atmLognormalFromNormalApprox(f, tt, approx, shift), WithinRel(s, 2.0 * u2 * u2 / 24.0));
        }
    }
    CHECK_THROWS_AS(atmNormalFromLognormalApprox(-0.03, 1.0, 0.2, 0.01), rl::InputError);
    CHECK_THROWS_AS(atmLognormalFromNormalApprox(-0.03, 1.0, 0.01, 0.01), rl::InputError);
}

TEST_CASE("Normal -> shifted lognormal -> normal round-trips to 1e-12 and the prices agree", "[vol][conversion]") {
    for (double shift : {0.0, 0.01, 0.03}) {
        for (double fwd : {0.035, 0.002, -0.004}) {
            if (fwd + shift <= 0.0) continue;
            for (double tt : {0.25, 2.0, 10.0}) {
                for (double sn : {0.005, 0.009, 0.015}) {
                    for (double off : {-0.01, -0.0025, 0.0, 0.005, 0.01}) {
                        const double k = fwd + off;
                        if (k + shift <= 1e-4) continue;
                        const OptionType t = otmType(fwd, k);
                        const double pn = bachelierPrice(fwd, k, tt, sn, t);
                        if (pn >= 0.99 * std::min(fwd + shift, k + shift)) continue; // not attainable in lognormal: separate case
                        INFO("shift " << shift << " F " << fwd << " T " << tt << " sigmaN " << sn << " K " << k);
                        const double sl = normalToShiftedLognormalVol(fwd, k, tt, sn, shift);
                        CHECK(sl > 0.0);
                        CHECK_THAT(blackPrice(fwd, k, tt, sl, shift, t), WithinRel(pn, 1e-12));
                        CHECK_THAT(shiftedLognormalToNormalVol(fwd, k, tt, sl, shift), WithinRel(sn, 1e-12));
                    }
                }
            }
        }
    }
}

TEST_CASE("Shifted lognormal -> normal -> shifted lognormal round-trips to 1e-12", "[vol][conversion]") {
    for (double shift : {0.0, 0.02}) {
        for (double tt : {0.5, 5.0}) {
            for (double sl : {0.10, 0.30, 0.80}) {
                for (double off : {-0.01, 0.0, 0.01}) {
                    const double fwd = 0.03, k = fwd + off;
                    const double sn = shiftedLognormalToNormalVol(fwd, k, tt, sl, shift);
                    CHECK(sn > 0.0);
                    CHECK_THAT(normalToShiftedLognormalVol(fwd, k, tt, sn, shift), WithinRel(sl, 1e-12));
                    CHECK_THAT(bachelierPrice(fwd, k, tt, sn, otmType(fwd, k)), WithinRel(blackPrice(fwd, k, tt, sl, shift, otmType(fwd, k)), 1e-12));
                }
            }
        }
    }
}

TEST_CASE("Conversion: shift conventions, monotonicity, and domain errors", "[vol][conversion]") {
    // A larger shift needs a smaller lognormal vol for the same normal vol: σ_LN ≈ σ_N/(F + s).
    const double s0 = normalToShiftedLognormalVol(0.02, 0.025, 2.0, 0.009, 0.0);
    const double s3 = normalToShiftedLognormalVol(0.02, 0.025, 2.0, 0.009, 0.03);
    CHECK(s3 < s0);
    CHECK_THAT(s0, WithinRel(0.009 / 0.0225, 0.05));
    // The normal vol is a property of the price, not of the quoting shift.
    CHECK_THAT(shiftedLognormalToNormalVol(0.02, 0.025, 2.0, s3, 0.03), WithinRel(0.009, 1e-12));
    // Monotone in the input vol.
    CHECK(shiftedLognormalToNormalVol(0.03, 0.035, 1.0, 0.31, 0.0) > shiftedLognormalToNormalVol(0.03, 0.035, 1.0, 0.30, 0.0));
    // Negative forward with a shift: the case this function exists for.
    CHECK_THAT(shiftedLognormalToNormalVol(-0.004, -0.002, 5.0, normalToShiftedLognormalVol(-0.004, -0.002, 5.0, 0.006, 0.02), 0.02), WithinRel(0.006, 1e-12));
    // Zero vol maps to zero vol.
    CHECK(normalToShiftedLognormalVol(0.03, 0.035, 1.0, 0.0, 0.0) == 0.0);
    // Domain.
    CHECK_THROWS_AS(normalToShiftedLognormalVol(-0.02, 0.01, 1.0, 0.01, 0.01), rl::InputError);  // F + s ≤ 0
    CHECK_THROWS_AS(shiftedLognormalToNormalVol(-0.02, 0.01, 1.0, 0.30, 0.01), rl::InputError);  // F + s ≤ 0
    CHECK_THROWS_AS(normalToShiftedLognormalVol(0.01, -0.02, 1.0, 0.01, 0.01), rl::InputError);  // K + s ≤ 0
    CHECK_THROWS_AS(normalToShiftedLognormalVol(0.03, 0.03, 0.0, 0.01, 0.0), rl::InputError);    // T = 0
    // Unattainable: an OTM put whose normal price exceeds K + s cannot be a lognormal price.
    CHECK_THROWS_AS(normalToShiftedLognormalVol(0.012, 0.0005, 10.0, 0.03, 0.0), rl::InputError);
}
