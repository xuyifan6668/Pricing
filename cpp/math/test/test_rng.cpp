#include "rl/math/rng.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <set>

using rl::math::Rng;
using rl::math::SeedPolicy;
using Catch::Matchers::WithinAbs;

TEST_CASE("mt19937_64 stream is the standard one", "[math][rng]") {
    // [rand.predef]: the 10000th value of a default-seeded (5489) mt19937_64 is 9981545732273789042.
    Rng r(5489);
    r.discard(9999);
    CHECK(r.nextU64() == 9981545732273789042ull);
}

TEST_CASE("Uniforms live strictly inside (0, 1) and use the top 53 bits", "[math][rng]") {
    Rng r(1);
    double mn = 1.0, mx = 0.0;
    for (int i = 0; i < 100000; ++i) {
        double u = r.uniform01();
        mn = std::min(mn, u);
        mx = std::max(mx, u);
    }
    CHECK(mn > 0.0);
    CHECK(mx < 1.0);
    CHECK(mn < 1e-3);
    CHECK(mx > 0.999);
    // The offset makes every value an odd multiple of 2^-54.
    Rng r2(1);
    double u = r2.uniform01();
    double k = u * 0x1.0p54;
    CHECK(k == std::floor(k));
    CHECK(std::fmod(k, 2.0) == 1.0);
}

TEST_CASE("Seed policy is deterministic, order-sensitive, and spreads", "[math][rng]") {
    CHECK(SeedPolicy::derive(1, 2) == SeedPolicy::derive(1, 2));
    CHECK(SeedPolicy::derive(1, 2) != SeedPolicy::derive(2, 1));
    CHECK(SeedPolicy::derive(1, 2, 0) != SeedPolicy::derive(1, 2, 1));
    CHECK(SeedPolicy::derive(0, 0) != 0);
    std::set<std::uint64_t> seen;
    for (std::uint64_t t = 0; t < 100; ++t)
        for (std::uint64_t s = 0; s < 100; ++s) seen.insert(SeedPolicy::derive(t, s));
    CHECK(seen.size() == 10000);
    // splitmix64 reference: mix(0) is the well-known first output of the zero-seeded generator.
    CHECK(SeedPolicy::mix(0) == 0xE220A8397B1DCDAFull);
}

TEST_CASE("Gaussian draws have the right moments and are reproducible", "[math][rng]") {
    Rng a(SeedPolicy::derive(7, 3)), b(SeedPolicy::derive(7, 3));
    const int n = 200000;
    double s1 = 0, s2 = 0, s3 = 0, s4 = 0;
    for (int i = 0; i < n; ++i) {
        double z = a.gaussian();
        CHECK(z == b.gaussian());
        s1 += z;
        s2 += z * z;
        s3 += z * z * z;
        s4 += z * z * z * z;
    }
    CHECK_THAT(s1 / n, WithinAbs(0.0, 0.01));
    CHECK_THAT(s2 / n, WithinAbs(1.0, 0.02));
    CHECK_THAT(s3 / n, WithinAbs(0.0, 0.03));
    CHECK_THAT(s4 / n, WithinAbs(3.0, 0.1));
}
