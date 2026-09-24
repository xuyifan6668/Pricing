#include "rl/market/hash.hpp"

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <limits>
#include <string>

using namespace rl::market;

TEST_CASE("FNV-1a 64 reference values", "[market][hash]") {
    // Published test vectors (Landon Curt Noll's FNV page): fnv1a64("") = offset basis,
    // fnv1a64("a") = 0xaf63dc4c8601ec8c, fnv1a64("foobar") = 0x85944171f73967e8.
    CHECK(fnv1a64("", 0) == Hasher::kOffsetBasis);
    CHECK(fnv1a64("a", 1) == 0xaf63dc4c8601ec8cull);
    CHECK(fnv1a64("foobar", 6) == 0x85944171f73967e8ull);
    CHECK(toHex(0xaf63dc4c8601ec8cull) == "0xaf63dc4c8601ec8c");
}

TEST_CASE("Hasher is deterministic and order/length sensitive", "[market][hash]") {
    Hasher a, b;
    a.mix("USD-SOFR").mix(1.5).mix(rl::core::Date(2026, 9, 4));
    b.mix("USD-SOFR").mix(1.5).mix(rl::core::Date(2026, 9, 4));
    CHECK(a.value() == b.value());
    Hasher c;
    c.mix(1.5).mix("USD-SOFR").mix(rl::core::Date(2026, 9, 4));
    CHECK(c.value() != a.value());
    // Length prefix: ("ab","c") ≠ ("a","bc").
    Hasher d, e;
    d.mix("ab").mix("c");
    e.mix("a").mix("bc");
    CHECK(d.value() != e.value());
}

TEST_CASE("Hasher: one ulp changes the hash; -0 and +0 do not differ", "[market][hash]") {
    const double x = 0.0385;
    const double y = std::nextafter(x, 1.0);
    CHECK(Hasher().mix(x).value() != Hasher().mix(y).value());
    CHECK(Hasher().mix(0.0).value() == Hasher().mix(-0.0).value());
    // Bumping by 1e-12 on a rate is many ulps: must differ.
    CHECK(Hasher().mix(x).value() != Hasher().mix(x + 1e-12).value());
    // NaN hashes to something (no throw) and differs from any number.
    CHECK(Hasher().mix(std::numeric_limits<double>::quiet_NaN()).value() != Hasher().mix(x).value());
}

TEST_CASE("Hasher: vector mix includes the length", "[market][hash]") {
    std::vector<double> v1{1.0, 2.0}, v2{1.0, 2.0, 0.0};
    CHECK(Hasher().mix(v1).value() != Hasher().mix(v2).value());
    CHECK(Hasher().mix(v1).value() == Hasher().mix(std::vector<double>{1.0, 2.0}).value());
}
