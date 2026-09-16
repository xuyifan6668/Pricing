#include "rl/math/sobol.hpp"

#include "rl/core/error.hpp"
#include "rl/math/rng.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <vector>

using rl::math::DirectionNumbers;
using rl::math::Sobol;
using Catch::Matchers::WithinAbs;

TEST_CASE("Joe-Kuo direction numbers load", "[math][sobol]") {
    const auto& dn = DirectionNumbers::joeKuo();
    CHECK(dn.maxDimension() >= 1111);
    REQUIRE(dn.dims.size() >= 3);
    // Dimension 2: s=1, a=0, m={1}. Dimension 3: s=2, a=1, m={1,3}. Dimension 4: s=3, a=1, m={1,3,1}.
    CHECK(dn.dims[0].d == 2);
    CHECK(dn.dims[0].s == 1);
    CHECK(dn.dims[0].m == std::vector<unsigned>{1});
    CHECK(dn.dims[1].s == 2);
    CHECK(dn.dims[1].a == 1);
    CHECK(dn.dims[1].m == std::vector<unsigned>{1, 3});
    CHECK(dn.dims[2].s == 3);
    CHECK(dn.dims[2].m == std::vector<unsigned>{1, 3, 1});
    CHECK_THROWS_AS(DirectionNumbers::load("does/not/exist.txt"), rl::InputError);
}

TEST_CASE("First eight points match the Joe-Kuo reference (dimensions 1-3)", "[math][sobol]") {
    // Output of Joe & Kuo's sobol.cc for N = 8, D = 3 (Gray-code order, origin first).
    const double ref[8][3] = {
        {0.0, 0.0, 0.0},       {0.5, 0.5, 0.5},      {0.75, 0.25, 0.25},   {0.25, 0.75, 0.75},
        {0.375, 0.375, 0.625}, {0.875, 0.875, 0.125}, {0.625, 0.125, 0.875}, {0.125, 0.625, 0.375},
    };
    Sobol s(3);
    for (int n = 0; n < 8; ++n) {
        auto p = s.next();
        REQUIRE(p.size() == 3);
        for (int j = 0; j < 3; ++j) CHECK(p[static_cast<std::size_t>(j)] == ref[n][j]);
    }
    CHECK(s.index() == 8);
}

TEST_CASE("Dimension 1 is the van der Corput sequence in base 2", "[math][sobol]") {
    Sobol s(1);
    std::vector<double> got;
    for (int n = 0; n < 8; ++n) got.push_back(s.next()[0]);
    CHECK(got == std::vector<double>{0.0, 0.5, 0.75, 0.25, 0.375, 0.875, 0.625, 0.125});
}

TEST_CASE("Skip and jump are consistent with sequential generation", "[math][sobol]") {
    for (std::uint64_t k : {0ull, 1ull, 7ull, 8ull, 1000ull, 65535ull}) {
        Sobol fresh(5);
        for (std::uint64_t i = 0; i < k; ++i) fresh.next();
        Sobol jumped(5, k);
        CHECK(jumped.index() == k);
        auto a = fresh.next(), b = jumped.next();
        CHECK(a == b);
        Sobol skipLater(5);
        skipLater.skip(k);
        CHECK(skipLater.next() == a);
    }
}

TEST_CASE("Stratification: the first 2^m points hit every dyadic interval exactly once per coordinate", "[math][sobol]") {
    const unsigned m = 10, dim = 8;
    const std::size_t N = std::size_t(1) << m;
    Sobol s(dim);
    std::vector<std::vector<int>> counts(dim, std::vector<int>(N, 0));
    std::vector<double> p(dim);
    for (std::size_t n = 0; n < N; ++n) {
        s.next(p);
        for (unsigned j = 0; j < dim; ++j) {
            REQUIRE(p[j] >= 0.0);
            REQUIRE(p[j] < 1.0);
            ++counts[j][static_cast<std::size_t>(p[j] * static_cast<double>(N))];
        }
    }
    for (unsigned j = 0; j < dim; ++j)
        for (std::size_t b = 0; b < N; ++b) CHECK(counts[j][b] == 1);
}

TEST_CASE("Two-dimensional (0,m,2)-net property for dimensions 1 and 2", "[math][sobol]") {
    const unsigned m = 8;
    const std::size_t N = std::size_t(1) << m;
    for (unsigned a = 0; a <= m; ++a) {
        const unsigned b = m - a;
        const std::size_t na = std::size_t(1) << a, nb = std::size_t(1) << b;
        std::vector<int> counts(na * nb, 0);
        Sobol s(2);
        std::vector<double> p(2);
        for (std::size_t n = 0; n < N; ++n) {
            s.next(p);
            auto ia = static_cast<std::size_t>(p[0] * static_cast<double>(na));
            auto ib = static_cast<std::size_t>(p[1] * static_cast<double>(nb));
            ++counts[ia * nb + ib];
        }
        for (int c : counts) CHECK(c == 1);
    }
}

TEST_CASE("Sobol beats pseudo-random on exp(-|x|_1) over [0,1]^4", "[math][sobol]") {
    const double exact = std::pow(1.0 - std::exp(-1.0), 4);
    auto integrand = [](const std::vector<double>& x) { return std::exp(-(x[0] + x[1] + x[2] + x[3])); };
    const std::size_t N = std::size_t(1) << 16;
    Sobol s(4, 1); // skip the origin
    std::vector<double> p(4);
    double acc = 0.0;
    for (std::size_t n = 0; n < N; ++n) {
        s.next(p);
        acc += integrand(p);
    }
    const double sobolErr = std::fabs(acc / static_cast<double>(N) - exact);

    double mtErr = 0.0;
    const int seeds = 8;
    for (int k = 0; k < seeds; ++k) {
        rl::math::Rng rng(rl::math::SeedPolicy::derive(99, static_cast<std::uint64_t>(k)));
        double a2 = 0.0;
        for (std::size_t n = 0; n < N; ++n) {
            for (auto& v : p) v = rng.uniform01();
            a2 += integrand(p);
        }
        mtErr += std::fabs(a2 / static_cast<double>(N) - exact);
    }
    mtErr /= seeds;
    INFO("sobol error " << sobolErr << " vs mean MT error " << mtErr);
    CHECK(sobolErr < 0.25 * mtErr);
    CHECK(sobolErr < 1e-4);
}

TEST_CASE("Validation", "[math][sobol]") {
    CHECK_THROWS_AS(Sobol(0), rl::InputError);
    CHECK_THROWS_AS(Sobol(100000), rl::InputError);
    Sobol s(2);
    std::vector<double> tooSmall(1);
    CHECK_THROWS_AS(s.next(tooSmall), rl::InputError);
}
