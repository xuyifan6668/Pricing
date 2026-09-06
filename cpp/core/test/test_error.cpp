#include "rl/core/error.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_exception.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <string>

using namespace rl;
using Catch::Matchers::ContainsSubstring;

TEST_CASE("Error carries a message and an empty context", "[core][error]") {
    Error e("boom");
    CHECK(e.message() == "boom");
    CHECK(std::string(e.what()) == "boom");
    CHECK(e.context().empty());
}

TEST_CASE("withContext appends context innermost-first and prints outermost-first", "[core][error]") {
    try {
        withContext("bootstrap", [] {
            return withContext("USD-SOFR", [] { return withContext("node 2y", []() -> int { throw NumericalError("no convergence"); }); });
        });
        FAIL("should have thrown");
    } catch (const NumericalError& e) {
        REQUIRE(e.context().size() == 3);
        CHECK(e.context()[0] == "node 2y");
        CHECK(e.context()[2] == "bootstrap");
        CHECK(std::string(e.what()) == "no convergence [context: bootstrap > USD-SOFR > node 2y]");
    }
}

TEST_CASE("withContext passes non-rl exceptions through untouched", "[core][error]") {
    CHECK_THROWS_AS(withContext("x", [] { throw std::logic_error("std"); }), std::logic_error);
}

TEST_CASE("withContext returns the callable's value", "[core][error]") {
    CHECK(withContext("x", [] { return 42; }) == 42);
}

TEST_CASE("RL_REQUIRE throws InputError with the expression and location", "[core][error]") {
    auto f = [](int x) { RL_REQUIRE(x > 0, "x must be positive"); return x; };
    CHECK(f(1) == 1);
    CHECK_THROWS_MATCHES(f(-1), InputError, Catch::Matchers::MessageMatches(ContainsSubstring("x must be positive") && ContainsSubstring("x > 0") && ContainsSubstring("test_error.cpp")));
}

TEST_CASE("NotImplemented names the function", "[core][error]") {
    auto f = []() -> double { RL_NOT_IMPLEMENTED("Foo::bar"); };
    CHECK_THROWS_MATCHES(f(), NotImplemented, Catch::Matchers::MessageMatches(ContainsSubstring("not implemented: Foo::bar")));
}

TEST_CASE("CalibrationError keeps residuals", "[core][error]") {
    CalibrationError e("fit failed", {1e-3, -2e-3});
    CHECK(e.residuals().size() == 2);
    CHECK(e.residuals()[1] == -2e-3);
}

TEST_CASE("Hierarchy: everything is an rl::Error and a std::exception", "[core][error]") {
    auto thrower = [](int kind) {
        if (kind == 0) throw InputError("a");
        if (kind == 1) throw NumericalError("a");
        if (kind == 2) throw NotImplemented("a");
        throw CalibrationError("a", {});
    };
    CHECK_THROWS_AS(thrower(0), Error);
    CHECK_THROWS_AS(thrower(1), Error);
    CHECK_THROWS_AS(thrower(2), Error);
    CHECK_THROWS_AS(thrower(3), std::exception);
}
