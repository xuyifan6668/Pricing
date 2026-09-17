#include "rl/aad/tape.hpp"

#include "rl/core/error.hpp"

#include <catch2/catch_test_macros.hpp>

#include <thread>

using namespace rl::aad;

TEST_CASE("Tape records leaves and nodes with their parents", "[aad][tape]") {
    Tape& t = Tape::current();
    t.reset();
    REQUIRE(t.empty());
    auto a = t.newLeaf();
    auto b = t.newLeaf();
    auto c = t.record2(a, 2.0, b, 3.0);
    auto d = t.record1(c, 0.5);
    CHECK(t.size() == 4);
    CHECK(t.argumentCount() == 3);
    CHECK(t.arity(a) == 0);
    CHECK(t.arity(c) == 2);
    CHECK(t.parent(c, 0) == a);
    CHECK(t.partial(c, 1) == 3.0);
    CHECK(t.arity(d) == 1);
    CHECK(t.parent(d, 0) == c);
    t.reset();
}

TEST_CASE("Passive parents are dropped", "[aad][tape]") {
    Tape& t = Tape::current();
    t.reset();
    auto a = t.newLeaf();
    auto c = t.record2(Tape::kPassive, 1.0, a, 7.0);
    CHECK(t.arity(c) == 1);
    CHECK(t.parent(c, 0) == a);
    CHECK(t.partial(c, 0) == 7.0);
    auto d = t.record1(Tape::kPassive, 1.0);
    CHECK(t.arity(d) == 0);
    t.reset();
}

TEST_CASE("mark/rewind truncates", "[aad][tape]") {
    Tape& t = Tape::current();
    t.reset();
    auto a = t.newLeaf();
    auto b = t.newLeaf();
    auto pos = t.mark();
    t.record2(a, 1.0, b, 1.0);
    t.record1(a, 1.0);
    CHECK(t.size() == 4);
    t.rewindTo(pos);
    CHECK(t.size() == 2);
    CHECK(t.argumentCount() == 0);
    // Recording after rewind continues from the mark.
    auto c = t.record2(a, 1.0, b, 1.0);
    CHECK(c == 2);
    t.reset();
}

TEST_CASE("Adjoint slots exist per node; seeding and clearing", "[aad][tape]") {
    Tape& t = Tape::current();
    t.reset();
    auto a = t.newLeaf();
    auto b = t.record1(a, 2.0);
    CHECK(t.adjoint(b) == 0.0);
    t.seedAdjoint(b, 1.0);
    CHECK(t.adjoint(b) == 1.0);
    t.clearAdjoints();
    CHECK(t.adjoint(b) == 0.0);
    CHECK_THROWS_AS(t.seedAdjoint(99, 1.0), rl::InputError);
    t.reset();
}

TEST_CASE("propagate: hand-built tape y = 2a, z = 3y + 5a gives dz/da = 11", "[aad][tape]") {
    Tape& t = Tape::current();
    t.reset();
    auto a = t.newLeaf();
    auto y = t.record1(a, 2.0);
    auto z = t.record2(y, 3.0, a, 5.0);
    t.seedAdjoint(z, 1.0);
    t.propagate();
    CHECK(t.adjoint(y) == 3.0);
    CHECK(t.adjoint(a) == 11.0);
    // A range sweep that stops above the leaf leaves the leaf untouched.
    t.clearAdjoints();
    t.seedAdjoint(z, 1.0);
    t.propagate(z, z);
    CHECK(t.adjoint(y) == 3.0);
    CHECK(t.adjoint(a) == 5.0);
    t.reset();
}

TEST_CASE("Tapes are thread-local", "[aad][tape]") {
    Tape& t = Tape::current();
    t.reset();
    t.newLeaf();
    t.newLeaf();
    std::size_t otherSize = 999;
    Tape* otherPtr = nullptr;
    std::thread th([&] {
        Tape& o = Tape::current();
        otherPtr = &o;
        o.reset();
        otherSize = o.size();
        o.newLeaf();
    });
    th.join();
    CHECK(otherPtr != &t);
    CHECK(otherSize == 0);
    CHECK(t.size() == 2);
    t.reset();
}
