// MarketSnapshot<T> / SnapshotBuilder<T>: build, hash, immutability.
#include "rl/market/snapshot.hpp"
#include "rl/market/snapshot_builder.hpp"

#include "rl/aad/real.hpp"
#include "rl/core/error.hpp"
#include "rl/market/fixings.hpp"
#include "rl/market/quote.hpp"
#include "rl/market/swaption_cube.hpp"

#include <catch2/catch_test_macros.hpp>

#include <memory>
#include <type_traits>
#include <vector>

using namespace rl::market;
using rl::aad::Real;
using rl::core::Date;

namespace {
const Date kRef(2026, 9, 4);
template <class T>
CurveSet<T> curves(double bump = 0.0) {
    std::vector<double> ts{1.0, 5.0, 30.0};
    std::vector<T> ys{T(0.035), T(0.17 + bump), T(1.02)};
    return CurveSet<T>({YieldCurve<T>(kRef, "USD-SOFR-OIS", ts, ys, CurveInterpolation::LinearLogDf)}, "USD-SOFR-OIS", {{"USD-SOFR", "USD-SOFR-OIS"}});
}
FixingStore fixings() {
    FixingStore fs;
    fs.add("USD-SOFR", Date(2026, 9, 3), 0.0385);
    return fs;
}
QuoteSet quotes(double v = 3.5) {
    Quote q;
    q.id = "SOFR-OIS-1Y";
    q.type = QuoteType::OisSwap;
    q.index = "USD-SOFR";
    q.tenor = rl::core::Tenor::parse("1Y");
    q.value = v;
    return QuoteSet(kRef, rl::core::Currency::USD, {q});
}
} // namespace

TEST_CASE("Snapshot: build, read, hash determinism", "[market][snapshot]") {
    auto s1 = SnapshotBuilder<double>().asOf(kRef).curves(curves<double>()).fixings(fixings()).quotes(quotes()).build();
    auto s2 = SnapshotBuilder<double>().asOf(kRef).curves(curves<double>()).fixings(fixings()).quotes(quotes()).build();
    static_assert(std::is_same_v<decltype(s1), std::shared_ptr<const MarketSnapshot<double>>>);
    CHECK(s1->asOf() == kRef);
    CHECK(s1->curves().discount().name() == "USD-SOFR-OIS");
    CHECK(s1->fixings().fixing("USD-SOFR", Date(2026, 9, 3)) == 0.0385);
    CHECK(s1->quotes().find("SOFR-OIS-1Y").value == 3.5);
    CHECK(s1->hash() == s2->hash());
    CHECK(s1->hash() != 0);
}

TEST_CASE("Snapshot: hash is sensitive to every component", "[market][snapshot][hash]") {
    const auto base = SnapshotBuilder<double>().asOf(kRef).curves(curves<double>()).fixings(fixings()).quotes(quotes()).build()->hash();
    CHECK(SnapshotBuilder<double>().asOf(kRef).curves(curves<double>(1e-12)).fixings(fixings()).quotes(quotes()).build()->hash() != base);
    CHECK(SnapshotBuilder<double>().asOf(kRef).curves(curves<double>()).fixings(FixingStore{}).quotes(quotes()).build()->hash() != base);
    CHECK(SnapshotBuilder<double>().asOf(kRef).curves(curves<double>()).fixings(fixings()).quotes(quotes(3.51)).build()->hash() != base);
    CHECK(SnapshotBuilder<double>().asOf(kRef).curves(curves<double>()).fixings(fixings()).quotes(QuoteSet{}).build()->hash() != base);
}

TEST_CASE("Snapshot: hash is of the market, not of the tape (Real == double)", "[market][snapshot][aad]") {
    const auto hd = SnapshotBuilder<double>().asOf(kRef).curves(curves<double>()).fixings(fixings()).build()->hash();
    std::vector<double> ts{1.0, 5.0, 30.0};
    std::vector<Real> ys{Real::input(0.035), Real::input(0.17), Real::input(1.02)};
    CurveSet<Real> cr({YieldCurve<Real>(kRef, "USD-SOFR-OIS", ts, ys, CurveInterpolation::LinearLogDf)}, "USD-SOFR-OIS", {{"USD-SOFR", "USD-SOFR-OIS"}});
    const auto hr = SnapshotBuilder<Real>().asOf(kRef).curves(cr).fixings(fixings()).build()->hash();
    CHECK(hd == hr);
    rl::aad::Tape::current().reset();
}

TEST_CASE("Snapshot: immutable and value-safe", "[market][snapshot]") {
    static_assert(!std::is_assignable_v<decltype(std::declval<const MarketSnapshot<double>&>().curves()), CurveSet<double>>,
                  "curves() must return a const reference");
    static_assert(!std::is_assignable_v<decltype(std::declval<const MarketSnapshot<double>&>().fixings()), FixingStore>,
                  "fixings() must return a const reference");
    MarketSnapshot<double> copy = [] {
        auto s = SnapshotBuilder<double>().asOf(kRef).curves(curves<double>()).build();
        return *s; // the shared object dies here
    }();
    CHECK(copy.curves().discount().df(1.0) == curves<double>().discount().df(1.0));
    CHECK(copy.asOf() == kRef);
}

SwaptionCube<double> tinyCube(Date asOf = kRef, double nu = 0.4) {
    using rl::core::Tenor;
    std::vector<CubeNode<double>> ns(4);
    for (auto& n : ns) {
        n.forward = 0.035;
        n.atmNormalVol = 0.009;
        n.sabr = SabrParams<double>{0.05, -0.1, nu, 0.5, 0.0};
    }
    return SwaptionCube<double>(asOf, "USD-SOFR", {Tenor::years(1), Tenor::years(5)}, {Tenor::years(2), Tenor::years(10)}, ns, 0.5, 0.0);
}

TEST_CASE("Snapshot: optional swaption cube, hashed only when present", "[market][snapshot][cube]") {
    const auto without = SnapshotBuilder<double>().asOf(kRef).curves(curves<double>()).fixings(fixings()).quotes(quotes()).build();
    CHECK(without->swaptionCube() == nullptr);
    const auto with = SnapshotBuilder<double>().asOf(kRef).curves(curves<double>()).fixings(fixings()).quotes(quotes()).swaptionCube(tinyCube()).build();
    REQUIRE(with->swaptionCube() != nullptr);
    CHECK(with->swaptionCube()->size() == 4);
    CHECK(with->hash() != without->hash());
    const auto bumped = SnapshotBuilder<double>().asOf(kRef).curves(curves<double>()).fixings(fixings()).quotes(quotes()).swaptionCube(tinyCube(kRef, 0.41)).build();
    CHECK(bumped->hash() != with->hash());
    CHECK_THROWS_AS(SnapshotBuilder<double>().asOf(kRef).curves(curves<double>()).swaptionCube(tinyCube(Date(2026, 9, 3))).build(), rl::InputError);
}

TEST_CASE("SnapshotBuilder validation", "[market][snapshot]") {
    CHECK_THROWS_AS(SnapshotBuilder<double>().curves(curves<double>()).build(), rl::InputError);          // no asOf
    CHECK_THROWS_AS(SnapshotBuilder<double>().asOf(kRef).build(), rl::InputError);                         // no curves
    CHECK_THROWS_AS(SnapshotBuilder<double>().asOf(Date(2026, 9, 5)).curves(curves<double>()).build(), rl::InputError); // asOf ≠ ref
    CHECK_NOTHROW(SnapshotBuilder<double>().asOf(kRef).curves(curves<double>()).build());                 // fixings/quotes optional
}
