#pragma once
#include "rl/core/error.hpp"
#include "rl/market/hash.hpp"
#include "rl/market/snapshot.hpp"
#include "rl/market/snapshot_builder.hpp"

#include <utility>

namespace rl::market {

template <math::Scalar T>
MarketSnapshot<T>::MarketSnapshot(core::Date asOf, CurveSet<T> curves, FixingStore fixings, QuoteSet quotes, std::optional<SwaptionCube<T>> cube)
    : asOf_(asOf), curves_(std::move(curves)), fixings_(std::move(fixings)), quotes_(std::move(quotes)), cube_(std::move(cube)),
      hash_(computeHash(asOf_, curves_, fixings_, quotes_, cube_)) {}

template <math::Scalar T>
std::uint64_t MarketSnapshot<T>::computeHash(core::Date asOf, const CurveSet<T>& curves, const FixingStore& fixings, const QuoteSet& quotes,
                                             const std::optional<SwaptionCube<T>>& cube) noexcept {
    Hasher h;
    h.mix(asOf);
    h.mix(static_cast<std::uint64_t>(curves.size()));
    for (const auto& c : curves.curves()) {
        h.mix(c.name());
        h.mix(static_cast<int>(c.interpolation()));
        h.mix(c.nodeTimes());
        h.mix(static_cast<std::uint64_t>(c.nodeValues().size()));
        for (const T& y : c.nodeValues()) h.mix(math::value(y)); // value(): the hash is of the market, not the tape
        h.mix(static_cast<std::uint64_t>(c.turns().size()));
        for (const auto& tj : c.turns()) h.mix(tj.start).mix(tj.end).mix(math::value(tj.jump));
    }
    h.mix(curves.discountName());
    h.mix(static_cast<std::uint64_t>(curves.projections().size()));
    for (const auto& [index, curve] : curves.projections()) h.mix(index).mix(curve);
    fixings.mixInto(h);
    quotes.mixInto(h);
    // DESIGN: mixed only when present, so snapshot hashes (and anything recorded against them) are unchanged.
    if (cube) {
        h.mix(std::string_view("swaption_cube"));
        cube->mixInto(h);
    }
    return h.value();
}

template <math::Scalar T>
std::shared_ptr<const MarketSnapshot<T>> SnapshotBuilder<T>::build() const {
    RL_REQUIRE(asOf_.has_value(), "SnapshotBuilder: asOf is required");
    RL_REQUIRE(curves_.has_value(), "SnapshotBuilder: curves are required");
    RL_REQUIRE(curves_->referenceDate() == *asOf_, "SnapshotBuilder: asOf " + asOf_->toIso() + " differs from the curves' reference date " +
                                                        curves_->referenceDate().toIso());
    if (cube_) RL_REQUIRE(cube_->asOf() == *asOf_, "SnapshotBuilder: swaption cube as-of " + cube_->asOf().toIso() + " differs from " + asOf_->toIso());
    // Private constructor: make_shared cannot reach it, so construct then wrap.
    return std::shared_ptr<const MarketSnapshot<T>>(new MarketSnapshot<T>(*asOf_, *curves_, fixings_, quotes_, cube_));
}

} // namespace rl::market

#define RL_MARKET_INSTANTIATE_SNAPSHOT(T)                                                              \
    template class rl::market::MarketSnapshot<T>;                                                      \
    template class rl::market::SnapshotBuilder<T>;
