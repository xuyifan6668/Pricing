#pragma once
// SnapshotBuilder<T>: the only way to make a MarketSnapshot<T>. Fluent setters, then build(), which
// validates (as-of present and equal to the curves' reference date and to the cube's as-of; curves present),
// computes the hash once, and returns shared_ptr<const MarketSnapshot<T>>.
//
// The builder is reusable: build() copies its parts, so a risk driver can set the base once and call
// curves(bumped).build() per scenario. Setting a part twice replaces it (last wins).
//
// Where T = Real leaves come from: the caller builds a CurveSet<Real> whose node values are
// Real::input(...) in a documented order — curve by curve in CurveSet::names() order, node by node,
// then turn jumps — and hands it to curves(). DESIGN: the builder does not call Real::input itself,
// so the snapshot module stays independent of the tape and the adjoint driver owns the leaf order.

#include "rl/core/error.hpp"
#include "rl/market/snapshot.hpp"

#include <memory>
#include <optional>

namespace rl::market {

template <math::Scalar T>
class SnapshotBuilder {
public:
    SnapshotBuilder& asOf(core::Date d) {
        asOf_ = d;
        return *this;
    }
    SnapshotBuilder& curves(CurveSet<T> cs) {
        curves_ = std::move(cs);
        return *this;
    }
    SnapshotBuilder& fixings(FixingStore f) {
        fixings_ = std::move(f);
        return *this;
    }
    SnapshotBuilder& quotes(QuoteSet q) {
        quotes_ = std::move(q);
        return *this;
    }
    /// Optional. A vega-by-node scenario sets a bumped cube here and rebuilds, exactly like a curve bump.
    SnapshotBuilder& swaptionCube(SwaptionCube<T> c) {
        cube_ = std::move(c);
        return *this;
    }

    /// Throws InputError if asOf or curves are missing, asOf != curves.referenceDate(), or a cube's as-of differs.
    std::shared_ptr<const MarketSnapshot<T>> build() const;

private:
    std::optional<core::Date> asOf_;
    std::optional<CurveSet<T>> curves_;
    FixingStore fixings_;
    QuoteSet quotes_;
    std::optional<SwaptionCube<T>> cube_;
};

} // namespace rl::market
