#pragma once
// MarketSnapshot<T>: one immutable object answering "what did the market look like when this number
// was computed" — curves, historical fixings, the quotes the curves came from, an as-of date and a
// content hash.
//
// Immutability is enforced by the type, not by convention: there is no non-const member function,
// the only constructor is private to SnapshotBuilder, and build() hands out shared_ptr<const>.
// Copies are independent values (no shared mutable state). Thread-safe to read from any number of
// risk threads.
//
// Active set: curve node values are T — they may be tape leaves. Fixings are
// passive double: a realised SOFR print is a number, not something a price is differentiated
// against. Quotes are stored as data; the bootstrap Jacobian, not the tape, carries sensitivities to
// them.
//
// The hash covers the *market*: as-of, every curve's name, kind, node times, value(node values) and
// turns, every fixing, every quote, and the swaption cube when one is present (an absent cube contributes
// nothing, so a curves-only snapshot hashes the same with or without the cube code). It is the same under T = double and T = Real for the same values,
// and identical across compilers (IEEE bit patterns, market/hash.hpp).

#include "rl/core/date.hpp"
#include "rl/market/curve_set.hpp"
#include "rl/market/fixings.hpp"
#include "rl/market/quote.hpp"
#include "rl/market/swaption_cube.hpp"

#include <cstdint>
#include <memory>
#include <optional>

namespace rl::market {

template <math::Scalar T>
class SnapshotBuilder;

template <math::Scalar T>
class MarketSnapshot {
public:
    core::Date asOf() const noexcept { return asOf_; }
    const CurveSet<T>& curves() const noexcept { return curves_; }
    const FixingStore& fixings() const noexcept { return fixings_; }
    const QuoteSet& quotes() const noexcept { return quotes_; }
    std::uint64_t hash() const noexcept { return hash_; }
    /// The calibrated swaption cube, or nullptr for a curves-only snapshot. DESIGN: a pointer into
    /// the snapshot rather than std::optional by reference, so callers write `if (auto* c = s.swaptionCube())`
    /// and cannot copy a cube by accident; it lives exactly as long as the snapshot.
    const SwaptionCube<T>* swaptionCube() const noexcept { return cube_ ? &*cube_ : nullptr; }

    // DESIGN: no with*() modifiers on the snapshot itself. A bumped market is a new build through the
    // builder (the bump driver does exactly that), so there is exactly one place that computes a
    // hash and exactly one place that validates as-of against the curves.

private:
    friend class SnapshotBuilder<T>;
    MarketSnapshot(core::Date asOf, CurveSet<T> curves, FixingStore fixings, QuoteSet quotes, std::optional<SwaptionCube<T>> cube);
    static std::uint64_t computeHash(core::Date asOf, const CurveSet<T>& curves, const FixingStore& fixings, const QuoteSet& quotes,
                                     const std::optional<SwaptionCube<T>>& cube) noexcept;

    core::Date asOf_;
    CurveSet<T> curves_;
    FixingStore fixings_;
    QuoteSet quotes_;
    std::optional<SwaptionCube<T>> cube_;
    std::uint64_t hash_ = 0;
};

} // namespace rl::market
