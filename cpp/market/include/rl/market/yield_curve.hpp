#pragma once
// YieldCurve<T>: a dumb container — pillars, ordinates, an interpolation rule, a time convention, and
// an optional turn-of-year jump layer. It knows nothing about swaps, quotes or "par"; that is linear/
// and bootstrap/ (docs/DESIGN.md §1, "market | linear"; docs/DESIGN.md "Curve object model").
//
// Value type, immutable after construction. Copy freely; share as shared_ptr<const> when a container
// (snapshot, Excel handle) owns it. Thread-safe to read.
//
// Parametrisation: the ordinates are y(t) = −log P(0,t) at the pillar times (ACT/365F years from the
// reference date), interpolated by the scheme named in `interpolation`. t = 0 with y = 0 is an implicit
// anchor, not a free node: the bootstrap has one parameter per instrument and t = 0 is not an instrument.
//
// Turn-of-year (jump layer): df(t) = df_smooth(t) · exp(−Σ_k J_k · overlap(t, [a_k, b_k])), with
// overlap the length of [0, t] ∩ [a_k, b_k] in years. The overnight forward is J_k higher over exactly
// the turn window and nothing else moves: the smooth interpolator never sees the turn, and the
// bootstrap sees one extra parameter per turn.

#include "rl/core/date.hpp"
#include "rl/market/curve_interpolation.hpp"
#include "rl/math/interpolation.hpp"
#include "rl/math/scalar.hpp"

#include <string>
#include <vector>

namespace rl::market {

/// A forward-rate jump of `jump` (rate units, active) over the curve-time window [start, end).
template <math::Scalar T>
struct TurnJump {
    double start = 0.0;
    double end = 0.0;
    T jump = T(0.0);
};

template <math::Scalar T>
class YieldCurve {
public:
    /// Throws InputError on: empty nodes, non-increasing or non-positive times, size mismatch,
    /// non-finite values, or a turn window that is empty or starts before the reference date.
    YieldCurve(core::Date referenceDate, std::string name, std::vector<double> nodeTimes, std::vector<T> nodeValues,
               CurveInterpolation interpolation, std::vector<TurnJump<T>> turns = {});

    core::Date referenceDate() const noexcept { return ref_; }
    const std::string& name() const noexcept { return name_; }
    CurveInterpolation interpolation() const noexcept { return kind_; }

    /// ACT/365F year fraction from the reference date; negative before it.
    double time(core::Date d) const;

    /// P(0,t). Exactly T(1) at t == 0; flat-forward extrapolation beyond the last node and before the
    /// first (the implicit anchor makes that the first interval's forward). Never throws for t ≥ 0.
    T df(double t) const;
    T df(core::Date d) const { return df(time(d)); }
    /// The smooth part only (no turns). Equals df() when there are no turns.
    T dfSmooth(double t) const;
    /// y(t) = −log P(0,t) including turns.
    T logDf(double t) const;

    const std::vector<double>& nodeTimes() const noexcept { return ts_; }   ///< free nodes only, t > 0, strictly increasing
    const std::vector<T>& nodeValues() const noexcept { return ys_; }       ///< y = −log P at the nodes
    const std::vector<TurnJump<T>>& turns() const noexcept { return turns_; }

    /// Same reference date, name, times, kind and turns; new ordinates. Re-validates.
    YieldCurve withNodeValues(std::vector<T> values) const;
    /// Same everything; new turn jumps (count must match).
    YieldCurve withTurnJumps(const std::vector<T>& jumps) const;

private:
    void validateAndBuild();
    T turnAdjustment(double t) const; ///< Σ_k J_k · overlap(t, [a_k, b_k])

    core::Date ref_;
    std::string name_;
    std::vector<double> ts_;
    std::vector<T> ys_;
    CurveInterpolation kind_;
    std::vector<TurnJump<T>> turns_;
    // DESIGN: the interpolator variant is a member rebuilt in the constructor, not a lazily built
    // cache: an immutable object must not carry a `mutable` (it would need a mutex under the parallel
    // risk driver). Rebuilding is O(n) per withNodeValues, ~50 times per bootstrap — negligible.
    math::Interpolator<T> interp_;
};

} // namespace rl::market
