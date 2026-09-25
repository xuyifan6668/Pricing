#pragma once
#include "rl/core/daycount.hpp"
#include "rl/core/error.hpp"
#include "rl/market/curve_set.hpp"
#include "rl/market/yield_curve.hpp"

#include <algorithm>
#include <cmath>
#include <utility>

namespace rl::market {

using std::exp;
using std::log;

// ------------------------------------------------------------------------------------------------
// YieldCurve
// ------------------------------------------------------------------------------------------------

template <math::Scalar T>
YieldCurve<T>::YieldCurve(core::Date referenceDate, std::string name, std::vector<double> nodeTimes, std::vector<T> nodeValues,
                          CurveInterpolation interpolation, std::vector<TurnJump<T>> turns)
    : ref_(referenceDate), name_(std::move(name)), ts_(std::move(nodeTimes)), ys_(std::move(nodeValues)), kind_(interpolation),
      turns_(std::move(turns)), interp_(math::Linear<T>({0.0, 1.0}, {T(0.0), T(0.0)})) {
    validateAndBuild();
}

template <math::Scalar T>
void YieldCurve<T>::validateAndBuild() {
    RL_REQUIRE(!ts_.empty(), "YieldCurve " + name_ + ": no nodes");
    RL_REQUIRE(ts_.size() == ys_.size(), "YieldCurve " + name_ + ": node times and values differ in size");
    RL_REQUIRE(ts_.front() > 0.0 && std::isfinite(ts_.front()), "YieldCurve " + name_ + ": node times must be positive");
    for (std::size_t i = 1; i < ts_.size(); ++i) RL_REQUIRE(ts_[i] > ts_[i - 1] && std::isfinite(ts_[i]), "YieldCurve " + name_ + ": node times not strictly increasing");
    for (const T& y : ys_) RL_REQUIRE(math::isFinite(y), "YieldCurve " + name_ + ": non-finite node value");
    for (const auto& tj : turns_) {
        RL_REQUIRE(tj.start >= 0.0 && tj.end > tj.start, "YieldCurve " + name_ + ": turn window must be non-empty and after the reference date");
        RL_REQUIRE(math::isFinite(tj.jump), "YieldCurve " + name_ + ": non-finite turn jump");
    }
    // The implicit (0, 0) anchor is prepended for the interpolator only; nodeTimes() never shows it.
    std::vector<double> ts;
    ts.reserve(ts_.size() + 1);
    ts.push_back(0.0);
    ts.insert(ts.end(), ts_.begin(), ts_.end());
    std::vector<T> ys;
    ys.reserve(ys_.size() + 1);
    ys.push_back(T(0.0));
    ys.insert(ys.end(), ys_.begin(), ys_.end());
    switch (kind_) {
    case CurveInterpolation::LinearLogDf:
        // DESIGN: Extrapolation::Linear on y continues the last segment's slope = flat instantaneous
        // forward beyond the last pillar, which is what a 30y curve must do when a 30y swap's payment
        // delay asks for P(0, 30y + 2d).
        interp_ = math::Linear<T>(std::move(ts), std::move(ys), math::Extrapolation::Linear);
        break;
    case CurveInterpolation::MonotoneConvex:
        interp_ = math::MonotoneConvex<T>(std::move(ts), std::move(ys));
        break;
    }
}

template <math::Scalar T>
double YieldCurve<T>::time(core::Date d) const {
    return core::yearFraction(core::DayCount::Act365Fixed, ref_, d);
}

template <math::Scalar T>
T YieldCurve<T>::turnAdjustment(double t) const {
    T adj(0.0);
    for (const auto& tj : turns_) {
        const double overlap = std::max(0.0, std::min(t, tj.end) - tj.start);
        if (overlap > 0.0) adj += tj.jump * overlap;
    }
    return adj;
}

template <math::Scalar T>
T YieldCurve<T>::dfSmooth(double t) const {
    if (t == 0.0) return T(1.0);
    return exp(-math::interpolate(interp_, t));
}

template <math::Scalar T>
T YieldCurve<T>::logDf(double t) const {
    if (t == 0.0) return T(0.0);
    T y = math::interpolate(interp_, t);
    if (!turns_.empty()) y += turnAdjustment(t);
    return y;
}

template <math::Scalar T>
T YieldCurve<T>::df(double t) const {
    if (t == 0.0) return T(1.0);
    return exp(-logDf(t));
}

template <math::Scalar T>
YieldCurve<T> YieldCurve<T>::withNodeValues(std::vector<T> values) const {
    return YieldCurve(ref_, name_, ts_, std::move(values), kind_, turns_);
}

template <math::Scalar T>
YieldCurve<T> YieldCurve<T>::withTurnJumps(const std::vector<T>& jumps) const {
    RL_REQUIRE(jumps.size() == turns_.size(), "YieldCurve " + name_ + ": turn jump count mismatch");
    std::vector<TurnJump<T>> turns = turns_;
    for (std::size_t k = 0; k < turns.size(); ++k) turns[k].jump = jumps[k];
    return YieldCurve(ref_, name_, ts_, ys_, kind_, std::move(turns));
}

// ------------------------------------------------------------------------------------------------
// CurveSet
// ------------------------------------------------------------------------------------------------

template <math::Scalar T>
CurveSet<T>::CurveSet(std::vector<YieldCurve<T>> curves, std::string discountCurveName, std::vector<std::pair<std::string, std::string>> projections)
    : curves_(std::move(curves)), projections_(std::move(projections)) {
    RL_REQUIRE(!curves_.empty(), "CurveSet: no curves");
    for (std::size_t i = 0; i < curves_.size(); ++i) {
        for (std::size_t j = 0; j < i; ++j) RL_REQUIRE(curves_[i].name() != curves_[j].name(), "CurveSet: duplicate curve name " + curves_[i].name());
        RL_REQUIRE(curves_[i].referenceDate() == curves_[0].referenceDate(), "CurveSet: curves have different reference dates");
    }
    discount_ = indexOf(discountCurveName);
    for (auto& [index, curve] : projections_) {
        RL_REQUIRE(!index.empty(), "CurveSet: empty index name in projection map");
        (void)indexOf(curve);
        for (const auto& [other, c2] : projections_)
            RL_REQUIRE(&other == &index || other != index, "CurveSet: index " + index + " mapped twice");
    }
}

template <math::Scalar T>
std::size_t CurveSet<T>::indexOf(std::string_view curveName) const {
    for (std::size_t i = 0; i < curves_.size(); ++i)
        if (curves_[i].name() == curveName) return i;
    throw InputError("CurveSet: no curve named " + std::string(curveName));
}

template <math::Scalar T>
const YieldCurve<T>& CurveSet<T>::projection(std::string_view indexName) const {
    for (const auto& [index, curve] : projections_)
        if (index == indexName) return curves_[indexOf(curve)];
    throw InputError("CurveSet: no projection curve for index " + std::string(indexName));
}

template <math::Scalar T>
bool CurveSet<T>::hasProjection(std::string_view indexName) const noexcept {
    for (const auto& [index, curve] : projections_)
        if (index == indexName) return true;
    return false;
}

template <math::Scalar T>
std::vector<std::string> CurveSet<T>::names() const {
    std::vector<std::string> out;
    out.reserve(curves_.size());
    for (const auto& c : curves_) out.push_back(c.name());
    return out;
}

template <math::Scalar T>
CurveSet<T> CurveSet<T>::withCurve(std::string_view curveName, YieldCurve<T> replacement) const {
    const std::size_t i = indexOf(curveName);
    RL_REQUIRE(replacement.name() == curveName, "CurveSet::withCurve: replacement must keep the name " + std::string(curveName));
    RL_REQUIRE(replacement.referenceDate() == referenceDate(), "CurveSet::withCurve: reference date mismatch");
    CurveSet out = *this;
    out.curves_[i] = std::move(replacement);
    return out;
}

} // namespace rl::market

#define RL_MARKET_INSTANTIATE_CURVES(T)                                                                \
    template class rl::market::YieldCurve<T>;                                                          \
    template class rl::market::CurveSet<T>;
