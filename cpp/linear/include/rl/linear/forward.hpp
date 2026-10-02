#pragma once
// Forward rates off one curve.
//   simpleForward(t1, t2, τ)     = (P(t1)/P(t2) − 1)/τ            money-market (and OIS-compounded) forward
//   continuousForward(t1, t2)    = log(P(t1)/P(t2)) / (t2 − t1)    average instantaneous forward in curve time
// The date overloads compute τ under the given day count and t under the curve's ACT/365F convention.

#include "rl/core/date.hpp"
#include "rl/core/daycount.hpp"
#include "rl/market/yield_curve.hpp"
#include "rl/math/scalar.hpp"

namespace rl::linear {

template <math::Scalar T>
T simpleForward(const market::YieldCurve<T>& curve, double t1, double t2, double tau);
template <math::Scalar T>
T simpleForward(const market::YieldCurve<T>& curve, core::Date d1, core::Date d2, core::DayCount dc);
template <math::Scalar T>
T continuousForward(const market::YieldCurve<T>& curve, double t1, double t2);

} // namespace rl::linear
