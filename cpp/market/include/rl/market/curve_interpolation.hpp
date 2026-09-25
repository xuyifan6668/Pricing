#pragma once
// Runtime-selectable interpolation kind for a yield curve (config plumbing). YieldCurve maps
// each kind to a math::Interpolator<T> alternative applied to y(t) = −log P(0,t):
//   LinearLogDf     math::Linear<T> on y      → piecewise-flat instantaneous forwards; the OIS front-end default
//   MonotoneConvex  math::MonotoneConvex<T>   → Hagan–West 2006; continuous forwards, local, positivity-preserving
// The strings are what data/quotes/*.json and the Excel layer use.

#include <string_view>

namespace rl::market {

enum class CurveInterpolation { LinearLogDf, MonotoneConvex };

/// "linear_logdf", "monotone_convex"
std::string_view toString(CurveInterpolation k);
/// Accepts the strings above plus "linear", "loglinear", "log_linear", "mc", "monotoneconvex"
/// (case-insensitive). Throws InputError.
CurveInterpolation curveInterpolationFromString(std::string_view s);

} // namespace rl::market
