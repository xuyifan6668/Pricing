#pragma once
// CurveSet<T>: the multi-curve bundle. Named curves plus two *roles* — which curve discounts, and which
// curve projects each index. In USD the SOFR curve plays both the discount role and the USD-SOFR
// projection role and is stored once, so a bootstrap parameter that moves it moves both roles,
// which is what the market means. A cross-currency or CSA curve (v2) is a new role, not a new copy.
//
// Value type; copies the curves (a few KB). Immutable after construction; thread-safe to read.
// DESIGN: a small vector with linear search rather than a map — at most five curves, and the lookup
// happens per cashflow in the bootstrap (1e6 times per solve), where a string hash would cost more
// than the comparison. Lookups by index name return a reference into this object; hold the set, not
// the reference, across a withCurve().

#include "rl/market/yield_curve.hpp"

#include <cstddef>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace rl::market {

template <math::Scalar T>
class CurveSet {
public:
    /// Throws InputError on: no curves, duplicate curve names, unknown discount name, a projection
    /// mapped to an unknown curve, an index mapped twice, or reference dates that differ.
    CurveSet(std::vector<YieldCurve<T>> curves, std::string discountCurveName,
             std::vector<std::pair<std::string /*index name*/, std::string /*curve name*/>> projections);

    core::Date referenceDate() const noexcept { return curves_[0].referenceDate(); }
    const YieldCurve<T>& discount() const noexcept { return curves_[discount_]; }
    const YieldCurve<T>& projection(std::string_view indexName) const; ///< throws InputError if unknown
    bool hasProjection(std::string_view indexName) const noexcept;
    const YieldCurve<T>& curve(std::string_view curveName) const { return curves_[indexOf(curveName)]; }
    const std::vector<YieldCurve<T>>& curves() const noexcept { return curves_; }
    const std::string& discountName() const noexcept { return curves_[discount_].name(); }
    const std::vector<std::pair<std::string, std::string>>& projections() const noexcept { return projections_; }
    std::vector<std::string> names() const; ///< curve names in construction order
    std::size_t size() const noexcept { return curves_.size(); }

    /// Replace one curve (same name, same reference date); roles unchanged. Throws if unknown.
    CurveSet withCurve(std::string_view curveName, YieldCurve<T> replacement) const;

private:
    std::size_t indexOf(std::string_view curveName) const; // throws InputError
    std::vector<YieldCurve<T>> curves_;
    std::vector<std::pair<std::string, std::string>> projections_;
    std::size_t discount_ = 0;
};

} // namespace rl::market
