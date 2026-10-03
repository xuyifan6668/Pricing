#pragma once
// Pillar selection: one free node per instrument, at the instrument's maturity, in ACT/365F time
// from the as-of date. Turn-of-year instruments do not add nodes; they add a jump parameter on the
// curve's jump layer (YieldCurve::turns), which is why CurveSpec carries them separately.
//
// The CurveSpec is the immutable part of a curve — everything except the ordinates — and is what
// BootstrapProblem::curves(z) needs to rebuild a YieldCurve<T> from a slice of z.

#include "rl/bootstrap/instrument_set.hpp"
#include "rl/core/date.hpp"
#include "rl/market/curve_interpolation.hpp"

#include <map>
#include <string>
#include <vector>

namespace rl::bootstrap {

struct CurveSpec {
    std::string name;                          ///< curve name (== the index it is named after)
    std::vector<std::string> projects;         ///< index names this curve projects (at least its own)
    bool isDiscount = false;
    market::CurveInterpolation interpolation = market::CurveInterpolation::LinearLogDf;
    std::vector<core::Date> pillarDates;       ///< strictly increasing, all after as-of
    std::vector<double> pillarTimes;           ///< ACT/365F from as-of
    std::vector<std::string> pillarIds;        ///< the instrument pinning each pillar
    std::vector<TurnOfYear> turns;             ///< jump parameters appended after the nodes in z (if any)

    std::size_t parameterCount() const noexcept { return pillarTimes.size() + turns.size(); }
};

struct PillarOptions {
    std::string discountCurve = "USD-SOFR";
    market::CurveInterpolation interpolation = market::CurveInterpolation::LinearLogDf;
    std::map<std::string, market::CurveInterpolation> perCurve; ///< overrides by curve name
};

/// One CurveSpec per curve named in the instrument set, in first-appearance order. Throws InputError if
/// the discount curve is not among them or a curve has no pillars.
std::vector<CurveSpec> makeCurveSpecs(const InstrumentSet& instruments, core::Date asOf, const PillarOptions& options = {});

/// Total parameter count over all specs (== instrument count when every instrument pins one pillar).
std::size_t parameterCount(const std::vector<CurveSpec>& specs) noexcept;

} // namespace rl::bootstrap
