#include "rl/bootstrap/pillars.hpp"

#include "rl/core/daycount.hpp"
#include "rl/core/error.hpp"

#include <algorithm>

namespace rl::bootstrap {

std::vector<CurveSpec> makeCurveSpecs(const InstrumentSet& instruments, core::Date asOf, const PillarOptions& options) {
    std::vector<CurveSpec> specs;
    bool haveDiscount = false;
    for (const std::string& name : curveNames(instruments)) {
        CurveSpec s;
        s.name = name;
        s.projects = {name}; // a curve projects the index it is named after; instruments assigned to it may be on other indices
        s.isDiscount = (name == options.discountCurve);
        haveDiscount = haveDiscount || s.isDiscount;
        auto it = options.perCurve.find(name);
        s.interpolation = it == options.perCurve.end() ? options.interpolation : it->second;
        for (const BootstrapInstrument* bi : instrumentsOn(instruments, name)) {
            if (const auto* turn = std::get_if<TurnOfYear>(&bi->trade)) {
                s.turns.push_back(*turn);
                continue;
            }
            s.pillarDates.push_back(bi->maturity);
            s.pillarTimes.push_back(core::yearFraction(core::DayCount::Act365Fixed, asOf, bi->maturity));
            s.pillarIds.push_back(bi->id);
        }
        RL_REQUIRE(!s.pillarTimes.empty(), "curve " + name + " has no pillar instruments");
        for (std::size_t i = 0; i < s.pillarTimes.size(); ++i) {
            RL_REQUIRE(s.pillarTimes[i] > 0.0, "curve " + name + ": pillar " + s.pillarIds[i] + " is not after as-of");
            if (i > 0) RL_REQUIRE(s.pillarTimes[i] > s.pillarTimes[i - 1], "curve " + name + ": pillars not strictly increasing at " + s.pillarIds[i]);
        }
        specs.push_back(std::move(s));
    }
    RL_REQUIRE(haveDiscount, "discount curve " + options.discountCurve + " has no instruments");
    return specs;
}

std::size_t parameterCount(const std::vector<CurveSpec>& specs) noexcept {
    std::size_t n = 0;
    for (const auto& s : specs) n += s.parameterCount();
    return n;
}

} // namespace rl::bootstrap
