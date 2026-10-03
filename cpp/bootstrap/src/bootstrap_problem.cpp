#include "rl/bootstrap/bootstrap_problem.hpp"

#include "rl/aad/real.hpp"
#include "rl/bootstrap/residuals.hpp"
#include "rl/core/error.hpp"

#include <algorithm>

namespace rl::bootstrap {

BootstrapProblem::BootstrapProblem(core::Date asOf, std::vector<CurveSpec> curveSpecs, InstrumentSet instruments, market::FixingStore fixings,
                                   market::QuoteSet quotes)
    : asOf_(asOf), specs_(std::move(curveSpecs)), instruments_(std::move(instruments)), fixings_(std::move(fixings)), quotes_(std::move(quotes)) {
    RL_REQUIRE(!specs_.empty(), "BootstrapProblem: no curve specs");
    for (const auto& s : specs_) nParams_ += s.parameterCount();
    RL_REQUIRE(nParams_ == instruments_.size(), "BootstrapProblem: " + std::to_string(nParams_) + " parameters for " + std::to_string(instruments_.size()) +
                                                    " instruments (the system must be square)");
    for (const auto& bi : instruments_) {
        RL_REQUIRE(std::any_of(specs_.begin(), specs_.end(), [&](const CurveSpec& s) { return s.name == bi.curveName; }),
                   "BootstrapProblem: instrument " + bi.id + " pins curve " + bi.curveName + " which has no spec");
        ids_.push_back(bi.id);
    }
}

std::pair<std::size_t, std::size_t> BootstrapProblem::parameterRange(std::string_view curveName) const {
    std::size_t offset = 0;
    for (const auto& s : specs_) {
        if (s.name == curveName) return {offset, s.parameterCount()};
        offset += s.parameterCount();
    }
    throw InputError("BootstrapProblem: no curve named " + std::string(curveName));
}

template <math::Scalar T>
market::CurveSet<T> BootstrapProblem::curves(const std::vector<T>& z) const {
    RL_REQUIRE(z.size() == nParams_, "BootstrapProblem::curves: z has " + std::to_string(z.size()) + " entries, expected " + std::to_string(nParams_));
    std::vector<market::YieldCurve<T>> cs;
    cs.reserve(specs_.size());
    std::string discount;
    std::vector<std::pair<std::string, std::string>> projections;
    std::size_t off = 0;
    for (const auto& s : specs_) {
        const std::size_t nNodes = s.pillarTimes.size();
        std::vector<T> ys(z.begin() + static_cast<std::ptrdiff_t>(off), z.begin() + static_cast<std::ptrdiff_t>(off + nNodes));
        std::vector<market::TurnJump<T>> turns;
        for (std::size_t k = 0; k < s.turns.size(); ++k) {
            const TurnOfYear& t = s.turns[k];
            turns.push_back({core::yearFraction(core::DayCount::Act365Fixed, asOf_, t.start), core::yearFraction(core::DayCount::Act365Fixed, asOf_, t.end),
                             z[off + nNodes + k]});
        }
        off += s.parameterCount();
        cs.emplace_back(asOf_, s.name, s.pillarTimes, std::move(ys), s.interpolation, std::move(turns));
        if (s.isDiscount) discount = s.name;
        for (const auto& index : s.projects) projections.emplace_back(index, s.name);
    }
    return market::CurveSet<T>(std::move(cs), std::move(discount), std::move(projections));
}

market::CurveSet<double> BootstrapProblem::curves(const math::Vector& z) const {
    return curves(std::vector<double>(z.data(), z.data() + z.size()));
}

math::Vector BootstrapProblem::quoteValues() const {
    math::Vector q(static_cast<Eigen::Index>(instruments_.size()));
    for (std::size_t i = 0; i < instruments_.size(); ++i) q[static_cast<Eigen::Index>(i)] = instruments_[i].quote;
    return q;
}

BootstrapProblem BootstrapProblem::withQuotes(const math::Vector& q) const {
    RL_REQUIRE(static_cast<std::size_t>(q.size()) == instruments_.size(), "withQuotes: size mismatch");
    std::vector<market::Quote> qs = quotes_.quotes();
    for (std::size_t i = 0; i < instruments_.size(); ++i) {
        const auto& bi = instruments_[i];
        qs[quotes_.indexOf(bi.id)].value = bi.source.fromRate(q[static_cast<Eigen::Index>(i)]);
    }
    market::QuoteSet newQuotes(quotes_.asOf(), quotes_.currency(), std::move(qs));
    InstrumentSetOptions io;
    io.includeTurns = std::any_of(specs_.begin(), specs_.end(), [](const CurveSpec& s) { return !s.turns.empty(); });
    if (!instruments_.empty()) {
        // Keep the notional the original set was built with (per-unit by default).
        if (const auto* s = std::get_if<instruments::Swap>(&instruments_.front().trade))
            if (!s->fixedLeg().empty())
                if (const auto* c = std::get_if<instruments::FixedCoupon>(&s->fixedLeg().cashflows().front())) io.notional = c->notional;
    }
    InstrumentSet rebuilt = buildInstrumentSet(newQuotes, io);
    // Specs are unchanged (same pillars, same turns); only the quotes and trades moved.
    return BootstrapProblem(asOf_, specs_, std::move(rebuilt), fixings_, std::move(newQuotes));
}

math::Vector BootstrapProblem::residuals(const math::Vector& z) const {
    const market::CurveSet<double> cs = curves(z);
    math::Vector r(static_cast<Eigen::Index>(instruments_.size()));
    for (std::size_t i = 0; i < instruments_.size(); ++i) r[static_cast<Eigen::Index>(i)] = residual(instruments_[i], cs, fixings_);
    return r;
}

math::Vector BootstrapProblem::scaling(const math::Vector& z) const {
    const market::CurveSet<double> cs = curves(z);
    math::Vector s(static_cast<Eigen::Index>(instruments_.size()));
    for (std::size_t i = 0; i < instruments_.size(); ++i) s[static_cast<Eigen::Index>(i)] = residualScale(instruments_[i], cs, fixings_);
    return s;
}

template market::CurveSet<double> BootstrapProblem::curves<double>(const std::vector<double>&) const;
template market::CurveSet<rl::aad::Real> BootstrapProblem::curves<rl::aad::Real>(const std::vector<rl::aad::Real>&) const;

} // namespace rl::bootstrap
