#pragma once
// BootstrapProblem: the calibration-problem abstraction instantiated for the global multi-curve
// bootstrap (docs/DESIGN.md "Calibration-problem abstraction").
//
// Every calibration in this library has the same shape: a parameter vector z, a residual vector
// F(z; q) that vanishes when the model reprices the instruments quoted at q, a Jacobian ∂F/∂z, a
// scaling that makes residuals commensurable, and instrument ids for the report. Three instances:
//    BootstrapProblem   z = all free nodes (+ turn jumps) of all curves; one residual per
//                               instrument; square; Newton (global_bootstrap.hpp).
//    SabrNodeProblem    z = (α, ρ, ν) or (ρ, ν); vol residuals at 5–7 strikes; m > n; LM (vol/sabr_fit.hpp).
//    LgmStripProblem    z = ζ increments; one swaption each; sequential 1-D roots.
// DESIGN: concrete classes satisfying the CalibrationProblem concept (math/calibration_problem.hpp), not a virtual base. The
// solvers are double-only free functions templated on the problem; nothing is dispatched per
// path, and Python binds the concrete BootstrapProblem, so a virtual base would exist only
// to be subclassed from Python.
//
// DESIGN: the quotes are an explicit input (q lives in the problem, not in closures), so that ∂F/∂q
// exists and par risk can be computed by the implicit-function theorem. withQuotes() rebuilds
// the trades whose quote is a contract term (a swap's fixed rate is written into its coupons).
//
// Layout of z: curves in curveSpecs() order; within a curve the pillar ordinates y = −log P in pillar
// order, then that curve's turn jumps. parameterRange(name) gives the slice. reads adjoints
// back by this layout.
//
// double-only apart from curves<T>(), which exists so a caller can rebuild CurveSet<Real> from
// Real::input leaves for the AAD-versus-analytic Jacobian check.

#include "rl/bootstrap/instrument_set.hpp"
#include "rl/bootstrap/pillars.hpp"
#include "rl/core/date.hpp"
#include "rl/market/curve_set.hpp"
#include "rl/market/fixings.hpp"
#include "rl/market/quote.hpp"
#include "rl/math/calibration_problem.hpp"
#include "rl/math/lm.hpp"

#include <concepts>
#include <cstddef>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace rl::bootstrap {

// The concept moved to math (2026-09-13) so vol and models can satisfy it without linking bootstrap.
using math::CalibrationProblem;

class BootstrapProblem {
public:
    /// Throws InputError unless Σ parameters over the specs equals the instrument count (square), and
    /// every instrument's curve has a spec.
    BootstrapProblem(core::Date asOf, std::vector<CurveSpec> curveSpecs, InstrumentSet instruments, market::FixingStore fixings, market::QuoteSet quotes);

    std::size_t size() const noexcept { return instruments_.size(); }
    std::size_t parameterCount() const noexcept { return nParams_; }
    std::size_t residualCount() const noexcept { return instruments_.size(); }
    const std::vector<CurveSpec>& curveSpecs() const noexcept { return specs_; }
    const InstrumentSet& instruments() const noexcept { return instruments_; }
    const market::QuoteSet& quotes() const noexcept { return quotes_; }
    const market::FixingStore& fixings() const noexcept { return fixings_; }
    core::Date asOf() const noexcept { return asOf_; }
    const std::vector<std::string>& instrumentIds() const noexcept { return ids_; }
    /// (offset, count) of a curve's parameters in z (nodes then turns). Throws InputError if unknown.
    std::pair<std::size_t, std::size_t> parameterRange(std::string_view curveName) const;

    /// Rebuild every curve from z. Pure function of z and the fixed specs.
    market::CurveSet<double> curves(const math::Vector& z) const;
    template <math::Scalar T>
    market::CurveSet<T> curves(const std::vector<T>& z) const;

    /// q, aligned with instrumentIds() (decimal rates / spreads / jumps; prices for futures).
    math::Vector quoteValues() const;
    /// Same instruments and specs, new quote values: trades are rebuilt so that a contract term (fixed
    /// rate, spread, price) reflects the new quote.
    BootstrapProblem withQuotes(const math::Vector& q) const;

    /// F(z; q): one call to bootstrap::residual per instrument, in instrument order.
    math::Vector residuals(const math::Vector& z) const;
    /// Per-residual scale to basis points of quote (bootstrap::residualScale), evaluated at z.
    math::Vector scaling(const math::Vector& z) const;

private:
    core::Date asOf_;
    std::vector<CurveSpec> specs_;
    InstrumentSet instruments_;
    market::FixingStore fixings_;
    market::QuoteSet quotes_;
    std::vector<std::string> ids_;
    std::size_t nParams_ = 0;
};

static_assert(CalibrationProblem<BootstrapProblem>);

} // namespace rl::bootstrap
