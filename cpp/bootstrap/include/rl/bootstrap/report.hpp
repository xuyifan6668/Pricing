#pragma once
// Repricing report: the calibration report's first table. For every instrument: quote, model
// value in the same units, residual in quote units and in basis points. Stamped with the snapshot hash
// so the report can be tied to the market it came from .

#include "rl/bootstrap/bootstrap_problem.hpp"
#include "rl/core/date.hpp"
#include "rl/market/quote.hpp"
#include "rl/math/lm.hpp"

#include <cstdint>
#include <string>
#include <vector>

namespace rl::bootstrap {

struct RepricingRow {
    std::string id;
    market::QuoteType type;
    std::string curve;
    core::Date maturity;
    double quote;        ///< as quoted (pct / bp / price)
    double model;        ///< model value in the same units
    double residual;     ///< model − quote, quoted units
    double residualBp;   ///< residual in basis points of rate (price residual × 100)
    double rawResidual;  ///< the solver's residual r_i(z) (whatever form you chose)
};

struct RepricingReport {
    core::Date asOf;
    std::uint64_t snapshotHash = 0;
    std::vector<RepricingRow> rows;
    double maxAbsBp = 0.0;
    double meanAbsBp = 0.0;

    std::string toCsv() const;   ///< header: id,type,curve,maturity,quote,model,residual,residual_bp,raw_residual
    std::string toJson() const;
};

/// Evaluate the report at parameter vector z. `snapshotHash` is passed through (0 if unknown).
RepricingReport reprice(const BootstrapProblem& problem, const math::Vector& z, std::uint64_t snapshotHash = 0);

} // namespace rl::bootstrap
