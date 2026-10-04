#include "rl/bootstrap/report.hpp"

#include "rl/bootstrap/residuals.hpp"
#include "rl/market/hash.hpp"

#include <nlohmann/json.hpp>

#include <cmath>
#include <cstdio>
#include <sstream>

namespace rl::bootstrap {

RepricingReport reprice(const BootstrapProblem& problem, const math::Vector& z, std::uint64_t snapshotHash) {
    RepricingReport rep;
    rep.asOf = problem.asOf();
    rep.snapshotHash = snapshotHash;
    const market::CurveSet<double> curves = problem.curves(z);
    const math::Vector raw = problem.residuals(z);
    double sum = 0.0;
    std::size_t i = 0;
    for (const BootstrapInstrument& bi : problem.instruments()) {
        RepricingRow row;
        row.id = bi.id;
        row.type = bi.type;
        row.curve = bi.curveName;
        row.maturity = bi.maturity;
        row.quote = bi.source.value;
        const double model = modelQuote(bi, curves, problem.fixings());
        row.model = bi.source.fromRate(model);
        row.residual = row.model - row.quote;
        // In bp of rate: price residuals are in price points (1 point = 100 bp of rate).
        row.residualBp = bi.source.unit == market::QuoteUnit::Price ? row.residual * 100.0 : (model - bi.source.asRate()) * 1e4;
        row.rawResidual = i < static_cast<std::size_t>(raw.size()) ? raw[static_cast<Eigen::Index>(i)] : std::nan("");
        rep.maxAbsBp = std::max(rep.maxAbsBp, std::fabs(row.residualBp));
        sum += std::fabs(row.residualBp);
        rep.rows.push_back(std::move(row));
        ++i;
    }
    rep.meanAbsBp = rep.rows.empty() ? 0.0 : sum / static_cast<double>(rep.rows.size());
    return rep;
}

std::string RepricingReport::toCsv() const {
    std::ostringstream os;
    os << "id,type,curve,maturity,quote,model,residual,residual_bp,raw_residual\n";
    char buf[256];
    for (const auto& r : rows) {
        std::snprintf(buf, sizeof buf, "%s,%s,%s,%s,%.10g,%.12g,%.6e,%.6e,%.6e\n", r.id.c_str(), std::string(market::toString(r.type)).c_str(), r.curve.c_str(),
                      r.maturity.toIso().c_str(), r.quote, r.model, r.residual, r.residualBp, r.rawResidual);
        os << buf;
    }
    return os.str();
}

std::string RepricingReport::toJson() const {
    nlohmann::ordered_json j;
    j["asOf"] = asOf.toIso();
    j["snapshotHash"] = market::toHex(snapshotHash);
    j["maxAbsBp"] = maxAbsBp;
    j["meanAbsBp"] = meanAbsBp;
    auto arr = nlohmann::ordered_json::array();
    for (const auto& r : rows) {
        nlohmann::ordered_json e;
        e["id"] = r.id;
        e["type"] = std::string(market::toString(r.type));
        e["curve"] = r.curve;
        e["maturity"] = r.maturity.toIso();
        e["quote"] = r.quote;
        e["model"] = r.model;
        e["residual"] = r.residual;
        e["residualBp"] = r.residualBp;
        e["rawResidual"] = r.rawResidual;
        arr.push_back(std::move(e));
    }
    j["rows"] = std::move(arr);
    return j.dump(2);
}

} // namespace rl::bootstrap
