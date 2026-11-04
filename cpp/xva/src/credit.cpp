#include "rl/xva/credit.hpp"

#include "rl/core/error.hpp"
#include "rl/math/rootfind.hpp"

#include <algorithm>
#include <cmath>

namespace rl::xva {

namespace {
struct Legs {
    double premiumPerSpread = 0.0, protectionPerLoss = 0.0;
};
/// Both legs on the quarterly grid out to `maturity`, per unit spread and per unit loss given default.
Legs cdsLegs(const HazardCurve& curve, double maturity, const std::function<double(double)>& discount) {
    Legs l;
    const int n = std::max(1, static_cast<int>(std::lround(maturity * 4.0)));
    const double dt = maturity / n;
    double qPrev = 1.0;
    for (int j = 1; j <= n; ++j) {
        const double t = dt * j, q = curve.survival(t), d = discount(t);
        l.premiumPerSpread += dt * d * q + 0.5 * dt * d * (qPrev - q);
        l.protectionPerLoss += d * (qPrev - q);
        qPrev = q;
    }
    return l;
}
} // namespace

HazardCurve::HazardCurve(std::vector<double> times, std::vector<double> hazards) : times_(std::move(times)), hazards_(std::move(hazards)) {
    RL_REQUIRE(!times_.empty() && times_.size() == hazards_.size(), "HazardCurve: need matching, non-empty nodes");
    for (std::size_t i = 0; i < times_.size(); ++i) {
        RL_REQUIRE(std::isfinite(times_[i]) && times_[i] > (i ? times_[i - 1] : 0.0), "HazardCurve: times must be positive and strictly increasing");
        RL_REQUIRE(std::isfinite(hazards_[i]) && hazards_[i] >= 0.0, "HazardCurve: hazard rates must be non-negative");
    }
}

HazardCurve HazardCurve::flat(double hazard) { return HazardCurve({1.0}, {hazard}); }

double HazardCurve::survival(double t) const {
    RL_REQUIRE(t >= 0.0, "HazardCurve::survival: time must be non-negative");
    double integral = 0.0, from = 0.0;
    for (std::size_t i = 0; i < times_.size() && from < t; ++i) {
        const double to = i + 1 == times_.size() ? t : std::min(t, times_[i]);
        integral += hazards_[i] * (to - from);
        from = to;
    }
    return std::exp(-integral);
}

double cdsParSpread(const HazardCurve& curve, double maturity, double recovery, const std::function<double(double)>& discount) {
    RL_REQUIRE(maturity > 0.0, "cdsParSpread: maturity must be positive");
    const Legs l = cdsLegs(curve, maturity, discount);
    return (1.0 - recovery) * l.protectionPerLoss / l.premiumPerSpread;
}

HazardCurve HazardCurve::bootstrap(const std::vector<CdsQuote>& quotes, double recovery, const std::function<double(double)>& discount) {
    RL_REQUIRE(!quotes.empty(), "HazardCurve::bootstrap: no quotes");
    RL_REQUIRE(recovery >= 0.0 && recovery < 1.0, "HazardCurve::bootstrap: recovery must be in [0, 1)");
    std::vector<double> times, hazards;
    for (const CdsQuote& q : quotes) {
        RL_REQUIRE(q.maturityYears > (times.empty() ? 0.0 : times.back()), "HazardCurve::bootstrap: quotes must be in increasing maturity");
        RL_REQUIRE(q.spread > 0.0, "HazardCurve::bootstrap: spreads must be positive");
        times.push_back(q.maturityYears);
        hazards.push_back(0.0);
        // The CDS's value to the protection buyer rises with the hazard rate on the last interval.
        auto value = [&](double h) {
            hazards.back() = h;
            const Legs l = cdsLegs(HazardCurve(times, hazards), q.maturityYears, discount);
            return (1.0 - recovery) * l.protectionPerLoss - q.spread * l.premiumPerSpread;
        };
        if (value(0.0) > 0.0)
            throw CalibrationError("HazardCurve::bootstrap: the " + std::to_string(q.maturityYears) + "y spread implies a negative hazard rate after the earlier maturities", {value(0.0)});
        double hi = 2.0 * q.spread / (1.0 - recovery) + 0.01;
        while (value(hi) < 0.0) {
            hi *= 2.0;
            RL_REQUIRE(hi < 100.0, "HazardCurve::bootstrap: spread is not attainable");
        }
        math::RootOptions ro;
        ro.xTol = 1e-15;
        ro.maxIter = 200;
        hazards.back() = math::brent(value, 0.0, hi, ro).x;
    }
    return HazardCurve(std::move(times), std::move(hazards));
}

} // namespace rl::xva
