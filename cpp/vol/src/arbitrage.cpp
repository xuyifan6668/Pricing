// Breeden–Litzenberger density diagnostic (sources and numerical choices in arbitrage.hpp).
#include "rl/vol/arbitrage.hpp"

#include "rl/core/error.hpp"
#include "rl/vol/bachelier.hpp"
#include "rl/vol/sabr.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace rl::vol {

std::vector<double> strikeGrid(double kMin, double kMax, std::size_t n) {
    RL_REQUIRE(n >= 2 && std::isfinite(kMin) && std::isfinite(kMax) && kMax > kMin, "strikeGrid: need n >= 2 and kMax > kMin");
    std::vector<double> ks(n);
    const double step = (kMax - kMin) / static_cast<double>(n - 1);
    for (std::size_t i = 0; i < n; ++i) ks[i] = kMin + step * static_cast<double>(i);
    ks.back() = kMax;
    return ks;
}

std::pair<double, double> diagnosticStrikeRange(const market::SabrParams<double>& p, double forward, double expiry, double atmNormalVol, double nStdDev,
                                                double floor) {
    RL_REQUIRE(expiry > 0.0 && atmNormalVol > 0.0 && nStdDev > 0.0 && floor > 0.0, "diagnosticStrikeRange: inputs must be positive");
    const double width = nStdDev * atmNormalVol * std::sqrt(expiry);
    double lo = forward - width;
    const double hi = forward + width;
    if (p.beta > 0.0) lo = std::max(lo, -p.shift + floor);
    RL_REQUIRE(hi > lo, "diagnosticStrikeRange: empty range (forward + shift below the floor?)");
    return {lo, hi};
}

std::vector<DensityPoint> riskNeutralDensity(const market::SabrParams<double>& p, double forward, double expiry, const std::vector<double>& strikes, double h) {
    RL_REQUIRE(std::isfinite(expiry) && expiry > 0.0, "riskNeutralDensity: expiry must be positive");
    RL_REQUIRE(std::isfinite(h) && h > 0.0, "riskNeutralDensity: step must be positive");
    const double nan = std::numeric_limits<double>::quiet_NaN();
    std::vector<DensityPoint> out;
    out.reserve(strikes.size());
    for (double k : strikes) {
        if (p.beta > 0.0) RL_REQUIRE(k - h + p.shift > 0.0, "riskNeutralDensity: strike " + std::to_string(k) + " is within h of -shift with beta > 0");
        const double vm = sabrNormalVol<double>(p, forward, k - h, expiry);
        const double v0 = sabrNormalVol<double>(p, forward, k, expiry);
        const double vp = sabrNormalVol<double>(p, forward, k + h, expiry);
        DensityPoint d{k, v0, nan, nan, nan};
        if (vm > 0.0 && v0 > 0.0 && vp > 0.0) {
            const double cm = bachelierPrice<double>(forward, k - h, expiry, vm, OptionType::Call);
            const double c0 = bachelierPrice<double>(forward, k, expiry, v0, OptionType::Call);
            const double cp = bachelierPrice<double>(forward, k + h, expiry, vp, OptionType::Call);
            d.callPrice = c0;
            d.dCalldK = (cp - cm) / (2.0 * h);
            d.density = (cp - 2.0 * c0 + cm) / (h * h);
        }
        out.push_back(d);
    }
    return out;
}

ArbitrageCheck checkDensity(const std::vector<DensityPoint>& density, double forward, double densityTol) {
    RL_REQUIRE(density.size() >= 2, "checkDensity: need at least two points");
    for (std::size_t i = 1; i < density.size(); ++i) RL_REQUIRE(density[i].strike > density[i - 1].strike, "checkDensity: strikes must be ascending");
    ArbitrageCheck c;
    c.kMin = density.front().strike;
    c.kMax = density.back().strike;
    c.points = density.size();
    c.minDensity = std::numeric_limits<double>::infinity();
    c.minDensityStrike = std::numeric_limits<double>::quiet_NaN();
    auto bad = [&](const DensityPoint& d) { return !(d.normalVol > 0.0) || std::isnan(d.density) || d.density < -densityTol; };

    for (std::size_t i = 0; i < density.size(); ++i) {
        const auto& d = density[i];
        if (!(d.normalVol > 0.0)) c.nonPositiveVol = true;
        if (!std::isnan(d.density) && d.density < c.minDensity) {
            c.minDensity = d.density;
            c.minDensityStrike = d.strike;
        }
        if (!std::isnan(d.dCalldK) && (d.dCalldK > 1e-10 || d.dCalldK < -1.0 - 1e-10)) c.callSpreadViolated = true;
        if (i > 0) {
            const auto& a = density[i - 1];
            const double w = d.strike - a.strike;
            const double pa = std::isnan(a.density) ? 0.0 : a.density, pb = std::isnan(d.density) ? 0.0 : d.density;
            c.integratedMass += 0.5 * w * (pa + pb);
            c.impliedMean += 0.5 * w * (pa * a.strike + pb * d.strike);
        }
    }
    // Scan away from the forward in both directions.
    const auto firstAbove = std::lower_bound(density.begin(), density.end(), forward, [](const DensityPoint& d, double f) { return d.strike < f; });
    for (auto it = std::make_reverse_iterator(firstAbove); it != density.rend(); ++it) {
        if (bad(*it)) {
            c.firstNegativeDensityBelowForward = it->strike;
            break;
        }
    }
    for (auto it = firstAbove; it != density.end(); ++it) {
        if (bad(*it)) {
            c.firstNegativeDensityAboveForward = it->strike;
            break;
        }
    }
    return c;
}

ArbitrageCheck checkSabrSmile(const market::SabrParams<double>& p, double forward, double expiry, double kMin, double kMax, std::size_t n) {
    return checkDensity(riskNeutralDensity(p, forward, expiry, strikeGrid(kMin, kMax, n)), forward);
}

} // namespace rl::vol
