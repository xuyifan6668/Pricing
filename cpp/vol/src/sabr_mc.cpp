#include "rl/vol/sabr_mc.hpp"

#include "rl/core/error.hpp"
#include "rl/math/rng.hpp"
#include "rl/vol/bachelier.hpp"
#include "rl/vol/sabr.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace rl::vol {

SabrMcResult sabrMonteCarlo(const market::SabrParams<double>& p, double forward, double expiry, const std::vector<double>& strikes, const SabrMcOptions& opt) {
    market::validateSabrParams(p);
    RL_REQUIRE(expiry > 0.0, "sabrMonteCarlo: expiry must be positive");
    RL_REQUIRE(opt.paths >= 1000 && opt.paths % 2 == 0 && opt.stepsPerYear >= 1, "sabrMonteCarlo: need an even number of at least 1000 paths and one step a year");
    const double f0 = forward + p.shift;
    RL_REQUIRE(p.beta == 0.0 || f0 > 0.0, "sabrMonteCarlo: forward + shift must be positive when beta > 0");

    const std::size_t n = opt.paths, half = n / 2;
    const int steps = std::max(1, static_cast<int>(std::ceil(expiry * opt.stepsPerYear - 1e-9)));
    const double dt = expiry / steps, root = std::sqrt(dt), rhoBar = std::sqrt(1.0 - p.rho * p.rho);
    const double volDrift = -0.5 * p.nu * p.nu * dt;
    std::vector<double> f(n, f0), a(n, p.alpha);
    std::vector<char> dead(n, 0);
    auto local = [&](double x) {
        if (p.beta == 0.0) return 1.0;
        if (p.beta == 0.5) return std::sqrt(x);
        if (p.beta == 1.0) return x;
        return std::pow(x, p.beta);
    };
    math::Rng rng(opt.seed);
    for (int k = 0; k < steps; ++k) {
        for (std::size_t q = 0; q < half; ++q) {
            const double z1 = rng.gaussian(), z2 = rng.gaussian();
            for (int side = 0; side < 2; ++side) {
                const std::size_t i = q + static_cast<std::size_t>(side) * half;
                const double s = side ? -1.0 : 1.0;
                const double dz = s * z2 * root, dw = p.rho * dz + rhoBar * s * z1 * root;
                if (!dead[i]) {
                    double x = f[i] + a[i] * local(f[i]) * dw;
                    if (p.beta > 0.0 && x <= 0.0) {
                        x = 0.0;
                        dead[i] = 1;
                    }
                    f[i] = x;
                }
                a[i] *= std::exp(p.nu * dz + volDrift);
            }
        }
    }

    SabrMcResult out;
    double meanF = 0.0, varF = 0.0;
    for (double x : f) meanF += x;
    meanF /= static_cast<double>(n);
    for (double x : f) varF += (x - meanF) * (x - meanF);
    varF /= static_cast<double>(n);
    out.forwardError = meanF - f0;
    out.absorbed = static_cast<double>(std::count(dead.begin(), dead.end(), 1)) / static_cast<double>(n);

    for (double strikeRaw : strikes) {
        SabrMcPoint pt;
        pt.strike = strikeRaw;
        const double k = strikeRaw + p.shift;
        const bool call = strikeRaw >= forward;
        // Payoff per path, and its covariance with F_T for the control variate.
        double meanY = 0.0, cov = 0.0;
        std::vector<double> y(n);
        for (std::size_t i = 0; i < n; ++i) {
            y[i] = call ? std::max(f[i] - k, 0.0) : std::max(k - f[i], 0.0);
            meanY += y[i];
        }
        meanY /= static_cast<double>(n);
        for (std::size_t i = 0; i < n; ++i) cov += (y[i] - meanY) * (f[i] - meanF);
        cov /= static_cast<double>(n);
        const double c = varF > 0.0 ? cov / varF : 0.0;
        // Pairwise averages of the controlled payoff give the standard error under antithetic sampling.
        double sum = 0.0, sumSq = 0.0;
        for (std::size_t q = 0; q < half; ++q) {
            const double v = 0.5 * ((y[q] - c * (f[q] - f0)) + (y[q + half] - c * (f[q + half] - f0)));
            sum += v;
            sumSq += v * v;
        }
        const double pairs = static_cast<double>(half);
        pt.price = sum / pairs;
        pt.standardError = std::sqrt(std::max(sumSq / pairs - pt.price * pt.price, 0.0) / pairs);
        pt.haganNormalVol = sabrNormalVol<double>(p, forward, strikeRaw, expiry);
        pt.normalVol = std::numeric_limits<double>::quiet_NaN();
        try {
            const OptionType type = call ? OptionType::Call : OptionType::Put;
            pt.normalVol = bachelierImpliedVol(pt.price, forward, strikeRaw, expiry, type);
            const double bump = 1e-6;
            const double vega = (bachelierPrice<double>(forward, strikeRaw, expiry, pt.normalVol + bump, type) - bachelierPrice<double>(forward, strikeRaw, expiry, pt.normalVol - bump, type)) / (2.0 * bump);
            pt.normalVolError = vega > 0.0 ? pt.standardError / vega : 0.0;
        } catch (const Error&) {
        }
        out.points.push_back(pt);
    }
    return out;
}

} // namespace rl::vol
