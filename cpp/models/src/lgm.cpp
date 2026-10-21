#include "rl/models/lgm.hpp"

#include "rl/core/error.hpp"

#include <algorithm>
#include <cmath>

namespace rl::models {

namespace {
/// G(t) = ∫₀ᵗ e^{2κs} ds: ζ grows by σ² · (G(b) − G(a)) over an interval with constant σ.
double G(double kappa, double t) noexcept {
    const double a = 2.0 * kappa * t;
    if (std::fabs(a) < 1e-8) return t * (1.0 + 0.5 * a);
    return std::expm1(a) / (2.0 * kappa);
}
} // namespace

Lgm1F::Lgm1F(double kappa, std::vector<double> zetaTimes, std::vector<double> zetaValues)
    : kappa_(kappa), times_(std::move(zetaTimes)), zetas_(std::move(zetaValues)) {
    RL_REQUIRE(std::isfinite(kappa_), "Lgm1F: kappa must be finite");
    RL_REQUIRE(!times_.empty() && times_.size() == zetas_.size(), "Lgm1F: need matching, non-empty zeta nodes");
    for (std::size_t i = 0; i < times_.size(); ++i) {
        RL_REQUIRE(std::isfinite(times_[i]) && times_[i] > (i ? times_[i - 1] : 0.0), "Lgm1F: zeta times must be positive and strictly increasing");
        RL_REQUIRE(std::isfinite(zetas_[i]) && zetas_[i] > 0.0 && zetas_[i] >= (i ? zetas_[i - 1] : 0.0),
                   "Lgm1F: zeta must be positive and non-decreasing (it is a variance accumulated over time)");
    }
}

Lgm1F Lgm1F::constantSigma(double kappa, double sigma, double horizon) {
    RL_REQUIRE(sigma > 0.0 && horizon > 0.0, "Lgm1F::constantSigma: sigma and horizon must be positive");
    return Lgm1F(kappa, {horizon}, {sigma * sigma * G(kappa, horizon)});
}

double Lgm1F::H(double kappa, double t) noexcept {
    const double a = kappa * t;
    if (std::fabs(a) < 1e-8) return t * (1.0 - 0.5 * a);
    return -std::expm1(-a) / kappa;
}

double Lgm1F::zeta(double t) const {
    RL_REQUIRE(t >= 0.0, "Lgm1F::zeta: time must be non-negative");
    if (t == 0.0) return 0.0;
    // Interval i runs from t_{i−1} to t_i (t_{−1} = 0, ζ_{−1} = 0); beyond the last node use the last interval.
    const std::size_t n = times_.size();
    std::size_t i = static_cast<std::size_t>(std::lower_bound(times_.begin(), times_.end(), t) - times_.begin());
    if (i >= n) i = n - 1;
    const double t0 = i ? times_[i - 1] : 0.0, z0 = i ? zetas_[i - 1] : 0.0;
    const double sigma2 = (zetas_[i] - z0) / (G(kappa_, times_[i]) - G(kappa_, t0));
    return z0 + sigma2 * (G(kappa_, t) - G(kappa_, t0));
}

std::vector<double> Lgm1F::hullWhiteSigmas() const {
    std::vector<double> out(times_.size());
    for (std::size_t i = 0; i < times_.size(); ++i) {
        const double t0 = i ? times_[i - 1] : 0.0, z0 = i ? zetas_[i - 1] : 0.0;
        out[i] = std::sqrt((zetas_[i] - z0) / (G(kappa_, times_[i]) - G(kappa_, t0)));
    }
    return out;
}

double Lgm1F::deflatedBond(double dfT, double T, double zetaT, double x) const noexcept {
    const double h = H(T);
    return dfT * std::exp(-h * x - 0.5 * h * h * zetaT);
}

double Lgm1F::numeraire(double dfT, double t, double zetaT, double x) const noexcept {
    const double h = H(t);
    return std::exp(h * x + 0.5 * h * h * zetaT) / dfT;
}

double Lgm1F::bond(double dfT0, double t, double dfT1, double T, double zetaT, double x) const noexcept {
    return deflatedBond(dfT1, T, zetaT, x) * numeraire(dfT0, t, zetaT, x);
}

} // namespace rl::models
