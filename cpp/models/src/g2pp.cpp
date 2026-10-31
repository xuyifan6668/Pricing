#include "rl/models/g2pp.hpp"

#include "rl/core/error.hpp"

#include <algorithm>
#include <cmath>

namespace rl::models {

namespace {
/// ∫₀ᵗ e^{ku} du.
double E(double k, double t) noexcept {
    const double a = k * t;
    if (std::fabs(a) < 1e-8) return t * (1.0 + 0.5 * a);
    return std::expm1(a) / k;
}
double H(double kappa, double t) noexcept {
    const double a = kappa * t;
    if (std::fabs(a) < 1e-8) return t * (1.0 - 0.5 * a);
    return -std::expm1(-a) / kappa;
}
} // namespace

G2pp::G2pp(double kappa1, double kappa2, double sigma1, double sigma2, double rho, std::vector<double> scaleTimes, std::vector<double> scales)
    : k1_(kappa1), k2_(kappa2), s1_(sigma1), s2_(sigma2), rho_(rho), times_(std::move(scaleTimes)), scales_(std::move(scales)) {
    RL_REQUIRE(std::isfinite(k1_) && std::isfinite(k2_), "G2pp: mean reversions must be finite");
    RL_REQUIRE(std::isfinite(s1_) && s1_ > 0.0, "G2pp: sigma1 must be positive");
    RL_REQUIRE(std::isfinite(s2_) && s2_ >= 0.0, "G2pp: sigma2 must be non-negative");
    RL_REQUIRE(std::isfinite(rho_) && std::fabs(rho_) <= 1.0, "G2pp: rho must be in [-1, 1]");
    RL_REQUIRE(times_.size() == scales_.size(), "G2pp: scale times and values must match");
    for (std::size_t i = 0; i < times_.size(); ++i) {
        RL_REQUIRE(std::isfinite(times_[i]) && times_[i] > (i ? times_[i - 1] : 0.0), "G2pp: scale times must be positive and strictly increasing");
        RL_REQUIRE(std::isfinite(scales_[i]) && scales_[i] >= 0.0, "G2pp: scales must be non-negative");
    }
}

G2pp G2pp::withScale(std::vector<double> scaleTimes, std::vector<double> scales) const {
    return G2pp(k1_, k2_, s1_, s2_, rho_, std::move(scaleTimes), std::move(scales));
}

double G2pp::H1(double t) const noexcept { return H(k1_, t); }
double G2pp::H2(double t) const noexcept { return H(k2_, t); }

G2pp::Loading G2pp::loading(double T, double reference) const noexcept {
    auto diff = [&](double kappa) {
        return T >= reference ? std::exp(-kappa * reference) * H(kappa, T - reference) : -std::exp(-kappa * T) * H(kappa, reference - T);
    };
    return {diff(k1_), diff(k2_)};
}

G2pp::Covariance G2pp::covariance(double t) const {
    RL_REQUIRE(t >= 0.0, "G2pp::covariance: time must be non-negative");
    // Σ_j s_j² ∫ over the part of interval j below t. With no scale given there is one interval, s = 1.
    double i11 = 0.0, i12 = 0.0, i22 = 0.0, from = 0.0;
    const std::size_t n = times_.size();
    for (std::size_t j = 0; from < t; ++j) {
        const bool last = j + 1 >= n;
        const double to = last ? t : std::min(t, times_[j]);
        const double s = n ? scales_[std::min(j, n - 1)] : 1.0;
        const double w = s * s;
        i11 += w * (E(2.0 * k1_, to) - E(2.0 * k1_, from));
        i12 += w * (E(k1_ + k2_, to) - E(k1_ + k2_, from));
        i22 += w * (E(2.0 * k2_, to) - E(2.0 * k2_, from));
        from = to;
    }
    return {s1_ * s1_ * i11, rho_ * s1_ * s2_ * i12, s2_ * s2_ * i22};
}

double G2pp::bondVariance(double T, const Covariance& c) const noexcept {
    const double h1 = H1(T), h2 = H2(T);
    return h1 * h1 * c.c11 + 2.0 * h1 * h2 * c.c12 + h2 * h2 * c.c22;
}

double G2pp::deflatedBond(double dfT, double T, const Covariance& c, double x1, double x2) const noexcept {
    return dfT * std::exp(-H1(T) * x1 - H2(T) * x2 - 0.5 * bondVariance(T, c));
}

double G2pp::forwardRateVolatility(double tau) const noexcept {
    // df(t, t + τ) = σ_1 e^{−κ_1 τ} dW_1 + σ_2 e^{−κ_2 τ} dW_2 (+ drift).
    const double a = s1_ * std::exp(-k1_ * tau), b = s2_ * std::exp(-k2_ * tau);
    return std::sqrt(a * a + 2.0 * rho_ * a * b + b * b);
}

double G2pp::forwardRateCorrelation(double tau1, double tau2) const noexcept {
    const double a1 = s1_ * std::exp(-k1_ * tau1), b1 = s2_ * std::exp(-k2_ * tau1);
    const double a2 = s1_ * std::exp(-k1_ * tau2), b2 = s2_ * std::exp(-k2_ * tau2);
    const double cov = a1 * a2 + rho_ * (a1 * b2 + a2 * b1) + b1 * b2;
    const double v = forwardRateVolatility(tau1) * forwardRateVolatility(tau2);
    return v > 0.0 ? std::clamp(cov / v, -1.0, 1.0) : 1.0;
}

} // namespace rl::models
