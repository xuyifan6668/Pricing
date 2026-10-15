#pragma once
// definitions of the Hagan SABR expansion (formulas and sources in sabr.hpp).
#include "rl/core/error.hpp"
#include "rl/math/dual.hpp"
#include "rl/math/rootfind.hpp"
#include "rl/vol/sabr.hpp"

#include <algorithm>
#include <cmath>
#include <string>

namespace rl::vol {

namespace detail {

/// exp(x) − 1 without cancellation: Taylor series to x⁷ for |x| < 1e-2 (truncation < x⁸/8! ≈ 2.5e-21).
template <math::Scalar T>
T expm1(const T& x) {
    using std::exp;
    if (std::fabs(math::value(x)) < 1e-2)
        return x * (1.0 + x * (1.0 / 2.0 + x * (1.0 / 6.0 + x * (1.0 / 24.0 + x * (1.0 / 120.0 + x * (1.0 / 720.0 + x * (1.0 / 5040.0)))))));
    return exp(x) - 1.0;
}

/// log(1 + x) without cancellation: alternating series to x⁹ for |x| < 1e-2 (truncation < x¹⁰/10 = 1e-21).
template <math::Scalar T>
T log1p(const T& x) {
    using std::log;
    if (std::fabs(math::value(x)) < 1e-2)
        return x * (1.0 - x * (1.0 / 2.0 - x * (1.0 / 3.0 - x * (1.0 / 4.0 - x * (1.0 / 5.0 - x * (1.0 / 6.0 - x * (1.0 / 7.0 - x * (1.0 / 8.0 - x / 9.0))))))));
    return log(1.0 + x);
}

/// ζ / x̂(ζ), x̂(ζ) = ln((√(1 − 2ρζ + ζ²) + ζ − ρ)/(1 − ρ)).
/// Written with log1p on an argument that has no cancellation: √D − 1 = ζ(ζ − 2ρ)/(√D + 1), so
///   ζ ≥ 0:  x̂ =  log1p( ζ (ζ − 2ρ + √D + 1) / ((√D + 1)(1 − ρ)) )
///   ζ < 0:  x̂ = −log1p( ζ (ζ − 2ρ − √D − 1) / ((√D + 1)(1 + ρ)) )    [from x̂ = ln((1+ρ)/(√D + ρ − ζ))]
/// and the series 1 − ρζ/2 + (2 − 3ρ²)ζ²/12 for |ζ| < 1e-6 (truncation O(ζ³) ≈ 1e-18), where x̂ → 0.
template <math::Scalar T>
T zetaOverXhat(const T& zeta, const T& rho) {
    using std::sqrt;
    if (std::fabs(math::value(zeta)) < 1e-6) return 1.0 - 0.5 * rho * zeta + (2.0 - 3.0 * rho * rho) * zeta * zeta / 12.0;
    const T sd = sqrt(1.0 - 2.0 * rho * zeta + zeta * zeta);
    const T xhat = math::value(zeta) >= 0.0 ? log1p(T(zeta * (zeta - 2.0 * rho + sd + 1.0) / ((sd + 1.0) * (1.0 - rho))))
                                            : T(-log1p(T(zeta * (zeta - 2.0 * rho - sd - 1.0) / ((sd + 1.0) * (1.0 + rho)))));
    return zeta / xhat;
}

inline void requireSabrQuery(double expiry, double strike) {
    RL_REQUIRE(std::isfinite(expiry) && expiry >= 0.0, "SABR: expiry must be non-negative");
    RL_REQUIRE(std::isfinite(strike), "SABR: strike must be finite");
}

} // namespace detail

template <math::Scalar T>
T sabrNormalVol(const market::SabrParams<T>& p, const T& forward, double strike, double expiry) {
    using std::exp;
    using std::log;
    market::validateSabrParams(p);
    detail::requireSabrQuery(expiry, strike);
    const double beta = p.beta;
    const T& alpha = p.alpha;
    const T& rho = p.rho;
    const T& nu = p.nu;

    if (beta == 0.0) {
        // Normal SABR: first factor is α exactly, ζ = (ν/α)(F − K); the shift cancels.
        const T zeta = (nu / alpha) * (forward - strike);
        const T corr = (2.0 - 3.0 * rho * rho) * nu * nu / 24.0;
        return alpha * detail::zetaOverXhat(zeta, rho) * (1.0 + corr * expiry);
    }

    const T f = forward + p.shift;
    const double k = strike + p.shift;
    RL_REQUIRE(math::value(f) > 0.0 && k > 0.0, "sabrNormalVol: forward + shift and strike + shift must be positive when beta > 0 (F+s = " +
                                                    std::to_string(math::value(f)) + ", K+s = " + std::to_string(k) + ")");
    const double c = 1.0 - beta;
    const T r = log(f / k); // log-moneyness; f → k is r → 0
    const T em1 = detail::expm1(r);
    const bool nearAtm = std::fabs(math::value(r)) < 1e-8;
    T factor0, q;
    if (beta == 1.0) {
        // (f − k)/ln(f/k) = k·expm1(r)/r
        const T ratio = nearAtm ? T(1.0 + r * (0.5 + r / 6.0)) : T(em1 / r);
        factor0 = alpha * k * ratio;
        q = r;
    } else {
        // (1−β)(f − k)/(f^{1−β} − k^{1−β}) = k^β · c·expm1(r)/expm1(c r),  c = 1 − β
        const T emc = detail::expm1(T(c * r));
        const T ratio = nearAtm ? T((1.0 + r * (0.5 + r / 6.0)) / (1.0 + c * r * (0.5 + c * r / 6.0))) : T(c * em1 / emc);
        factor0 = alpha * std::pow(k, beta) * ratio;
        q = std::pow(k, c) * emc / c; // (f^{1−β} − k^{1−β})/(1−β) without cancellation
    }
    const T zeta = (nu / alpha) * q;
    const T favC = exp(c * 0.5 * (log(f) + std::log(k))); // f_av^{1−β}
    const T corr = -beta * (2.0 - beta) * alpha * alpha / (24.0 * favC * favC) + rho * alpha * nu * beta / (4.0 * favC) + (2.0 - 3.0 * rho * rho) * nu * nu / 24.0;
    return factor0 * detail::zetaOverXhat(zeta, rho) * (1.0 + corr * expiry);
}

template <math::Scalar T>
T sabrLognormalVol(const market::SabrParams<T>& p, const T& forward, double strike, double expiry) {
    using std::exp;
    using std::log;
    market::validateSabrParams(p);
    detail::requireSabrQuery(expiry, strike);
    const T f = forward + p.shift;
    const double k = strike + p.shift;
    RL_REQUIRE(math::value(f) > 0.0 && k > 0.0, "sabrLognormalVol: forward + shift and strike + shift must be positive");
    const double beta = p.beta, c = 1.0 - beta;
    const T& alpha = p.alpha;
    const T& rho = p.rho;
    const T& nu = p.nu;
    const T r = log(f / k);
    const T fkC2 = exp(0.5 * c * (log(f) + std::log(k))); // (fk)^{(1−β)/2}
    const T zeta = (nu / alpha) * fkC2 * r;
    const T r2 = r * r;
    const T denom = fkC2 * (1.0 + c * c * r2 / 24.0 + c * c * c * c * r2 * r2 / 1920.0);
    const T corr = c * c * alpha * alpha / (24.0 * fkC2 * fkC2) + 0.25 * rho * beta * nu * alpha / fkC2 + (2.0 - 3.0 * rho * rho) * nu * nu / 24.0;
    return alpha / denom * detail::zetaOverXhat(zeta, rho) * (1.0 + corr * expiry);
}

template <math::Scalar T>
T sabrAtmNormalVol(const market::SabrParams<T>& p, const T& forward, double expiry) {
    using std::pow;
    market::validateSabrParams(p);
    detail::requireSabrQuery(expiry, 0.0);
    const T& alpha = p.alpha;
    const T& rho = p.rho;
    const T& nu = p.nu;
    const double beta = p.beta;
    if (beta == 0.0) return alpha * (1.0 + (2.0 - 3.0 * rho * rho) * nu * nu * expiry / 24.0);
    const T f = forward + p.shift;
    RL_REQUIRE(math::value(f) > 0.0, "sabrAtmNormalVol: forward + shift must be positive when beta > 0");
    const T fb = pow(f, beta);
    const T fc = pow(f, 1.0 - beta);
    const T corr = -beta * (2.0 - beta) * alpha * alpha / (24.0 * fc * fc) + rho * alpha * nu * beta / (4.0 * fc) + (2.0 - 3.0 * rho * rho) * nu * nu / 24.0;
    return alpha * fb * (1.0 + corr * expiry);
}

template <math::Scalar T>
T alphaFromAtmNormalVol(const T& atmNormalVol, const T& forward, double expiry, double beta, const T& rho, const T& nu, double shift) {
    using std::pow;
    RL_REQUIRE(math::isFinite(atmNormalVol) && math::value(atmNormalVol) > 0.0, "alphaFromAtmNormalVol: ATM normal vol must be positive");
    RL_REQUIRE(std::isfinite(expiry) && expiry >= 0.0, "alphaFromAtmNormalVol: expiry must be non-negative");
    RL_REQUIRE(std::isfinite(beta) && beta >= 0.0 && beta <= 1.0, "alphaFromAtmNormalVol: beta must lie in [0, 1]");
    RL_REQUIRE(math::value(rho) > -1.0 && math::value(rho) < 1.0, "alphaFromAtmNormalVol: rho must lie in (-1, 1)");
    RL_REQUIRE(math::value(nu) >= 0.0, "alphaFromAtmNormalVol: nu must be non-negative");

    // σ_ATM(α) = a1 α + a2 α² + a3 α³  (sabrAtmNormalVol expanded in α); a3 ≤ 0 < a1.
    T a1, a2, a3;
    const T lead = 1.0 + (2.0 - 3.0 * rho * rho) * nu * nu * expiry / 24.0;
    if (beta == 0.0) {
        a1 = lead;
        a2 = T(0.0);
        a3 = T(0.0);
    } else {
        const T f = forward + shift;
        RL_REQUIRE(math::value(f) > 0.0, "alphaFromAtmNormalVol: forward + shift must be positive when beta > 0");
        const T fb = pow(f, beta);
        const T fc = pow(f, 1.0 - beta);
        a1 = fb * lead;
        a2 = fb * rho * nu * beta * expiry / (4.0 * fc);
        a3 = -fb * beta * (2.0 - beta) * expiry / (24.0 * fc * fc);
    }
    const double A1 = math::value(a1), A2 = math::value(a2), A3 = math::value(a3), S = math::value(atmNormalVol);
    if (!(A1 > 0.0))
        throw NumericalError("alphaFromAtmNormalVol: 1 + (2 - 3 rho^2) nu^2 T / 24 <= 0; the expansion is meaningless for these rho, nu, T");

    double a0 = S / A1;
    if (A3 < 0.0) {
        // First positive stationary point of the cubic, written without cancellation for either sign of A2.
        const double disc = std::sqrt(4.0 * A2 * A2 - 12.0 * A3 * A1);
        const double aStar = A2 < 0.0 ? 2.0 * A1 / (-2.0 * A2 + disc) : (-2.0 * A2 - disc) / (6.0 * A3);
        auto g = [&](double a) { return ((A3 * a + A2) * a + A1) * a - S; };
        const double gMax = g(aStar);
        if (gMax < 0.0)
            throw NumericalError("alphaFromAtmNormalVol: ATM normal vol " + std::to_string(S * 1e4) + " bp exceeds the maximum " +
                                 std::to_string((gMax + S) * 1e4) + " bp the SABR expansion can produce for these (beta, rho, nu, T)");
        math::RootOptions ro;
        ro.xTol = 1e-15 * aStar;
        ro.maxIter = 200;
        auto fdf = [&](double a) { return std::pair<double, double>{g(a), (3.0 * A3 * a + 2.0 * A2) * a + A1}; };
        a0 = math::newtonSafeguarded(fdf, std::clamp(a0, 0.0, aStar), 0.0, aStar, ro).x;
    }
    // D20: one Newton step in T arithmetic at the double root (value unchanged; tangent = implicit-function derivative).
    const T g = ((a3 * a0 + a2) * a0 + a1) * a0 - atmNormalVol;
    const T dg = (3.0 * a3 * a0 + 2.0 * a2) * a0 + a1;
    return a0 - g / dg;
}

template <math::Scalar T>
SabrVolGradient<T> sabrNormalVolGradient(const market::SabrParams<T>& p, const T& forward, double strike, double expiry) {
    using D = math::Dual<T, 3>;
    market::SabrParams<D> pd;
    pd.alpha = D::variable(p.alpha, 0);
    pd.rho = D::variable(p.rho, 1);
    pd.nu = D::variable(p.nu, 2);
    pd.beta = p.beta;
    pd.shift = p.shift;
    const D v = sabrNormalVol(pd, D::constant(forward), strike, expiry);
    return {v.base(), v.tangent(0), v.tangent(1), v.tangent(2)};
}

template <math::Scalar T>
SabrAtmConstrainedGradient<T> sabrNormalVolGradientAtmConstrained(const T& atmNormalVol, const T& forward, double strike, double expiry, double beta,
                                                                  const T& rho, const T& nu, double shift) {
    using D = math::Dual<T, 2>;
    const D rd = D::variable(rho, 0);
    const D nd = D::variable(nu, 1);
    const D fd = D::constant(forward);
    const D a = alphaFromAtmNormalVol(D::constant(atmNormalVol), fd, expiry, beta, rd, nd, shift);
    market::SabrParams<D> pd;
    pd.alpha = a;
    pd.rho = rd;
    pd.nu = nd;
    pd.beta = beta;
    pd.shift = shift;
    const D v = sabrNormalVol(pd, fd, strike, expiry);
    return {v.base(), v.tangent(0), v.tangent(1), a.base(), a.tangent(0), a.tangent(1)};
}

} // namespace rl::vol

#define RL_VOL_INSTANTIATE_SABR(T)                                                                                                             \
    namespace rl::vol {                                                                                                                        \
    template T sabrNormalVol<T>(const market::SabrParams<T>&, const T&, double, double);                                                       \
    template T sabrLognormalVol<T>(const market::SabrParams<T>&, const T&, double, double);                                                    \
    template T sabrAtmNormalVol<T>(const market::SabrParams<T>&, const T&, double);                                                            \
    template T alphaFromAtmNormalVol<T>(const T&, const T&, double, double, const T&, const T&, double);                                       \
    template SabrVolGradient<T> sabrNormalVolGradient<T>(const market::SabrParams<T>&, const T&, double, double);                              \
    template SabrAtmConstrainedGradient<T> sabrNormalVolGradientAtmConstrained<T>(const T&, const T&, double, double, double, const T&,        \
                                                                                  const T&, double);                                           \
    }
