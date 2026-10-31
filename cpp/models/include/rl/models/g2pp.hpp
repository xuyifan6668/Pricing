#pragma once
// The two-factor Gaussian model G2++ (Brigo & Mercurio, "Interest Rate Models", ch. 4), written in
// the same deflated-bond form as the one-factor LGM in lgm.hpp so that the two share every pricer
// idea. It is also the two-factor Hull–White model; the three are the same model.
//
// Why a second factor. In a one-factor model every rate is a function of one number, so all rates
// move together: the 2-year and the 30-year rate are perfectly correlated over any short interval.
// Two consequences:
//   * one mean reversion cannot match swaption vols across BOTH expiry and tenor: it fits a strip
//     and misses the rest of the matrix;
//   * anything whose value depends on rates of different maturities moving apart, or on how
//     correlated the swap rates seen on different exercise dates are, is priced with a correlation
//     the model was never asked about.
// With two factors of different speeds the curve can steepen and flatten as well as shift.
//
// Two states with x(0) = 0 and, under the model's numeraire measure, no drift:
//
//     dx_i = s(t) σ_i e^{κ_i t} dW_i,     d<W_1, W_2> = ρ dt,     x(t) ~ N(0, C(t)),
//     C_ij(t) = ρ_ij σ_i σ_j ∫₀ᵗ s(u)² e^{(κ_i + κ_j) u} du.
//
// The deflated zero-coupon bond is one exponential, exactly as in LGM:
//
//     P(t, T; x) / N(t, x) = D(T) · exp( −H(T)·x − ½ H(T)ᵀ C(t) H(T) ),    H_i(T) = (1 − e^{−κ_i T})/κ_i.
//
// So the model fits today's curve by construction, deflated values are martingales in a driftless
// Gaussian state, and Monte Carlo steps between any two dates are exact Gaussian draws.
//
// Only differences of H matter: replacing H by H − H(R) for any fixed R is the same model under the
// numeraire P(·, R) (times a constant). Pricers use that freedom. With a fast factor C(t) grows like
// e^{2κt}, and under the textbook numeraire above deflated values become products of enormous and
// tiny numbers; measured from the exercise date (European) or the last cashflow (Monte Carlo) every
// quantity stays the size of a rate volatility. `loading` returns H(T) − H(R).
//
// In Brigo–Mercurio's notation r = x + y + φ(t) with dx = −a x dt + σ dW_1, dy = −b y dt + η dW_2:
// a = κ_1, b = κ_2, σ = σ_1, η = σ_2, ρ = ρ. With s ≡ 1 this class IS that model (tested against
// their zero-bond option formula). With σ_2 = 0 it is Lgm1F with constant Hull–White volatility.
//
// s(t) is a piecewise-constant scale common to both factors. It changes the LEVEL of volatility
// through time and leaves the SHAPE (the two speeds, the ratio σ_2/σ_1, ρ) alone, so the
// correlation structure can be chosen first and the level then bootstrapped to a swaption strip,
// one expiry at a time, as ζ is in the one-factor model (g2pp_calibration.hpp).
//
// Immutable after construction; safe to share across threads.

#include <cstddef>
#include <vector>

namespace rl::models {

class G2pp {
public:
    /// Covariance matrix of the state (x_1, x_2) at some time.
    struct Covariance {
        double c11 = 0.0, c12 = 0.0, c22 = 0.0;
    };

    /// Constant parameters (s ≡ 1) when no scale is given. Otherwise scale s_j applies on
    /// (t_{j−1}, t_j], t_{−1} = 0, and the last scale continues beyond the last time.
    /// Throws InputError unless κ_i are finite, σ_1 > 0, σ_2 ≥ 0, |ρ| ≤ 1, the scale times are
    /// positive and strictly increasing, and the scales are non-negative.
    G2pp(double kappa1, double kappa2, double sigma1, double sigma2, double rho, std::vector<double> scaleTimes = {}, std::vector<double> scales = {});

    double kappa1() const noexcept { return k1_; }
    double kappa2() const noexcept { return k2_; }
    double sigma1() const noexcept { return s1_; }
    double sigma2() const noexcept { return s2_; }
    double rho() const noexcept { return rho_; }
    const std::vector<double>& scaleTimes() const noexcept { return times_; }
    const std::vector<double>& scales() const noexcept { return scales_; }

    /// The same shape with a different volatility scale.
    G2pp withScale(std::vector<double> scaleTimes, std::vector<double> scales) const;

    double H1(double t) const noexcept;
    double H2(double t) const noexcept;

    /// H(T) − H(reference), computed without cancellation as ±e^{−κ·min} · H(|T − reference|).
    struct Loading {
        double h1 = 0.0, h2 = 0.0;
    };
    Loading loading(double T, double reference) const noexcept;

    /// Lᵀ C L.
    static double quadraticForm(const Loading& l, const Covariance& c) noexcept { return l.h1 * l.h1 * c.c11 + 2.0 * l.h1 * l.h2 * c.c12 + l.h2 * l.h2 * c.c22; }

    /// C(t). Requires t ≥ 0.
    Covariance covariance(double t) const;

    /// H(T)ᵀ C H(T): the variance of the log deflated bond maturing at T.
    double bondVariance(double T, const Covariance& c) const noexcept;

    /// P(t, T; x) / N(t, x) given D(T) and C(t).
    double deflatedBond(double dfT, double T, const Covariance& c, double x1, double x2) const noexcept;

    /// Instantaneous correlation between the forward rates with times to maturity tau1 and tau2.
    /// It is exactly 1 in a one-factor model; here it falls as the two maturities move apart.
    double forwardRateCorrelation(double tau1, double tau2) const noexcept;

    /// Instantaneous volatility of the forward rate with time to maturity tau, at scale 1.
    double forwardRateVolatility(double tau) const noexcept;

private:
    double k1_, k2_, s1_, s2_, rho_;
    std::vector<double> times_, scales_;
};

} // namespace rl::models
