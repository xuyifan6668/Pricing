#include "rl/models/g2pp_swaption.hpp"

#include "rl/core/error.hpp"
#include "rl/math/normal.hpp"
#include "rl/math/quadrature.hpp"
#include "rl/math/rootfind.hpp"

#include <algorithm>
#include <cmath>
#include <map>

namespace rl::models {

namespace {

struct Term {
    double h1, h2, c; ///< loadings L_k on the two states and the coefficient a_k D_k exp(−½ L_kᵀ C L_k)
};

/// The payoff under the numeraire P(·, t): loadings H(T_k) − H(t), so that c_k is a_k D_k times the
/// convexity of a bond option's variance and nothing in the sum is large.
std::vector<Term> terms(const G2pp& model, double t, const G2pp::Covariance& cov, const std::vector<BondFlow>& flows) {
    std::vector<Term> out;
    out.reserve(flows.size());
    for (const BondFlow& f : flows) {
        const G2pp::Loading l = model.loading(f.time, t);
        out.push_back({l.h1, l.h2, f.amount * f.df * std::exp(-0.5 * G2pp::quadraticForm(l, cov))});
    }
    return out;
}

/// E[(Σ b_k e^{−h_k y})⁺] for y ~ N(mean, sd²), the payoff changing sign at most once.
double singleCrossingExpectation(const std::vector<double>& b, const std::vector<double>& h, double mean, double sd) {
    auto f = [&](double y) {
        double v = 0.0;
        for (std::size_t k = 0; k < b.size(); ++k) v += b[k] * std::exp(-h[k] * y);
        return v;
    };
    auto whole = [&] {
        double v = 0.0;
        for (std::size_t k = 0; k < b.size(); ++k) v += b[k] * std::exp(-h[k] * mean + 0.5 * h[k] * h[k] * sd * sd);
        return v;
    };
    if (sd <= 0.0) return std::max(f(mean), 0.0);
    // Beyond 12 standard deviations the density is below 1e-32.
    const double lo = mean - 12.0 * sd, hi = mean + 12.0 * sd;
    const double flo = f(lo), fhi = f(hi);
    if ((flo < 0.0) == (fhi < 0.0)) return fhi > 0.0 ? whole() : 0.0;
    math::RootOptions ro;
    ro.xTol = 1e-13 * sd;
    ro.maxIter = 200;
    const double y = math::brent(f, lo, hi, ro).x;
    // ∫_{y > y*} e^{−h y} φ(y; m, s²) dy = e^{−h m + ½ h² s²} Φ(−(y* − m + h s²)/s); the other side flips the sign.
    const double side = fhi > 0.0 ? -1.0 : 1.0;
    double v = 0.0;
    for (std::size_t k = 0; k < b.size(); ++k)
        v += b[k] * std::exp(-h[k] * mean + 0.5 * h[k] * h[k] * sd * sd) * math::normalCdf(side * (y - mean + h[k] * sd * sd) / sd);
    return v;
}

const math::GaussHermite& rule(int n) {
    thread_local std::map<int, math::GaussHermite> cache;
    auto it = cache.find(n);
    if (it == cache.end()) it = cache.emplace(n, math::GaussHermite(static_cast<std::size_t>(n))).first;
    return it->second;
}

} // namespace

double g2ppDeflatedValue(const G2pp& model, const G2pp::Covariance& c, const std::vector<BondFlow>& flows, double x1, double x2) noexcept {
    double g = 0.0;
    for (const BondFlow& f : flows) g += f.amount * model.deflatedBond(f.df, f.time, c, x1, x2);
    return g;
}

double g2ppOptionValue(const G2pp& model, double t, const std::vector<BondFlow>& flows, int quadraturePoints) {
    RL_REQUIRE(t > 0.0, "g2ppOptionValue: the exercise time must be positive");
    RL_REQUIRE(quadraturePoints >= 2, "g2ppOptionValue: need at least two quadrature points");
    RL_REQUIRE(!flows.empty(), "g2ppOptionValue: no flows");
    const G2pp::Covariance cov = model.covariance(t);
    RL_REQUIRE(cov.c11 > 0.0, "g2ppOptionValue: the state has no variance at the exercise time");
    const std::vector<Term> ts = terms(model, t, cov, flows);
    const std::size_t n = ts.size();

    // Split the state into the variable the closed form is taken in, y = w·x, and what is independent
    // of it: x = g y + d ε with g = C w / (wᵀ C w), ε ~ N(0, 1) independent of y. Each bond is then
    // c_k exp(−(H_k·g) y − (H_k·d) ε).
    //
    // DESIGN: Brigo & Mercurio condition on x_1 (w = (0, 1)). That is exact, but the outer integrand
    // is only as smooth as the conditional variance of x_2 makes it: when the factors are close to
    // perfectly correlated it tends to the kinked payoff itself and Gauss–Hermite stops converging.
    // Taking y along the loadings that decide exercise puts almost all of the payoff's variance in
    // the closed-form direction, so the outer integrand is nearly flat, and exactly constant in the
    // one-factor limit. The closed form needs the payoff to cross zero once in y, which holds when the
    // loadings H_k·g are monotone in maturity; where they are not (strongly negative correlation can
    // do that), fall back to conditioning on one factor, whose loadings always are.
    auto split = [&](double w1, double w2, std::vector<double>& h, std::vector<double>& e, double& sdY) {
        const double cw1 = cov.c11 * w1 + cov.c12 * w2, cw2 = cov.c12 * w1 + cov.c22 * w2;
        const double vy = w1 * cw1 + w2 * cw2;
        const double g1 = cw1 / vy, g2 = cw2 / vy;
        // Covariance of the remainder, C − (C w)(C w)ᵀ / vy, has rank one: d dᵀ.
        const double m11 = std::max(cov.c11 - cw1 * g1, 0.0), m22 = std::max(cov.c22 - cw2 * g2, 0.0), m12 = cov.c12 - cw1 * g2;
        double d1 = 0.0, d2 = 0.0;
        if (m11 >= m22 && m11 > 1e-18 * cov.c11) {
            d1 = std::sqrt(m11);
            d2 = m12 / d1;
        } else if (m22 > 1e-18 * std::max(cov.c22, 1e-300)) {
            d2 = std::sqrt(m22);
            d1 = m12 / d2;
        }
        h.resize(n);
        e.resize(n);
        for (std::size_t k = 0; k < n; ++k) {
            h[k] = ts[k].h1 * g1 + ts[k].h2 * g2;
            e[k] = ts[k].h1 * d1 + ts[k].h2 * d2;
        }
        sdY = std::sqrt(vy);
        return d1 != 0.0 || d2 != 0.0;
    };
    auto monotone = [&](const std::vector<double>& h) {
        // Flows are in time order; equal times give equal loadings.
        bool up = true, down = true;
        for (std::size_t k = 1; k < n; ++k) {
            if (h[k] < h[k - 1]) up = false;
            if (h[k] > h[k - 1]) down = false;
        }
        return up || down;
    };

    std::vector<double> h, e, b(n);
    double sdY = 0.0;
    // Whether the payoff is positive depends on the state only through (H_k − H_first)·x, so the
    // direction that carries the exercise decision is the spread between the last and first loadings.
    const Term& first = ts.front();
    const Term& last = ts.back();
    bool twoDimensional = split(last.h1 - first.h1, last.h2 - first.h2, h, e, sdY);
    if (!monotone(h)) {
        // Condition on the factor that matters less to the exercise decision, solve in the other.
        const double w1 = last.h1 - first.h1, w2 = last.h2 - first.h2;
        const bool solveInFirst = w1 * w1 * cov.c11 >= w2 * w2 * cov.c22;
        twoDimensional = split(solveInFirst ? 1.0 : 0.0, solveInFirst ? 0.0 : 1.0, h, e, sdY);
    }
    if (!twoDimensional) {
        for (std::size_t k = 0; k < n; ++k) b[k] = ts[k].c;
        return singleCrossingExpectation(b, h, 0.0, sdY);
    }
    return rule(quadraturePoints).integrate([&](double eps) {
        for (std::size_t k = 0; k < n; ++k) b[k] = ts[k].c * std::exp(-e[k] * eps);
        return singleCrossingExpectation(b, h, 0.0, sdY);
    });
}

double g2ppOptionValueByIntegration(const G2pp& model, double t, const std::vector<BondFlow>& flows, int points) {
    RL_REQUIRE(t > 0.0 && points >= 3 && points % 2 == 1, "g2ppOptionValueByIntegration: need t > 0 and an odd point count");
    const G2pp::Covariance cov = model.covariance(t);
    const double beta = cov.c12 / cov.c11, sd1 = std::sqrt(cov.c11), sd2 = std::sqrt(std::max(cov.c22 - cov.c12 * beta, 0.0));
    RL_REQUIRE(sd2 > 0.0, "g2ppOptionValueByIntegration: the model is one-factor; integrate in one dimension instead");
    const std::vector<Term> ts = terms(model, t, cov, flows);
    // Independent standard normals (u, v): x_1 = sd1 u, x_2 = β x_1 + sd2 v.
    const double lo = -8.0, step = 16.0 / (points - 1);
    auto weight = [&](int i) { return (i == 0 || i == points - 1) ? 1.0 : (i % 2 ? 4.0 : 2.0); };
    double sum = 0.0;
    for (int i = 0; i < points; ++i) {
        const double u = lo + i * step, x1 = sd1 * u;
        double inner = 0.0;
        for (int j = 0; j < points; ++j) {
            const double v = lo + j * step, x2 = beta * x1 + sd2 * v;
            double g = 0.0;
            for (const Term& k : ts) g += k.c * std::exp(-k.h1 * x1 - k.h2 * x2);
            if (g > 0.0) inner += weight(j) * g * math::normalPdf(v);
        }
        sum += weight(i) * inner * math::normalPdf(u);
    }
    return sum * step * step / 9.0;
}

double g2ppSwaptionPrice(const G2pp& model, const ExerciseRight& right) { return g2ppOptionValue(model, right.time, right.flows); }

double g2ppZeroBondPut(const G2pp& model, double expiry, double maturity, double dfExpiry, double dfMaturity, double strikeAmount) {
    RL_REQUIRE(expiry > 0.0 && maturity > expiry && strikeAmount > 0.0, "g2ppZeroBondPut: need 0 < expiry < maturity and a positive strike");
    const G2pp::Covariance c = model.covariance(expiry);
    const double d1h = model.H1(maturity) - model.H1(expiry), d2h = model.H2(maturity) - model.H2(expiry);
    const double v = std::sqrt(d1h * d1h * c.c11 + 2.0 * d1h * d2h * c.c12 + d2h * d2h * c.c22);
    const double d1 = std::log(dfExpiry / (strikeAmount * dfMaturity)) / v + 0.5 * v;
    return dfExpiry * math::normalCdf(d1) - strikeAmount * dfMaturity * math::normalCdf(d1 - v);
}

} // namespace rl::models
