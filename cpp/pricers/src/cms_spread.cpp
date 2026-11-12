#include "rl/pricers/cms_spread.hpp"

#include "rl/core/error.hpp"
#include "rl/math/normal.hpp"

#include <cmath>
#include <vector>

namespace rl::pricers {

namespace {

using models::G2pp;

/// A swap's par rate as a function of the state: strike ± Σ c_k e^{−L_k·x} / Σ a_k e^{−L_k·x}.
struct RateFunction {
    struct Term {
        double h1, h2, c;
    };
    std::vector<Term> value, annuity;
    double strike = 0.0, sign = 1.0;

    RateFunction(const G2pp& model, const G2pp::Covariance& cov, const models::ExerciseRight& r) : strike(r.strike), sign(r.payer ? 1.0 : -1.0) {
        auto terms = [&](const std::vector<models::BondFlow>& flows) {
            std::vector<Term> out;
            for (const models::BondFlow& f : flows) {
                const G2pp::Loading l = model.loading(f.time, r.time);
                out.push_back({l.h1, l.h2, f.amount * f.df * std::exp(-0.5 * G2pp::quadraticForm(l, cov))});
            }
            return out;
        };
        value = terms(r.flows);
        annuity = terms(r.annuity);
    }

    double operator()(double x1, double x2) const noexcept {
        double v = 0.0, a = 0.0;
        for (const Term& t : value) v += t.c * std::exp(-t.h1 * x1 - t.h2 * x2);
        for (const Term& t : annuity) a += t.c * std::exp(-t.h1 * x1 - t.h2 * x2);
        return strike + sign * v / a;
    }
};

} // namespace

CmsSpreadResult priceCmsSpread(const G2pp& model, const models::ExerciseRight& longSwap, const models::ExerciseRight& shortSwap, double strike, bool call,
                               double dfExpiry, int points) {
    RL_REQUIRE(std::fabs(longSwap.time - shortSwap.time) < 1e-9, "priceCmsSpread: the two swaps must fix at the same time");
    RL_REQUIRE(points >= 3 && points % 2 == 1, "priceCmsSpread: need an odd number of points per axis");
    RL_REQUIRE(dfExpiry > 0.0, "priceCmsSpread: the discount factor must be positive");
    const double t = longSwap.time;
    const G2pp::Covariance cov = model.covariance(t);
    RL_REQUIRE(cov.c11 > 0.0, "priceCmsSpread: the state has no variance at expiry");
    const RateFunction longRate(model, cov, longSwap), shortRate(model, cov, shortSwap);

    // Independent standard normals (u, v): x_1 = sd1 u, x_2 = β x_1 + sd2 v. With no independent
    // second factor the v direction carries nothing and the integral is one-dimensional.
    const double sd1 = std::sqrt(cov.c11), beta = cov.c12 / cov.c11;
    const double cond = cov.c22 - cov.c12 * beta;
    const double sd2 = cond > 1e-14 * std::max(cov.c22, 1e-300) ? std::sqrt(cond) : 0.0;
    const double lo = -7.5, step = 15.0 / (points - 1);
    auto weight = [&](int i) { return ((i == 0 || i == points - 1) ? 1.0 : (i % 2 ? 4.0 : 2.0)) * step / 3.0; };
    const int inner = sd2 > 0.0 ? points : 1;

    double payoff = 0.0, eL = 0.0, eS = 0.0, eLL = 0.0, eSS = 0.0, eLS = 0.0, mass = 0.0;
    for (int i = 0; i < points; ++i) {
        const double u = lo + i * step, x1 = sd1 * u, wu = weight(i) * math::normalPdf(u);
        for (int j = 0; j < inner; ++j) {
            const double v = inner > 1 ? lo + j * step : 0.0;
            const double w = wu * (inner > 1 ? weight(j) * math::normalPdf(v) : 1.0);
            const double a = longRate(x1, beta * x1 + sd2 * v), b = shortRate(x1, beta * x1 + sd2 * v);
            const double spread = a - b - strike;
            payoff += w * std::max(call ? spread : -spread, 0.0);
            mass += w;
            eL += w * a;
            eS += w * b;
            eLL += w * a * a;
            eSS += w * b * b;
            eLS += w * a * b;
        }
    }
    // The grid's own total weight differs from 1 by the truncated tails; normalising removes it from the moments.
    eL /= mass;
    eS /= mass;
    const double varL = std::max(eLL / mass - eL * eL, 0.0), varS = std::max(eSS / mass - eS * eS, 0.0), covLS = eLS / mass - eL * eS;
    CmsSpreadResult r;
    r.price = dfExpiry * payoff;
    r.longRate = eL;
    r.shortRate = eS;
    r.forwardSpread = eL - eS;
    r.longVolBp = std::sqrt(varL / t) * 1e4;
    r.shortVolBp = std::sqrt(varS / t) * 1e4;
    r.spreadVolBp = std::sqrt(std::max(varL + varS - 2.0 * covLS, 0.0) / t) * 1e4;
    r.rateCorrelation = varL > 0.0 && varS > 0.0 ? covLS / std::sqrt(varL * varS) : 1.0;
    return r;
}

} // namespace rl::pricers
