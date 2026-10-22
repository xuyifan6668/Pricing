#include "rl/models/lgm_swaption.hpp"

#include "rl/core/error.hpp"
#include "rl/math/normal.hpp"
#include "rl/math/rootfind.hpp"
#include "rl/vol/bachelier.hpp"

#include <cmath>

namespace rl::models {

double deflatedValue(const Lgm1F& model, double zeta, const std::vector<BondFlow>& flows, double x) noexcept {
    double g = 0.0;
    for (const BondFlow& f : flows) g += f.amount * model.deflatedBond(f.df, f.time, zeta, x);
    return g;
}

ExerciseBoundary exerciseBoundary(const Lgm1F& model, double zeta, const std::vector<BondFlow>& flows) {
    RL_REQUIRE(zeta > 0.0, "exerciseBoundary: zeta must be positive");
    const double sd = std::sqrt(zeta);
    auto g = [&](double x) { return deflatedValue(model, zeta, flows, x); };

    // Scan ±12 standard deviations for sign changes. Beyond that the density is below 1e-32 and
    // whatever the payoff does there cannot matter to a price.
    const int n = 480;
    const double lo = -12.0 * sd, step = 24.0 * sd / n;
    ExerciseBoundary b;
    int crossings = 0;
    double a = lo, ga = g(a), bracketLo = 0.0, bracketHi = 0.0;
    for (int i = 1; i <= n; ++i) {
        const double c = lo + i * step, gc = g(c);
        if ((ga < 0.0) != (gc < 0.0)) {
            ++crossings;
            bracketLo = a;
            bracketHi = c;
            b.exerciseAbove = gc > 0.0;
        }
        a = c;
        ga = gc;
    }
    if (crossings == 0) {
        b.exerciseAbove = ga > 0.0;
        return b;
    }
    if (crossings > 1) throw NumericalError("exerciseBoundary: the payoff changes sign more than once; it is not a swap-like payoff");
    math::RootOptions ro;
    ro.xTol = 1e-14 * sd;
    ro.maxIter = 200;
    b.exists = true;
    b.x = math::brent(g, bracketLo, bracketHi, ro).x;
    return b;
}

double lgmOptionValue(const Lgm1F& model, double zeta, const std::vector<BondFlow>& flows) {
    const ExerciseBoundary b = exerciseBoundary(model, zeta, flows);
    if (!b.exists) {
        double v = 0.0;
        for (const BondFlow& f : flows) v += f.amount * f.df;
        return b.exerciseAbove ? v : 0.0;
    }
    // ∫_{x > x*} exp(−H x − ½H²ζ) φ(x; 0, ζ) dx = Φ(−(x* + Hζ)/√ζ): completing the square shifts the
    // Gaussian's mean by −Hζ. Every bond shares the one boundary x*.
    const double sd = std::sqrt(zeta);
    const double side = b.exerciseAbove ? -1.0 : 1.0;
    double v = 0.0;
    for (const BondFlow& f : flows) v += f.amount * f.df * math::normalCdf(side * (b.x + model.H(f.time) * zeta) / sd);
    return v;
}

double lgmOptionValueByIntegration(const Lgm1F& model, double zeta, const std::vector<BondFlow>& flows, int points) {
    RL_REQUIRE(zeta > 0.0 && points >= 3 && points % 2 == 1, "lgmOptionValueByIntegration: need zeta > 0 and an odd point count");
    const double sd = std::sqrt(zeta);
    const double lo = -10.0 * sd, h = 20.0 * sd / (points - 1);
    double sum = 0.0;
    for (int i = 0; i < points; ++i) {
        const double x = lo + i * h;
        const double g = deflatedValue(model, zeta, flows, x);
        const double w = (i == 0 || i == points - 1) ? 1.0 : (i % 2 ? 4.0 : 2.0);
        if (g > 0.0) sum += w * g * math::normalPdf(x / sd) / sd;
    }
    return sum * h / 3.0;
}

double lgmSwaptionPrice(const Lgm1F& model, const ExerciseRight& right) { return lgmOptionValue(model, model.zeta(right.time), right.flows); }

double marketSwaptionPrice(const ExerciseRight& right, double normalVol) {
    const vol::OptionType type = right.payer ? vol::OptionType::Call : vol::OptionType::Put;
    return right.annuityValue() * vol::bachelierPrice<double>(right.forwardRate(), right.strike, right.time, normalVol, type);
}

double impliedNormalVol(const ExerciseRight& right, double price) {
    const vol::OptionType type = right.payer ? vol::OptionType::Call : vol::OptionType::Put;
    return vol::bachelierImpliedVol(price / right.annuityValue(), right.forwardRate(), right.strike, right.time, type);
}

} // namespace rl::models
