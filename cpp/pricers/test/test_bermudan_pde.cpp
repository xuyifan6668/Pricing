// The PDE pricer, tested on the properties that define a Bermudan and on its numerical behaviour.
#include "rl/pricers/bermudan_pde.hpp"

#include "rl/models/lgm_swaption.hpp"
#include "test_support.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <algorithm>
#include <cmath>

using namespace rl::pricers;
using namespace rl::pricers::test;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {
PdeOptions fine() {
    PdeOptions o;
    o.spaceSteps = 1600;
    o.timeSteps = 1600;
    return o;
}
} // namespace

namespace {
// An independent pricer for cross-checking: backward induction with the EXACT Gaussian transition,
//   V_{i-1}(x) = ∫ max(V_i, h_i)(y) φ(y − x; ζ_i − ζ_{i−1}) dy,
// by the trapezoid rule on a fine grid. No time stepping at all, so it shares no discretisation
// error with the finite-difference scheme (it is O(grid²) per date, which is why it is not the pricer).
double priceByTransitionIntegration(const BermudanProblem& p, int steps) {
    const std::size_t n = p.size(), N = static_cast<std::size_t>(steps) + 1;
    const double L = 9.0 * std::sqrt(p.zeta(n - 1)), h = 2.0 * L / steps, norm = 1.0 / std::sqrt(2.0 * M_PI);
    std::vector<double> x(N), v(N, 0.0), w(N);
    for (std::size_t j = 0; j < N; ++j) x[j] = -L + h * static_cast<double>(j);
    for (std::size_t i = n; i-- > 0;) {
        for (std::size_t j = 0; j < N; ++j) v[j] = std::max(v[j], p.exerciseValue(i, x[j]));
        const double sd = std::sqrt(p.zeta(i) - (i ? p.zeta(i - 1) : 0.0));
        auto integrate = [&](double from) {
            double s = 0.0;
            for (std::size_t k = 0; k < N; ++k) {
                const double u = (x[k] - from) / sd;
                if (std::fabs(u) < 9.0) s += v[k] * std::exp(-0.5 * u * u);
            }
            return s * h * norm / sd;
        };
        if (i == 0) return integrate(0.0);
        for (std::size_t j = 0; j < N; ++j) w[j] = integrate(x[j]);
        v.swap(w);
    }
    return 0.0;
}
} // namespace

TEST_CASE("PDE agrees with exact-transition integration to 0.01 bp of notional", "[pricers][pde][crosscheck]") {
    for (bool payer : {true, false}) {
        const CalibratedBermudan b = bermudan10nc1(0.035, payer);
        PdeOptions o = fine();
        o.rannacher = false;
        CHECK_THAT(priceBermudanPde(b.problem, o).value, WithinAbs(priceByTransitionIntegration(b.problem, 1500), 1e-6));
    }
}

TEST_CASE("A Bermudan with one exercise date is the European: PDE matches the closed form", "[pricers][pde]") {
    for (double kappa : {0.0, 0.03, 0.10}) {
        for (bool payer : {true, false}) {
            const CalibratedBermudan b = bermudan10nc1(0.035, payer, kappa);
            const std::vector<double> european = b.problem.europeanValues();
            for (std::size_t i : {std::size_t{0}, std::size_t{4}, std::size_t{8}}) {
                const double pde = priceBermudanPde(subset(b.problem, {i}), fine()).value;
                CHECK_THAT(pde, WithinRel(european[i], 2e-5));
            }
        }
    }
}

TEST_CASE("Bermudan >= every co-terminal European, and <= their sum", "[pricers][pde]") {
    for (double strike : {0.025, 0.035, 0.045}) {
        for (bool payer : {true, false}) {
            const CalibratedBermudan b = bermudan10nc1(strike, payer);
            const double v = priceBermudanPde(b.problem).value;
            const std::vector<double> european = b.problem.europeanValues();
            double sum = 0.0;
            for (double e : european) {
                CHECK(v >= e - 1e-7);
                sum += e;
            }
            CHECK(v <= sum); // holding ALL the Europeans dominates holding the one-shot choice among them
            CHECK(v > *std::max_element(european.begin(), european.end())); // the right to wait is worth something
        }
    }
}

TEST_CASE("More exercise dates never reduce the value", "[pricers][pde]") {
    const CalibratedBermudan b = bermudan10nc1();
    double prev = 0.0;
    std::vector<std::size_t> keep;
    for (std::size_t i = b.problem.size(); i-- > 0;) { // add dates from the last backwards
        keep.insert(keep.begin(), i);
        const double v = priceBermudanPde(subset(b.problem, keep)).value;
        CHECK(v >= prev - 1e-9);
        prev = v;
    }
}

TEST_CASE("Moneyness limits: deep in the money -> the swap from the first call date; deep out -> zero", "[pricers][pde]") {
    // Deep ITM payer (strike 0.5%): exercise at the first date is certain, value = forward value of that swap.
    const CalibratedBermudan itm = bermudan10nc1(0.005, true);
    CHECK_THAT(priceBermudanPde(itm.problem).value, WithinRel(itm.problem.rights().front().forwardValue(), 2e-3));
    // Out of the money the option is worth less and less: value falls monotonically in the strike.
    // (Each strike is calibrated on its own point of the smile, as a desk would price it.)
    double prev = 1.0;
    for (double strike : {0.02, 0.03, 0.04, 0.05, 0.07, 0.09}) {
        const double v = priceBermudanPde(bermudan10nc1(strike, true).problem).value;
        CHECK(v >= 0.0);
        CHECK(v < prev);
        prev = v;
    }
    CHECK(prev < 0.1 * priceBermudanPde(bermudan10nc1(0.035, true).problem).value);
}

TEST_CASE("Mean reversion: the calibrated Europeans do not move, the Bermudan's switch value rises", "[pricers][pde][kappa]") {
    double prevSwitch = -1.0;
    std::vector<double> reference;
    for (double kappa : {-0.01, 0.0, 0.02, 0.05, 0.10}) {
        const CalibratedBermudan b = bermudan10nc1(0.035, true, kappa);
        const std::vector<double> european = b.problem.europeanValues();
        if (reference.empty()) reference = european;
        for (std::size_t i = 0; i < european.size(); ++i) CHECK_THAT(european[i], WithinRel(reference[i], 1e-8));
        const double sw = priceBermudanPde(b.problem).value - *std::max_element(european.begin(), european.end());
        CHECK(sw > prevSwitch);
        prevSwitch = sw;
    }
}

TEST_CASE("Convergence: second order over four doublings of the grid, with and without Rannacher", "[pricers][pde][convergence]") {
    const CalibratedBermudan b = bermudan10nc1();
    for (bool rannacher : {true, false}) {
        auto price = [&](int n) {
            PdeOptions o;
            o.spaceSteps = n;
            o.timeSteps = n;
            o.rannacher = rannacher;
            return priceBermudanPde(b.problem, o).value;
        };
        const double exact = price(6400);
        const double coarse = std::fabs(price(100) - exact), fine = std::fabs(price(1600) - exact);
        INFO("rannacher " << rannacher << " error at 100: " << coarse << ", at 1600: " << fine);
        // Order p means the error falls by 16^p over a 16-fold refinement. The error is not monotone
        // step to step (it depends on where the exercise boundary sits between nodes), so the order
        // is measured across the whole range: require p > 1.5.
        CHECK(fine < coarse / 64.0);
        CHECK(fine < 1e-6); // under 0.01 bp of notional
    }
}

TEST_CASE("Rannacher start-up: smoother gamma, at some cost in price accuracy", "[pricers][pde][rannacher]") {
    const CalibratedBermudan b = bermudan10nc1();
    PdeOptions ref;
    ref.spaceSteps = 6400;
    ref.timeSteps = 6400;
    const double exact = priceBermudanPde(b.problem, ref).value;
    // A fine space grid with few time steps: the regime where undamped Crank–Nicolson modes survive.
    auto run = [&](bool rannacher) {
        PdeOptions o;
        o.spaceSteps = 3200;
        o.timeSteps = 36;
        o.rannacher = rannacher;
        return priceBermudanPde(b.problem, o);
    };
    // Total variation of the second difference around today's state, relative to its peak: about 2
    // for a smooth single-humped gamma, larger when it oscillates.
    auto roughness = [](const PdeResult& r) {
        const std::size_t c = r.values.size() / 2;
        double tv = 0.0, peak = 0.0;
        for (std::size_t j = c - 400; j < c + 400; ++j) {
            const double g0 = r.values[j - 1] - 2.0 * r.values[j] + r.values[j + 1], g1 = r.values[j] - 2.0 * r.values[j + 1] + r.values[j + 2];
            tv += std::fabs(g1 - g0);
            peak = std::max(peak, std::fabs(g0));
        }
        return tv / peak;
    };
    const PdeResult with = run(true), without = run(false);
    INFO("roughness with " << roughness(with) << " without " << roughness(without));
    CHECK(roughness(with) < 3.0);
    CHECK(roughness(without) > 2.0 * roughness(with));
    CHECK(std::fabs(without.value - exact) < std::fabs(with.value - exact));
}

TEST_CASE("Exercise boundary: payer exercises above a rate that is above the strike", "[pricers][pde][boundary]") {
    const CalibratedBermudan b = bermudan10nc1(0.035, true);
    const PdeResult r = priceBermudanPde(b.problem);
    REQUIRE(r.boundaryRate.size() == 9);
    for (std::size_t i = 0; i + 1 < 9; ++i) {
        INFO("date " << i << " boundary rate " << r.boundaryRate[i]);
        // Exercising gives up the remaining optionality, so the swap must be IN the money by a margin.
        CHECK(r.boundaryRate[i] > 0.035);
        CHECK(r.boundaryRate[i] < 0.07);
    }
    // At the last date there is nothing to wait for: exercise as soon as the swap is worth anything.
    CHECK_THAT(r.boundaryRate[8], WithinAbs(0.035, 2e-5));
    // The receiver mirrors it: exercise below a rate that is below the strike.
    const PdeResult rr = priceBermudanPde(bermudan10nc1(0.035, false).problem);
    for (std::size_t i = 0; i + 1 < 9; ++i) CHECK(rr.boundaryRate[i] < 0.035);
}

TEST_CASE("Grid validation and insensitivity to the grid's width", "[pricers][pde]") {
    const CalibratedBermudan b = bermudan10nc1();
    PdeOptions narrow = fine(), wide = fine();
    narrow.stdDevs = 6.0;
    wide.stdDevs = 9.0;
    wide.spaceSteps = 2400;
    CHECK_THAT(priceBermudanPde(b.problem, narrow).value, WithinAbs(priceBermudanPde(b.problem, wide).value, 1e-8));
    PdeOptions bad;
    bad.spaceSteps = 2;
    CHECK_THROWS(priceBermudanPde(b.problem, bad));
}

TEST_CASE("Calibration strike: the trade's own strike is the default; other strikes reprice their own swaptions", "[pricers][bermudan][smile]") {
    const auto withRate = [](double k) { return spotSwap(10, k); };
    const auto trade = spotSwap(10, 0.035);
    const CalibratedBermudan base = bermudan10nc1();
    const std::vector<double> own(base.problem.size(), 0.035);
    const CalibratedBermudan same = calibrateBermudanAtStrikes(usdMarket(), trade, 1.0, 0, 0.03, withRate, own);
    CHECK_THAT(priceBermudanPde(same.problem).value, WithinAbs(priceBermudanPde(base.problem).value, 1e-14));

    // At the money: each calibration swaption is at its own forward, and is repriced there.
    std::vector<double> atm;
    for (const auto& r : base.problem.rights()) atm.push_back(r.forwardRate());
    const CalibratedBermudan atTheMoney = calibrateBermudanAtStrikes(usdMarket(), trade, 1.0, 0, 0.03, withRate, atm);
    CHECK(atTheMoney.calibration.ok);
    CHECK(atTheMoney.calibration.maxAbsErrorBp < 1e-6);
    for (std::size_t i = 0; i < atm.size(); ++i) CHECK_THAT(atTheMoney.calibration.rows[i].strike, WithinAbs(atm[i], 1e-15));
    // The trade priced is unchanged: same rights, same strike. Only the model's variance differs.
    for (std::size_t i = 0; i < atm.size(); ++i) CHECK_THAT(atTheMoney.problem.rights()[i].strike, WithinAbs(0.035, 1e-15));
    // The market has a smile, so the two calibrations give different models and different prices.
    CHECK(std::fabs(priceBermudanPde(atTheMoney.problem).value - priceBermudanPde(base.problem).value) > 1e-6);
    CHECK_THROWS_AS(calibrateBermudanAtStrikes(usdMarket(), trade, 1.0, 0, 0.03, withRate, {0.035}), rl::InputError);
}
