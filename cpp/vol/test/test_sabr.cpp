// Hagan SABR: analytic limits, ATM continuity, symmetry, shift identity, cross-checks against exact conversions,
// the Dual gradients against finite differences, the ATM constraint, and Real.
#include "rl/vol/sabr.hpp"

#include "rl/aad/real.hpp"
#include "rl/core/error.hpp"
#include "rl/vol/bachelier.hpp"
#include "rl/vol/black.hpp"
#include "test_support.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <algorithm>
#include <cmath>

using namespace rl::vol;
using rl::aad::Real;
using rl::market::SabrParams;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {
SabrParams<double> sp(double alpha, double beta, double rho, double nu, double shift = 0.0) {
    SabrParams<double> p;
    p.alpha = alpha;
    p.rho = rho;
    p.nu = nu;
    p.beta = beta;
    p.shift = shift;
    return p;
}
/// Exact normal vol of a shifted-Black price: pieces only.
double exactNormalOfBlack(double f, double k, double tt, double sl, double shift) {
    const OptionType t = k >= f ? OptionType::Call : OptionType::Put;
    return bachelierImpliedVol(blackPrice(f, k, tt, sl, shift, t), f, k, tt, t);
}
} // namespace

TEST_CASE("SABR beta = 0 at the money: alpha (1 + (2 - 3 rho^2) nu^2 T / 24)", "[vol][sabr][limits]") {
    for (double rho : {-0.6, 0.0, 0.3}) {
        const auto p = sp(0.0085, 0.0, rho, 0.45);
        const double expected = 0.0085 * (1.0 + (2.0 - 3.0 * rho * rho) * 0.45 * 0.45 * 3.0 / 24.0);
        CHECK_THAT(sabrNormalVol(p, 0.03, 0.03, 3.0), WithinRel(expected, 1e-15));
        CHECK_THAT(sabrAtmNormalVol(p, 0.03, 3.0), WithinRel(expected, 1e-15));
    }
}

TEST_CASE("SABR nu = 0: the local-vol leading term with its alpha^2 correction", "[vol][sabr][limits]") {
    const double f = 0.035, tt = 2.0, a = 0.05, b = 0.5, c = 1.0 - b;
    for (double off : {-0.01, -0.005, 0.005, 0.01}) {
        const double k = f + off;
        const double lead = a * c * (f - k) / (std::pow(f, c) - std::pow(k, c));
        const double corr = -b * (2.0 - b) * a * a / (24.0 * std::pow(f * k, c));
        CHECK_THAT(sabrNormalVol(sp(a, b, -0.3, 0.0), f, k, tt), WithinRel(lead * (1.0 + corr * tt), 1e-13));
    }
}

TEST_CASE("SABR is continuous and differentiable through the money for every beta", "[vol][sabr][limits]") {
    for (double beta : {0.0, 0.5, 1.0}) {
        const auto p = sp(beta == 0.0 ? 0.009 : (beta == 1.0 ? 0.25 : 0.047), beta, -0.25, 0.5);
        const double f = 0.036, tt = 1.0;
        const double v0 = sabrNormalVol(p, f, f, tt);
        const double up = sabrNormalVol(p, f, f + 1e-9, tt), dn = sabrNormalVol(p, f, f - 1e-9, tt);
        CHECK(std::fabs(up - v0) < 1e-9);
        CHECK(std::fabs(dn - v0) < 1e-9);
        const double skewTiny = (up - dn) / 2e-9;
        const double skewFd = (sabrNormalVol(p, f, f + 1e-5, tt) - sabrNormalVol(p, f, f - 1e-5, tt)) / 2e-5;
        INFO("beta " << beta);
        CHECK_THAT(skewTiny, WithinAbs(skewFd, 1e-6));
        CHECK_THAT(v0, WithinRel(sabrAtmNormalVol(p, f, tt), 1e-14));
    }
}

TEST_CASE("SABR zeta/x(zeta) is continuous across its series threshold on both sides", "[vol][sabr][limits]") {
    // β = 0 and α = ν make ζ = F − K, so the threshold |ζ| = 1e-6 sits at |F − K| = 1 bp·0.01.
    const double a = 0.01, rho = -0.3, tt = 1.0, f = 0.03;
    const auto p = sp(a, 0.0, rho, a);
    const double scale = a * (1.0 + (2.0 - 3.0 * rho * rho) * a * a * tt / 24.0);
    for (double sign : {1.0, -1.0}) {
        const double z1 = sign * 0.999e-6, z2 = sign * 1.001e-6;
        const double v1 = sabrNormalVol(p, f, f - z1, tt), v2 = sabrNormalVol(p, f, f - z2, tt);
        const double zm = 0.5 * (z1 + z2);
        const double expected = scale * (-0.5 * rho + (2.0 - 3.0 * rho * rho) * zm / 6.0) * (z2 - z1);
        CHECK_THAT(v2 - v1, WithinAbs(expected, 1e-16));
    }
}

TEST_CASE("SABR with rho = 0 is symmetric under F <-> K", "[vol][sabr]") {
    for (double beta : {0.0, 0.5, 1.0}) {
        const auto p = sp(beta == 0.0 ? 0.009 : (beta == 1.0 ? 0.25 : 0.047), beta, 0.0, 0.6);
        CHECK_THAT(sabrNormalVol(p, 0.030, 0.042, 2.0), WithinRel(sabrNormalVol(p, 0.042, 0.030, 2.0), 1e-14));
    }
    // And negative rho skews the smile down to the right.
    const auto q = sp(0.009, 0.0, -0.4, 0.5);
    CHECK(sabrNormalVol(q, 0.03, 0.025, 1.0) > sabrNormalVol(q, 0.03, 0.035, 1.0));
}

TEST_CASE("SABR shift is a displaced forward and strike; negative rates work with a shift or beta = 0", "[vol][sabr]") {
    const auto shifted = sp(0.08, 0.5, -0.2, 0.4, 0.02);
    const auto plain = sp(0.08, 0.5, -0.2, 0.4, 0.0);
    CHECK_THAT(sabrNormalVol(shifted, -0.002, 0.004, 3.0), WithinRel(sabrNormalVol(plain, 0.018, 0.024, 3.0), 1e-14));
    CHECK(sabrNormalVol(shifted, -0.002, -0.005, 3.0) > 0.0);
    CHECK(sabrNormalVol(sp(0.006, 0.0, -0.2, 0.4), -0.004, -0.009, 3.0) > 0.0);
    CHECK_THROWS_AS(sabrNormalVol(plain, -0.002, 0.004, 3.0), rl::InputError);
    CHECK_THROWS_AS(sabrNormalVol(sp(0.08, 0.5, 1.0, 0.4), 0.03, 0.03, 1.0), rl::InputError);
}

TEST_CASE("SABR beta = 1, nu = 0 is Black: the normal form matches the exact conversion to the expansion order", "[vol][sabr][crosscheck]") {
    const double f = 0.03, tt = 1.0, a = 0.2;
    for (double m : {0.8, 0.9, 1.0, 1.1, 1.25}) {
        const double k = f * m;
        INFO("K/F " << m);
        CHECK_THAT(sabrNormalVol(sp(a, 1.0, 0.0, 0.0), f, k, tt), WithinRel(exactNormalOfBlack(f, k, tt, a, 0.0), 1e-5));
        CHECK_THAT(sabrLognormalVol(sp(a, 1.0, 0.0, 0.0), f, k, tt), WithinRel(a, 1e-15));
    }
}

TEST_CASE("SABR lognormal (2.17a) and normal (B.69) forms agree to a small fraction of a bp at 1Y", "[vol][sabr][crosscheck]") {
    const auto p = sp(0.045, 0.5, -0.2, 0.4);
    const double f = 0.036, tt = 1.0;
    for (double off : {-0.01, -0.005, 0.0, 0.005, 0.01}) {
        const double k = f + off;
        const double viaLognormal = exactNormalOfBlack(f, k, tt, sabrLognormalVol(p, f, k, tt), 0.0);
        CHECK_THAT(sabrNormalVol(p, f, k, tt) * 1e4, WithinAbs(viaLognormal * 1e4, 0.1));
    }
}

TEST_CASE("SABR gradient in (alpha, rho, nu) matches central finite differences", "[vol][sabr][gradient]") {
    for (const auto& node : {test::kBenign, test::kStressed}) {
        const auto p = node.params();
        for (double off : test::kOffsets) {
            const double k = node.forward + off;
            const auto g = sabrNormalVolGradient(p, node.forward, k, node.expiry);
            CHECK(g.vol == sabrNormalVol(p, node.forward, k, node.expiry));
            auto fd = [&](int which) {
                const double h = 1e-6;
                auto up = p, dn = p;
                if (which == 0) { up.alpha += h; dn.alpha -= h; }
                if (which == 1) { up.rho += h; dn.rho -= h; }
                if (which == 2) { up.nu += h; dn.nu -= h; }
                return (sabrNormalVol(up, node.forward, k, node.expiry) - sabrNormalVol(dn, node.forward, k, node.expiry)) / (2.0 * h);
            };
            INFO("offset " << off);
            CHECK_THAT(g.dAlpha, WithinRel(fd(0), 1e-6) || WithinAbs(fd(0), 1e-10));
            CHECK_THAT(g.dRho, WithinRel(fd(1), 1e-6) || WithinAbs(fd(1), 1e-10));
            CHECK_THAT(g.dNu, WithinRel(fd(2), 1e-6) || WithinAbs(fd(2), 1e-10));
        }
    }
}

TEST_CASE("alphaFromAtmNormalVol inverts the ATM vol to 1e-14 and fails loudly when the ATM vol is unattainable", "[vol][sabr][atm]") {
    int solved = 0, unattainable = 0;
    for (double beta : {0.0, 0.3, 0.5, 1.0}) {
        for (double rho : {-0.7, 0.0, 0.5}) {
            for (double nu : {0.0, 0.3, 1.0}) {
                for (double tt : {1.0 / 12.0, 1.0, 10.0}) {
                    const double f = 0.036, atm = 0.0085;
                    INFO("beta " << beta << " rho " << rho << " nu " << nu << " T " << tt);
                    try {
                        const double a = alphaFromAtmNormalVol(atm, f, tt, beta, rho, nu, 0.0);
                        CHECK(a > 0.0);
                        CHECK_THAT(sabrAtmNormalVol(sp(a, beta, rho, nu), f, tt), WithinRel(atm, 1e-14));
                        ++solved;
                    } catch (const rl::NumericalError&) {
                        // Legitimate only if no α reaches the ATM vol: scan the cubic by brute force.
                        const double scale = atm / std::pow(f, beta);
                        double best = 0.0;
                        for (int i = 1; i <= 20000; ++i) best = std::max(best, sabrAtmNormalVol(sp(1e-3 * i * scale, beta, rho, nu), f, tt));
                        CHECK(best < atm * (1.0 + 1e-9));
                        ++unattainable;
                    }
                }
            }
        }
    }
    CHECK(solved > 3 * unattainable); // long expiry, β near 1 and large ν cap the ATM vol; most of the grid is reachable
    // β = 1, ν = 0, F = 3%, T = 30y: the ATM cubic peaks at ≈ 103 bp; 150 bp has no α.
    CHECK_THROWS_AS(alphaFromAtmNormalVol(0.0150, 0.03, 30.0, 1.0, 0.0, 0.0, 0.0), rl::NumericalError);
    CHECK_THROWS_AS(alphaFromAtmNormalVol(-0.01, 0.03, 1.0, 0.5, 0.0, 0.3, 0.0), rl::InputError);
}

TEST_CASE("ATM-constrained gradient: total derivatives through alpha(rho, nu) match re-solve-and-bump", "[vol][sabr][gradient][atm]") {
    const auto& node = test::kBenign;
    for (double off : test::kOffsets) {
        const double k = node.forward + off;
        const auto g = sabrNormalVolGradientAtmConstrained(node.atm, node.forward, k, node.expiry, node.beta, node.rho, node.nu, node.shift);
        auto vol = [&](double rho, double nu) {
            const double a = alphaFromAtmNormalVol(node.atm, node.forward, node.expiry, node.beta, rho, nu, node.shift);
            return sabrNormalVol(sp(a, node.beta, rho, nu), node.forward, k, node.expiry);
        };
        const double h = 1e-6;
        const double dRho = (vol(node.rho + h, node.nu) - vol(node.rho - h, node.nu)) / (2.0 * h);
        const double dNu = (vol(node.rho, node.nu + h) - vol(node.rho, node.nu - h)) / (2.0 * h);
        INFO("offset " << off);
        CHECK_THAT(g.vol, WithinRel(vol(node.rho, node.nu), 1e-15));
        CHECK_THAT(g.dRho, WithinRel(dRho, 1e-6) || WithinAbs(dRho, 1e-10));
        CHECK_THAT(g.dNu, WithinRel(dNu, 1e-6) || WithinAbs(dNu, 1e-10));
        if (off == 0.0) {
            CHECK_THAT(g.dRho, WithinAbs(0.0, 1e-12)); // the ATM vol is pinned
            CHECK_THAT(g.dNu, WithinAbs(0.0, 1e-12));
        }
    }
}

TEST_CASE("SABR under Real: bit-identical to double, active through the ATM constraint", "[vol][sabr][aad]") {
    const auto& node = test::kStressed;
    const Real atm = Real::input(node.atm), fwd = Real::input(node.forward), rho = Real::input(node.rho), nu = Real::input(node.nu);
    const Real a = alphaFromAtmNormalVol(atm, fwd, node.expiry, node.beta, rho, nu, node.shift);
    CHECK(a.isActive());
    CHECK(a.value() == alphaFromAtmNormalVol(node.atm, node.forward, node.expiry, node.beta, node.rho, node.nu, node.shift));
    SabrParams<Real> pr;
    pr.alpha = a;
    pr.rho = rho;
    pr.nu = nu;
    pr.beta = node.beta;
    pr.shift = node.shift;
    for (double off : test::kOffsets) {
        const Real v = sabrNormalVol(pr, fwd, node.forward + off, node.expiry);
        CHECK(v.isActive());
        CHECK(v.value() == sabrNormalVol(node.params(), node.forward, node.forward + off, node.expiry));
    }
    rl::aad::Tape::current().reset();
}
