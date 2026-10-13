#pragma once
// Shared fixtures for the vol tests: the shipped USD vol screen, and two reference nodes taken from the generating
// surface in tools/gen_vol_quotes.py (numbers rounded; the tests only need them to be representative).
#include "rl/core/calendar.hpp" // core::dataDir
#include "rl/market/sabr_params.hpp"
#include "rl/market/vol_quotes.hpp"
#include "rl/vol/sabr.hpp"

#include <vector>

namespace rl::vol::test {

inline market::SwaptionVolQuotes usdVolQuotes() {
    return market::SwaptionVolQuotes::fromFile(core::dataDir() / "vol" / "usd_swaption_cube_2026-09-04.json");
}

struct NodeSpec {
    double forward;
    double expiry;
    double atm;
    double beta;
    double rho;
    double nu;
    double shift;

    market::SabrParams<double> params() const {
        market::SabrParams<double> p;
        p.rho = rho;
        p.nu = nu;
        p.beta = beta;
        p.shift = shift;
        p.alpha = alphaFromAtmNormalVol<double>(atm, forward, expiry, beta, rho, nu, shift);
        return p;
    }
};

/// A benign 1Yx10Y USD node: 88 bp ATM, mild skew.
inline const NodeSpec kBenign{0.0360, 1.0, 0.0088, 0.5, -0.13, 0.41, 0.0};
/// The stressed 10Yx30Y node: nu = 0.70 on a 10-year expiry. Fits perfectly; density negative ~170 bp below F.
inline const NodeSpec kStressed{0.0383, 10.0, 0.00703, 0.5, -0.45, 0.70, 0.0};

/// The vol screen's strike offsets, ATM included.
inline const std::vector<double> kOffsets{-0.01, -0.005, -0.0025, 0.0, 0.0025, 0.005, 0.01};

} // namespace rl::vol::test
