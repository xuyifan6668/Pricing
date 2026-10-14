#include "rl/vol/conversion.hpp"

#include "rl/core/error.hpp"
#include "rl/vol/bachelier.hpp"
#include "rl/vol/black.hpp"

#include <cmath>
#include <string>

namespace rl::vol {

// Exact conversion between the two quoting conventions: the two vols are defined to give the SAME
// undiscounted option price, so each direction is "price under one model, invert under the other".
//
//   normal (Bachelier):            dF = σ_N dW           V = ω(F−K)Φ(ωd) + σ_N√T φ(d),  d = (F−K)/(σ_N√T)
//   shifted lognormal (Black):     dF = σ_LN (F+s) dW    V = ω[(F+s)Φ(ωd₁) − (K+s)Φ(ωd₂)]
//
// Both directions work on the OUT-OF-THE-MONEY option. By put–call parity the answer is the same for
// the call and the put, but the in-the-money price is intrinsic plus time value, and inverting it
// loses the digits the intrinsic occupies.
//
// The shift s is a quoting convention published with the screen (it exists so that lognormal vols can
// be quoted when rates are near or below zero); it is never a free parameter chosen to make a
// conversion succeed.

namespace {
OptionType otm(double forward, double strike) { return strike >= forward ? OptionType::Call : OptionType::Put; }

void requireDomain(double forward, double strike, double expiry, double vol, double shift, const char* fn) {
    RL_REQUIRE(std::isfinite(forward) && std::isfinite(strike) && std::isfinite(vol) && std::isfinite(shift), std::string(fn) + ": inputs must be finite");
    RL_REQUIRE(expiry > 0.0, std::string(fn) + ": expiry must be positive (no time value to match at T = 0)");
    RL_REQUIRE(vol >= 0.0, std::string(fn) + ": vol must be non-negative");
    RL_REQUIRE(forward + shift > 0.0, std::string(fn) + ": forward + shift must be positive, or the lognormal model does not exist");
    RL_REQUIRE(strike + shift > 0.0, std::string(fn) + ": strike + shift must be positive, or the lognormal price is intrinsic for every vol");
}
} // namespace

// Normal → shifted lognormal. This direction can FAIL: an out-of-the-money lognormal option is worth
// less than min(F+s, K+s) however large the vol (the shifted forward cannot go below zero), while the
// normal price is unbounded. A low-strike put with a high normal vol can therefore have a price no
// lognormal vol reproduces; that is reported, not clamped.
double normalToShiftedLognormalVol(double forward, double strike, double expiry, double normalVol, double shift) {
    requireDomain(forward, strike, expiry, normalVol, shift, "normalToShiftedLognormalVol");
    if (normalVol == 0.0) return 0.0;
    const OptionType type = otm(forward, strike);
    const double price = bachelierPrice(forward, strike, expiry, normalVol, type);
    const double bound = type == OptionType::Call ? forward + shift : strike + shift;
    if (!(price < bound))
        throw InputError("normalToShiftedLognormalVol: the normal price " + std::to_string(price) + " is not attainable by a lognormal model (bound " +
                         std::to_string(bound) + " at forward " + std::to_string(forward) + ", strike " + std::to_string(strike) + ", shift " +
                         std::to_string(shift) + ")");
    return blackImpliedVol(price, forward, strike, expiry, shift, type);
}

// Shifted lognormal → normal. Always attainable: the Bachelier price is unbounded in σ_N, so every
// lognormal price has a normal vol. That asymmetry is why normal vol is the robust quoting convention
// for rates, and why the cube and the LGM calibration in this library work in normal vol.
double shiftedLognormalToNormalVol(double forward, double strike, double expiry, double lognormalVol, double shift) {
    requireDomain(forward, strike, expiry, lognormalVol, shift, "shiftedLognormalToNormalVol");
    if (lognormalVol == 0.0) return 0.0;
    const OptionType type = otm(forward, strike);
    const double price = blackPrice(forward, strike, expiry, lognormalVol, shift, type);
    return bachelierImpliedVol(price, forward, strike, expiry, type);
}

// --- ATM series ------------------------------------------------------------------------------------------------

double atmNormalFromLognormalApprox(double forward, double expiry, double lognormalVol, double shift) {
    const double f = forward + shift;
    RL_REQUIRE(f > 0.0, "atmNormalFromLognormalApprox: forward + shift must be positive");
    RL_REQUIRE(expiry >= 0.0 && lognormalVol >= 0.0, "atmNormalFromLognormalApprox: expiry and vol must be non-negative");
    return lognormalVol * f * (1.0 - lognormalVol * lognormalVol * expiry / 24.0);
}

double atmLognormalFromNormalApprox(double forward, double expiry, double normalVol, double shift) {
    const double f = forward + shift;
    RL_REQUIRE(f > 0.0, "atmLognormalFromNormalApprox: forward + shift must be positive");
    RL_REQUIRE(expiry >= 0.0 && normalVol >= 0.0, "atmLognormalFromNormalApprox: expiry and vol must be non-negative");
    const double s0 = normalVol / f;
    return s0 * (1.0 + s0 * s0 * expiry / 24.0);
}

} // namespace rl::vol
