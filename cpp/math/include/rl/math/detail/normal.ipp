#pragma once
// standard normal pdf, cdf and the AS241 inverse. Coefficients from Wichura
// (1988), checked against the R sources (qnorm.c).
#include "rl/core/error.hpp"
#include "rl/math/normal.hpp"

#include <cmath>

namespace rl::math {

/**
 * φ(x) = (2π)^{−1/2} exp(−x²/2). Underflows to 0 for |x| > ~38.6; that is fine. Use a constexpr
 * 1/√(2π) = 0.3989422804014327; do not call sqrt at runtime. Test: φ(0) = 0.3989422804014327,
 * symmetry, and ∫φ = 1 by Gauss–Hermite.
 */
template <Scalar T>
T normalPdf(const T& x) {
    constexpr double invSqrt2Pi = 0.3989422804014327;
    return invSqrt2Pi * exp(-0.5 * x * x); // exp underflows to 0 for |x| > ~38.6, which is the right answer
}

/**
 * Φ(x) = ½ erfc(−x/√2).
 *
 * Why erfc and not ½(1 + erf(x/√2)): for x = −8, erf ≈ −1 + 1.2e-15 and the addition loses every
 * significant digit; erfc(8/√2 = 5.657) = 1.2e-15 keeps them all. Bachelier vega at deep OTM strikes
 * and the Andersen–Broadie exercise probabilities live in that tail. Both tails: Φ(−40) ≈ 3.7e-350
 * underflows to 0, Φ(40) = 1 exactly; neither is an error.
 * Test: Φ(0) = 0.5; Φ(x) + Φ(−x) = 1 to 1e-16; Φ(−1.959963984540054) = 0.025; AS241 round trip.
 */
template <Scalar T>
T normalCdf(const T& x) {
    // Branches on value(x): the tails are constants, so no derivative is lost there.
    if (value(x) >= 40.0) return T(1.0);
    if (value(x) <= -40.0) return T(0.0);
    constexpr double invSqrt2 = 0.70710678118654752440;
    return 0.5 * erfc(-x * invSqrt2);
}

/**
 * Φ⁻¹(p), AS241 PPND16 (Wichura 1988, Applied Statistics 37(3), 477–484). Three rational
 * approximations in q = p − ½ and r:
 *
 *   |q| ≤ 0.425:  r = 0.180625 − q²;               x = q · A(r)/B(r)            (central region)
 *   else:         r = min(p, 1 − p);  r = √(−log r)
 *       r ≤ 5:    r ← r − 1.6;                     x = C(r)/D(r)                (tail)
 *       r > 5:    r ← r − 5;                       x = E(r)/F(r)                (far tail)
 *       x ← −x if q < 0
 *
 * A..F are degree-7 polynomials; the 48 coefficients are in the paper and in every reference
 * implementation (e.g. the R source `qnorm.c`, or Acklam's page). Copy them exactly; a single typo
 * shows up only at the 1e-9 level in one region, which is why the test sweeps p on a log grid.
 * Domain: p ∈ (0, 1). Throw InputError for p ≤ 0 or p ≥ 1 (a Sobol point at exactly 0 is the usual
 * culprit — see sobol.hpp on skipping the origin). For p = 1e-300, x ≈ −37.05; the far-tail branch
 * handles it.
 * Hot loop: one call per gaussian draw in MC. Branch on value(p); polynomial evaluation in T.
 * Test: cdf(inv(p)) = p to 1e-15 relative for p ∈ [1e-300, 1 − 1e-16]; inv(0.975) = 1.959963984540054;
 * antisymmetry inv(1 − p) = −inv(p).
 */
template <Scalar T>
T normalInverseCdf(const T& p) {
    RL_REQUIRE(value(p) > 0.0 && value(p) < 1.0, "normalInverseCdf: p must be in (0, 1)");

    T q = p - T(0.5);
    T val;

    if (std::fabs(value(q)) <= 0.425) {
        // Central region: |q| <= 0.425
        T r = T(0.180625) - q * q;
        T num = ((((((r * T(2509.0809287301226727) +
                     T(33430.575583588128105))  * r + T(67265.770927008700853)) * r +
                     T(45921.953931549871457))  * r + T(13731.693765509461125)) * r +
                     T(1971.5909503065514427))  * r + T(133.14166789178437745)) * r +
                     T(3.387132872796366608);
        T den = ((((((r * T(5226.495278852854561) +
                     T(28729.085735721942674))  * r + T(39307.89580009271061))  * r +
                     T(21213.794301586595867))  * r + T(5394.1960214247511077)) * r +
                     T(687.1870074920579083))   * r + T(42.313330701600911252)) * r +
                     T(1.0);
        val = q * num / den;
    } else {
        // Tail regions
        T r = (value(q) < 0.0) ? p : (T(1.0) - p);
        r = sqrt(-log(r));

        if (value(r) <= 5.0) {
            // Tail: r <= 5
            r = r - T(1.6);
            T num = ((((((r * T(7.7454501427834140764e-4) +
                         T(0.0227238449892691845833)) * r + T(0.24178072517745061177)) * r +
                         T(1.27045825245236838258))    * r + T(3.64784832476320460504))  * r +
                         T(5.7694972214606914055))     * r + T(4.6303378461565452959))   * r +
                         T(1.42343711074968357734);
            T den = ((((((r * T(1.05075007164441684324e-9) +
                         T(5.475938084995344946e-4))   * r + T(0.0151986665636164571966)) * r +
                         T(0.14810397642748007459))    * r + T(0.68976733498510000455))   * r +
                         T(1.6763848301838038494))     * r + T(2.05319162663775882187))   * r +
                         T(1.0);
            val = num / den;
        } else {
            // Far tail: r > 5
            r = r - T(5.0);
            T num = ((((((r * T(2.01033439929228813265e-7) +
                         T(2.71155556874348757815e-5)) * r + T(0.0012426609473880784386)) * r +
                         T(0.026532189526576123093))   * r + T(0.29656057182850489123))    * r +
                         T(1.7848265399172913358))     * r + T(5.4637849111641143699))     * r +
                         T(6.6579046435011037772);
            T den = ((((((r * T(2.04426310338993978564e-15) +
                         T(1.4215117583164458887e-7))  * r + T(1.8463183175100546818e-5))  * r +
                         T(7.868691311456132591e-4))   * r + T(0.0148753612908506148525))  * r +
                         T(0.13692988092273580531))    * r + T(0.59983220655588793769))    * r +
                         T(1.0);
            val = num / den;
        }

        if (value(q) < 0.0) val = -val;
    }

    return val;
}

} // namespace rl::math

#define RL_MATH_INSTANTIATE_NORMAL(T)                                                                  \
    template T rl::math::normalPdf<T>(const T&);                                                       \
    template T rl::math::normalCdf<T>(const T&);                                                       \
    template T rl::math::normalInverseCdf<T>(const T&);
