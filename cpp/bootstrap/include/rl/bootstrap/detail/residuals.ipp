#pragma once
// One residual per instrument type. Instantiated for double in src/residuals.cpp and for the AAD
// scalar in src/instantiate_real.cpp (the Jacobian is taken by an adjoint sweep through these).
#include "rl/bootstrap/residuals.hpp"

#include "rl/core/error.hpp"
#include "rl/linear/annuity.hpp"
#include "rl/linear/forward.hpp"
#include "rl/linear/leg_pv.hpp"
#include "rl/linear/par_rate.hpp"

#include <cmath>
#include <variant>

namespace rl::bootstrap {

// Residual convention: RATE FORM throughout. r_i = (model quote) − (market quote), in rate units:
// implied rate for deposits and FRAs, par rate for swaps, par spread for basis swaps, the jump for a
// turn, and (model price − quoted price)/100 for a future.
//
// Why rate form rather than PV form:
//   * every residual is already "how far off is this instrument, in rate", so a 1-week deposit and a
//     30-year basis swap are commensurable without an annuity scaling (a PV residual for the same
//     1 bp error is ~3e-7 on the overnight deposit and ~1.8e-3 on the 30y basis swap);
//   * ∂F/∂q is −I (−1/100 for futures), which makes the implicit-function step for par risk trivial;
//   * the solver's convergence number is the number in the repricing report.
// The cost is one division by the annuity per swap, which is noise.

/// Deposit: F_impl = (P(t_s)/P(t_e) − 1)/τ on the index's projection curve, against the quoted rate.
template <math::Scalar T>
T residual(const instruments::Deposit& d, double quote, const market::CurveSet<T>& curves, const market::FixingStore&) {
    return linear::impliedDepositRate(d, curves) - quote;
}

/// FRA: the same shape on the term index's projection curve.
template <math::Scalar T>
T residual(const instruments::Fra& f, double quote, const market::CurveSet<T>& curves, const market::FixingStore&) {
    return linear::impliedFraRate(f, curves) - quote;
}

/**
 * SR3 future. The contract settles to 100·(1 − R) with R the compounded overnight rate over the
 * reference quarter. The curve gives the *forward* F over that quarter; the futures rate is F + c,
 * where the convexity adjustment c > 0 exists because the future is margined daily: the long is paid
 * when rates fall (and reinvests low) and pays when rates rise (and funds high), so the futures rate
 * must sit above the forward to compensate.
 * SIMPLIFICATION: c is an input from the quote file. A desk derives it from a short-rate model
 * (under Hull–White it grows roughly like ½σ²·t_s·t_e), which makes it a vol-dependent number.
 * Residual in rate units: (model price − quoted price)/100; a higher curve gives a lower price.
 */
template <math::Scalar T>
T residual(const instruments::Future& f, double quotedPrice, const market::CurveSet<T>& curves, const market::FixingStore&) {
    return (linear::impliedFuturePrice(f, curves) - quotedPrice) / 100.0;
}

/**
 * OIS swap: par rate − quote, with K_par = PV(float leg)/annuity. The float leg is compounded
 * overnight coupons (realised fixings for any elapsed part, then projection-curve DF ratios),
 * discounted on the discount curve. When projection = discount the float leg telescopes to
 * P(t_0) − P(t_n) and this is the textbook formula; for a Fed Funds OIS discounted at SOFR it is not.
 */
template <math::Scalar T>
T residual(const instruments::Swap& s, double quote, const market::CurveSet<T>& curves, const market::FixingStore& fixings) {
    return linear::parRate(s, curves, fixings) - quote;
}

/**
 * Basis swap (Fed Funds + s versus SOFR flat): par spread − quote. This is the residual that couples
 * the curves: it reads Fed Funds forwards, SOFR forwards and SOFR discount factors, so it has
 * non-zero derivatives with respect to the nodes of both curves. The quote file assigns the 25Y and
 * 30Y basis swaps to the SOFR curve (where no SOFR OIS is quoted), so the SOFR long end depends on
 * the Fed Funds curve while the Fed Funds curve depends on SOFR discounting: neither curve can be
 * built first, which is why the solve is global.
 */
template <math::Scalar T>
T residual(const instruments::BasisSwap& b, double quote, const market::CurveSet<T>& curves, const market::FixingStore& fixings) {
    return linear::parSpread(b, curves, fixings, 0) - quote;
}

/// Turn of year: the jump parameter J of the curve's jump layer against the quoted jump. Linear in J,
/// but every instrument whose accrual spans the turn window sees J through its forwards.
template <math::Scalar T>
T residual(const TurnOfYear& t, double quote, const market::CurveSet<T>& curves, const market::FixingStore&) {
    return modelTurn(t, curves) - quote;
}

/// The jump the curve carries over the turn window, matched by the window's start and end times.
template <math::Scalar T>
T modelTurn(const TurnOfYear& t, const market::CurveSet<T>& curves) {
    const market::YieldCurve<T>& curve = curves.projection(t.index);
    const double a = curve.time(t.start), b = curve.time(t.end);
    for (const auto& jump : curve.turns())
        if (std::abs(jump.start - a) < 1e-9 && std::abs(jump.end - b) < 1e-9) return jump.jump;
    throw InputError("turn " + t.id + " has no matching jump on curve " + curve.name());
}

// --- dispatchers -------------------------------------------------------------------------------

template <math::Scalar T>
T residual(const BootstrapInstrument& bi, const market::CurveSet<T>& curves, const market::FixingStore& fixings) {
    return withContext("residual " + bi.id, [&] {
        return std::visit([&](const auto& trade) -> T { return residual(trade, bi.quote, curves, fixings); }, bi.trade);
    });
}

template <math::Scalar T>
T modelQuote(const BootstrapInstrument& bi, const market::CurveSet<T>& curves, const market::FixingStore& fixings) {
    struct Visitor {
        const market::CurveSet<T>& curves;
        const market::FixingStore& fixings;
        T operator()(const instruments::Deposit& d) const { return linear::impliedDepositRate(d, curves); }
        T operator()(const instruments::Fra& f) const { return linear::impliedFraRate(f, curves); }
        T operator()(const instruments::Future& f) const { return linear::impliedFuturePrice(f, curves); }
        T operator()(const instruments::Swap& s) const { return linear::parRate(s, curves, fixings); }
        T operator()(const instruments::BasisSwap& b) const { return linear::parSpread(b, curves, fixings, 0); }
        T operator()(const TurnOfYear& t) const { return modelTurn(t, curves); }
    };
    return withContext("model quote " + bi.id, [&] { return std::visit(Visitor{curves, fixings}, bi.trade); });
}

} // namespace rl::bootstrap

#define RL_BOOTSTRAP_INSTANTIATE_RESIDUALS(T)                                                                                   \
    namespace rl::bootstrap {                                                                                                   \
    template T residual<T>(const instruments::Deposit&, double, const market::CurveSet<T>&, const market::FixingStore&);        \
    template T residual<T>(const instruments::Fra&, double, const market::CurveSet<T>&, const market::FixingStore&);            \
    template T residual<T>(const instruments::Future&, double, const market::CurveSet<T>&, const market::FixingStore&);         \
    template T residual<T>(const instruments::Swap&, double, const market::CurveSet<T>&, const market::FixingStore&);           \
    template T residual<T>(const instruments::BasisSwap&, double, const market::CurveSet<T>&, const market::FixingStore&);      \
    template T residual<T>(const TurnOfYear&, double, const market::CurveSet<T>&, const market::FixingStore&);                  \
    template T modelTurn<T>(const TurnOfYear&, const market::CurveSet<T>&);                                                     \
    template T residual<T>(const BootstrapInstrument&, const market::CurveSet<T>&, const market::FixingStore&);                 \
    template T modelQuote<T>(const BootstrapInstrument&, const market::CurveSet<T>&, const market::FixingStore&);               \
    }
