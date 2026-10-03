#pragma once
// Residuals for the global bootstrap. A residual r_i(z; q_i) is zero exactly when the curves rebuilt
// from z reprice instrument i at its quote q_i.
//
// This library uses RATE FORM: r = model quote − market quote (par rate, par spread, implied rate, or
// (model price − quoted price)/100 for a future). Every residual is then already in rate units, the
// scale to basis points is 1e4 for every instrument, and ∂F/∂q = −I. The alternative, PV form, needs
// each residual divided by its annuity before a 1-week deposit and a 30-year basis swap can be compared.
// The reasoning is at the top of detail/residuals.ipp.

#include "rl/bootstrap/instrument_set.hpp"
#include "rl/market/curve_set.hpp"
#include "rl/market/fixings.hpp"
#include "rl/math/scalar.hpp"

namespace rl::bootstrap {

// --- one residual per instrument type (bodies in detail/residuals.ipp) ------------------------------

template <math::Scalar T>
T residual(const instruments::Deposit& d, double quote, const market::CurveSet<T>& curves, const market::FixingStore& fixings);
template <math::Scalar T>
T residual(const instruments::Fra& f, double quote, const market::CurveSet<T>& curves, const market::FixingStore& fixings);
template <math::Scalar T>
T residual(const instruments::Future& f, double quotedPrice, const market::CurveSet<T>& curves, const market::FixingStore& fixings);
template <math::Scalar T>
T residual(const instruments::Swap& s, double quote, const market::CurveSet<T>& curves, const market::FixingStore& fixings);
template <math::Scalar T>
T residual(const instruments::BasisSwap& b, double quote, const market::CurveSet<T>& curves, const market::FixingStore& fixings);
template <math::Scalar T>
T residual(const TurnOfYear& t, double quote, const market::CurveSet<T>& curves, const market::FixingStore& fixings);

/// Factor that turns residual i into basis points of its quote: 1e4 in rate form.
double residualScale(const BootstrapInstrument& bi, const market::CurveSet<double>& curves, const market::FixingStore& fixings);

/// Model value of a turn in quote units (the jump the curve actually carries over [start, end]).
template <math::Scalar T>
T modelTurn(const TurnOfYear& t, const market::CurveSet<T>& curves);

// --- dispatchers ------------------------------------------------------------------------------------

/// std::visit over the trade variant into the residual() overloads above.
template <math::Scalar T>
T residual(const BootstrapInstrument& bi, const market::CurveSet<T>& curves, const market::FixingStore& fixings);

/// Model value in the *quoted* units for the repricing report: par rate, par spread, implied
/// deposit/FRA rate (decimal), futures price, turn jump. Uses linear:: directly (not the
/// solver's own residual), so the report is an independent check of the solve.
template <math::Scalar T>
T modelQuote(const BootstrapInstrument& bi, const market::CurveSet<T>& curves, const market::FixingStore& fixings);

} // namespace rl::bootstrap
