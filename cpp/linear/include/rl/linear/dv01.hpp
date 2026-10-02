#pragma once
// DV01s.
//   fixedRateDv01(swap)  = ∂PV/∂K × 1e-4 = sign_fix · A_fix · 1e-4          (analytic; the number a desk calls "the DV01")
//   parallelDv01(trade)  = PV(curves bumped +1bp in every zero rate) − PV(curves)   (one-sided bump-and-revalue)
//   parallelBump(curves, dr) returns a CurveSet with every curve's ordinates y_i = r_i t_i shifted by dr·t_i,
//   i.e. a parallel shift of dr in continuously-compounded zero rates. Used by tests and the demo; the risk
//   module has the real ladders.

#include "rl/instruments/basis_swap.hpp"
#include "rl/instruments/swap.hpp"
#include "rl/market/curve_set.hpp"
#include "rl/market/fixings.hpp"
#include "rl/math/scalar.hpp"

namespace rl::linear {

template <math::Scalar T>
T fixedRateDv01(const instruments::Swap& s, const market::CurveSet<T>& curves);

template <math::Scalar T>
market::CurveSet<T> parallelBump(const market::CurveSet<T>& curves, double rateBump);

template <math::Scalar T>
T parallelDv01(const instruments::Swap& s, const market::CurveSet<T>& curves, const market::FixingStore& fixings, double bumpBp = 1.0);
template <math::Scalar T>
T parallelDv01(const instruments::BasisSwap& b, const market::CurveSet<T>& curves, const market::FixingStore& fixings, double bumpBp = 1.0);

} // namespace rl::linear
