#pragma once
// Annuity (PV01 of a unit coupon): A = Σ_j τ_j N_j P_disc(t^pay_j) over the coupon cashflows of a leg with
// payment after the reference date. Unsigned (the leg's side is not applied). NotionalExchange
// cashflows do not count. Andersen & Piterbarg vol. 1 §5.3.

#include "rl/instruments/leg.hpp"
#include "rl/market/curve_set.hpp"
#include "rl/math/scalar.hpp"

namespace rl::linear {

template <math::Scalar T>
T annuity(const instruments::Leg& leg, const market::CurveSet<T>& curves);

} // namespace rl::linear
