#pragma once
// Par rate and par spread.
//   parRate(swap):   the single fixed rate K such that PV(swap) = 0, i.e. K = −sign_fix · PV(float leg) / A_fix
//                    (Andersen & Piterbarg vol. 1 eq. 5.10 in multi-curve form: the float leg is projected
//                    off its index curve, the annuity discounted off the discount curve).
//   parSpread(basis, i): the spread s on leg i such that PV = 0, assuming all coupons on that leg carry
//                    the same spread: s = s_current − PV(basis) / (sign_i · A_i).

#include "rl/instruments/basis_swap.hpp"
#include "rl/instruments/swap.hpp"
#include "rl/market/curve_set.hpp"
#include "rl/market/fixings.hpp"
#include "rl/math/scalar.hpp"

#include <cstddef>

namespace rl::linear {

template <math::Scalar T>
T parRate(const instruments::Swap& s, const market::CurveSet<T>& curves, const market::FixingStore& fixings);

/// `spreadLeg` is 0 or 1 (the leg whose spread is solved for; the builders put the spread leg first).
/// Throws InputError unless every overnight/term coupon on that leg carries the same spread.
template <math::Scalar T>
T parSpread(const instruments::BasisSwap& b, const market::CurveSet<T>& curves, const market::FixingStore& fixings,
            std::size_t spreadLeg = 0);

} // namespace rl::linear
