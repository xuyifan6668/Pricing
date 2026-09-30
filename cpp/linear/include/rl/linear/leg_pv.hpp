#pragma once
// Present value of cashflows, legs and linear trades off a CurveSet<T> (closed forms, implemented;
// sources in detail/linear.ipp). Called ~1e6 times per bootstrap solve, so: no allocation, no virtual
// calls, one std::visit per cashflow, two df() calls per projected coupon.
//
// Conventions (Andersen & Piterbarg, *Interest Rate Modeling* vol. 1, §5.1–5.3 and §6.2; ISDA 2021
// Definitions §7 for overnight compounding):
//   - Discounting is always off curves.discount(); projection off curves.projection(index.name).
//   - A cashflow whose payment date is on or before the curve's reference date is worth zero (settled).
//   - An overnight coupon's rate is compounded (or averaged, per Index::compounding) over its observation
//     period: realised fixings for observation days before the reference date, the projection curve's
//     discount-factor ratio for the rest. Lookback is applied as an *observation shift* (weights taken
//     from the observation dates); lockout freezes the last k fixings at the fixing of the k-th last day.
//   - A term coupon fixed before the reference date takes its fixing from the store; otherwise the simple
//     forward P(start)/P(end) − 1 over the coupon's tau.
//   - Legs carry the pay/receive sign; a swap PV is the signed sum of its legs.
//
// Historical fixings are passive doubles; everything else downstream of a curve node is T.

#include "rl/instruments/basis_swap.hpp"
#include "rl/instruments/cashflow.hpp"
#include "rl/instruments/deposit.hpp"
#include "rl/instruments/fra.hpp"
#include "rl/instruments/future.hpp"
#include "rl/instruments/leg.hpp"
#include "rl/instruments/swap.hpp"
#include "rl/market/curve_set.hpp"
#include "rl/market/fixings.hpp"
#include "rl/math/scalar.hpp"

namespace rl::linear {

using instruments::BasisSwap;
using instruments::Cashflow;
using instruments::Deposit;
using instruments::FixedCoupon;
using instruments::Fra;
using instruments::Future;
using instruments::Leg;
using instruments::NotionalExchange;
using instruments::OvernightCoupon;
using instruments::Swap;
using instruments::TermCoupon;
using market::CurveSet;
using market::FixingStore;

// --- coupon rates -------------------------------------------------------------------------------------

/// The compounded (or averaged) overnight rate R over the coupon's observation period, before spread:
///   compounded:  R = [ Π_{i∈past}(1 + r_i τ_i) · P(t_a)/P(t_e) − 1 ] / τ_obs
///   averaged:    R = [ Σ_{i∈past} r_i τ_i + log(P(t_a)/P(t_e)) ] / τ_obs     (continuous-average approximation)
/// where r_i are realised fixings, τ_i the ACT/360 (index day count) weights from observation day i to
/// the next business day, t_a the first observation day on or after the reference date and t_e the
/// observation end. Throws InputError if a required fixing is missing (context: the index and date).
template <math::Scalar T>
T overnightRate(const OvernightCoupon& c, const CurveSet<T>& curves, const FixingStore& fixings);

/// Realised fixing if the fixing date is before the reference date (or a fixing exists for it),
/// otherwise the simple forward (P(start)/P(end) − 1)/tau off the projection curve.
template <math::Scalar T>
T termRate(const TermCoupon& c, const CurveSet<T>& curves, const FixingStore& fixings);

// --- cashflow and leg PV (unsigned by leg side; the Leg overload applies sign()) ---------------------

template <math::Scalar T>
T pv(const FixedCoupon& c, const CurveSet<T>& curves);
template <math::Scalar T>
T pv(const OvernightCoupon& c, const CurveSet<T>& curves, const FixingStore& fixings);
template <math::Scalar T>
T pv(const TermCoupon& c, const CurveSet<T>& curves, const FixingStore& fixings);
template <math::Scalar T>
T pv(const NotionalExchange& c, const CurveSet<T>& curves);
/// Dispatches on the variant; zero for payment ≤ reference date.
template <math::Scalar T>
T pv(const Cashflow& cf, const CurveSet<T>& curves, const FixingStore& fixings);
/// sign() × Σ pv(cashflow).
template <math::Scalar T>
T pv(const Leg& leg, const CurveSet<T>& curves, const FixingStore& fixings);

// --- trade PV ------------------------------------------------------------------------------------------

template <math::Scalar T>
T pv(const Swap& s, const CurveSet<T>& curves, const FixingStore& fixings);
template <math::Scalar T>
T pv(const BasisSwap& b, const CurveSet<T>& curves, const FixingStore& fixings);

/// Rate implied for a deposit by its index's projection curve: (P(start)/P(end) − 1)/tau.
template <math::Scalar T>
T impliedDepositRate(const Deposit& d, const CurveSet<T>& curves);
/// PV of a deposit as a forward-rate agreement settled at the end: N·tau·(implied − rate)·P_disc(end).
/// SIMPLIFICATION: a real deposit exchanges principal; for curve building only the rate matters.
template <math::Scalar T>
T pv(const Deposit& d, const CurveSet<T>& curves);
/// Forward implied for a FRA: (P(start)/P(end) − 1)/tau off the projection curve (fixing ignored if past;
/// bootstrap instruments never fix before the reference date).
template <math::Scalar T>
T impliedFraRate(const Fra& f, const CurveSet<T>& curves);
/// N·tau·(implied − rate)·P_disc(end). SIMPLIFICATION: settlement at end, not the discounted-at-start FRA convention.
template <math::Scalar T>
T pv(const Fra& f, const CurveSet<T>& curves);
/// Compounded forward over the reference period off the projection curve: (P(start)/P(end) − 1)/tau.
template <math::Scalar T>
T impliedFutureForward(const Future& f, const CurveSet<T>& curves);
/// Futures rate = forward + convexity adjustment; price = 100·(1 − rate).
template <math::Scalar T>
T impliedFutureRate(const Future& f, const CurveSet<T>& curves);
template <math::Scalar T>
T impliedFuturePrice(const Future& f, const CurveSet<T>& curves);

} // namespace rl::linear
