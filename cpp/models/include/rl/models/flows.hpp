#pragma once
// The model's view of a swap: a list of zero-coupon bond positions.
//
// LGM prices anything that can be written as Σ_k a_k · P(t, T_k), because each deflated bond is one
// exponential in the state. A swap is exactly that:
//   fixed coupon   pays N K τ at its payment date                 → one bond, amount ∓N K τ
//   float coupon   pays N (P(s)/P(e) − 1) for accrual [s, e]      → +N at s, −N at e
// (receiving compounded overnight interest over [s, e] on N is the same as being given N at s and
// handing N back at e: invest at the overnight rate in between).
//
// SIMPLIFICATION: OIS coupons pay two business days after the accrual end. That delay is kept as a
// deterministic factor δ = D(pay)/D(e) on the float coupon, frozen at today's curve, so that the bond
// decomposition reprices the swap on today's curve EXACTLY while staying a pure sum of bonds. The
// neglected effect is the convexity of a two-day delay, far below a hundredth of a basis point.
// Single-curve: the float index must be projected off the discount curve (SOFR swaps under SOFR
// discounting). A Fed Funds swap would need a deterministic basis spread on top.

#include "rl/core/date.hpp"
#include "rl/instruments/bermudan.hpp"
#include "rl/instruments/swap.hpp"
#include "rl/market/curve_set.hpp"

#include <string>
#include <vector>

namespace rl::models {

struct BondFlow {
    double time = 0.0;   ///< T_k, years ACT/365F from the curve's reference date
    double amount = 0.0; ///< a_k, signed from the holder's side (positive = received)
    double df = 0.0;     ///< D(T_k) = P(0, T_k) on the discount curve
};

/// What the holder gets by exercising at `time`: the swap coupons that start accruing on or after
/// the exercise date, as bonds, plus the same swap's fixed-leg annuity as bonds (for reading a par
/// rate off a model value: S = K ± V/A).
struct ExerciseRight {
    core::Date date;
    double time = 0.0;             ///< exercise time in years
    std::vector<BondFlow> flows;   ///< the swap entered into, from the holder's side
    std::vector<BondFlow> annuity; ///< N τ_i at each fixed payment date (positive)
    double strike = 0.0;           ///< the swap's fixed rate
    bool payer = true;             ///< holder pays fixed on exercise
    double tenorYears = 0.0;       ///< length of the swap entered into, for the vol cube
    std::string label;

    /// Value today of entering the swap now at the fixed strike: Σ a_k D(T_k).
    double forwardValue() const noexcept;
    /// Σ N τ_i D(pay_i).
    double annuityValue() const noexcept;
    /// Forward par rate of the swap entered into: K ± forwardValue/annuityValue.
    double forwardRate() const noexcept;
};

/// The coupons of `swap` whose accrual starts on or after `from`, as bond flows from the side of
/// whoever holds the swap as written. Throws InputError if a float coupon is not an overnight coupon
/// on the discount curve's index, has a spread, or has already started accruing.
std::vector<BondFlow> swapAsBonds(const instruments::Swap& swap, const market::CurveSet<double>& curves, core::Date from);

/// The exercise right into `swap` at `exerciseDate`.
ExerciseRight makeExerciseRight(const instruments::Swap& swap, const market::CurveSet<double>& curves, core::Date exerciseDate);

/// One exercise right per exercise date of the Bermudan, in date order.
std::vector<ExerciseRight> makeExerciseRights(const instruments::BermudanSwaption& bermudan, const market::CurveSet<double>& curves);

/// Exercise dates for a co-terminal Bermudan on `swap`: `noticeDays` business days before each fixed
/// period start that is at least `firstCallYears` after the swap's own start (the "non-call" period).
std::vector<core::Date> coterminalExerciseDates(const instruments::Swap& swap, double firstCallYears, int noticeDays, const std::string& calendar);

} // namespace rl::models
