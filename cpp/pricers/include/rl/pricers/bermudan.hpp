#pragma once
// The Bermudan swaption as the pricers see it.
//
// The holder may, on any one of the exercise dates t_1 < … < t_n, enter the remaining part of a
// fixed swap. Exercising gives up the right to exercise later, so at each date the holder compares
//     exercise value   h_i(x) = value of the remaining swap in state x
//     continuation     C_i(x) = value of still holding the option
// and the option is worth max(h_i, C_i) there. Three facts follow without any model and are used as
// tests of every pricer:
//     Bermudan ≥ every European with one of its exercise dates  (never exercising elsewhere is allowed)
//     Bermudan is non-decreasing in the set of exercise dates   (more rights cannot hurt)
//     Bermudan with a single exercise date IS that European
//
// Everything is in DEFLATED units, V/N: under LGM these are martingales in a driftless Gaussian state,
// so no discounting appears anywhere in the pricers; the price today is the deflated value at x = 0.

#include "rl/instruments/swap.hpp"
#include "rl/models/flows.hpp"
#include "rl/models/lgm.hpp"
#include "rl/models/lgm_calibration.hpp"
#include "rl/models/usd_market.hpp"

#include <cstddef>
#include <vector>

namespace rl::pricers {

class BermudanProblem {
public:
    /// Rights must be in increasing exercise time. Throws InputError otherwise.
    BermudanProblem(models::Lgm1F model, std::vector<models::ExerciseRight> rights);

    const models::Lgm1F& model() const noexcept { return model_; }
    const std::vector<models::ExerciseRight>& rights() const noexcept { return rights_; }
    std::size_t size() const noexcept { return rights_.size(); }
    /// ζ at exercise date i: the variance of the state there.
    double zeta(std::size_t i) const noexcept { return zetas_[i]; }
    /// Deflated value of exercising at date i in state x: Σ_k c_k exp(−H_k x).
    double exerciseValue(std::size_t i, double x) const noexcept;
    /// The par rate of the swap entered at date i, in state x (deflation cancels in the ratio).
    double swapRate(std::size_t i, double x) const noexcept;
    /// Model prices of the Europeans with each single exercise date.
    std::vector<double> europeanValues() const;

private:
    // DESIGN: per exercise date, every bond's (H_k, a_k D_k exp(−½H_k²ζ)) is computed once here, with
    // bonds on the same date merged. The pricers call exerciseValue millions of times (every grid
    // node, every path at every date); this leaves one exp per distinct cashflow date in that loop.
    struct Term {
        double h;
        double c;
    };
    static std::vector<Term> terms(const models::Lgm1F& model, double zeta, const std::vector<models::BondFlow>& flows);

    models::Lgm1F model_;
    std::vector<models::ExerciseRight> rights_;
    std::vector<double> zetas_;
    std::vector<std::vector<Term>> value_, annuity_;
};

/// A Bermudan with its model calibrated to its own co-terminal strip.
struct CalibratedBermudan {
    BermudanProblem problem;
    models::LgmCalibration calibration;
};

/// Exercise rights from `swap` (annual call dates after `firstCallYears`, `noticeDays` business days
/// of notice), strip from the market's cube at the swap's fixed rate, LGM bootstrapped at `kappa`.
CalibratedBermudan calibrateBermudan(const models::UsdMarket& market, const instruments::Swap& swap, double firstCallYears, int noticeDays, double kappa);

} // namespace rl::pricers
