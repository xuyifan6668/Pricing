#pragma once
// The counterparty's credit curve: survival probabilities from CDS spreads.
//
// CVA needs the probability that the counterparty defaults in each future interval. The market
// quotes that through credit default swaps: the protection buyer pays a running spread until
// default or maturity, and receives (1 − recovery) of notional at default. A hazard rate λ(t) that
// is constant between CDS maturities is bootstrapped, shortest maturity first, so that each quoted
// spread makes its CDS worth zero:
//
//     premium leg     s · Σ_j [ Δ_j D(t_j) Q(t_j)  +  ½ Δ_j D(t_j) (Q(t_{j−1}) − Q(t_j)) ]
//     protection leg  (1 − R) · Σ_j D(t_j) (Q(t_{j−1}) − Q(t_j))
//
// with quarterly dates t_j, Q(t) = exp(−∫₀ᵗ λ) the survival probability and D the risk-free
// discount factor. The second premium term is the spread accrued up to a default inside the period.
//
// SIMPLIFICATION: a quarterly grid in year fractions (no IMM dates, no upfront plus standard
// coupon), default paid at the end of its quarter. Rates and credit are independent. To first order
// λ ≈ s / (1 − R), the "credit triangle", which the tests check.

#include <cstddef>
#include <functional>
#include <vector>

namespace rl::xva {

struct CdsQuote {
    double maturityYears = 0.0;
    double spread = 0.0; ///< running par spread, decimal (0.0120 = 120 bp)
};

class HazardCurve {
public:
    /// Hazard rate hazards[i] on (times[i−1], times[i]], times[−1] = 0; the last continues beyond.
    /// Throws InputError unless times are positive and strictly increasing and hazards non-negative.
    HazardCurve(std::vector<double> times, std::vector<double> hazards);

    /// A constant hazard rate.
    static HazardCurve flat(double hazard);

    /// Bootstrap to par CDS spreads in increasing maturity. `discount(t)` is the risk-free discount
    /// factor. Throws InputError on unordered quotes or a recovery outside [0, 1); CalibrationError
    /// if a spread implies a negative hazard rate (the curve is arbitrageable).
    static HazardCurve bootstrap(const std::vector<CdsQuote>& quotes, double recovery, const std::function<double(double)>& discount);

    const std::vector<double>& times() const noexcept { return times_; }
    const std::vector<double>& hazards() const noexcept { return hazards_; }

    /// Q(t) = P(no default before t).
    double survival(double t) const;
    /// P(default in (t0, t1]).
    double defaultProbability(double t0, double t1) const { return survival(t0) - survival(t1); }

private:
    std::vector<double> times_, hazards_;
};

/// Par spread of a CDS of the given maturity on `curve`: the spread that makes it worth zero.
double cdsParSpread(const HazardCurve& curve, double maturityYears, double recovery, const std::function<double(double)>& discount);

} // namespace rl::xva
