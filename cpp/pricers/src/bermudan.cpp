#include "rl/pricers/bermudan.hpp"

#include "rl/core/error.hpp"
#include "rl/models/lgm_swaption.hpp"

#include <cmath>

namespace rl::pricers {

BermudanProblem::BermudanProblem(models::Lgm1F model, std::vector<models::ExerciseRight> rights) : model_(std::move(model)), rights_(std::move(rights)) {
    for (std::size_t i = 0; i < rights_.size(); ++i) {
        RL_REQUIRE(i == 0 || rights_[i].time > rights_[i - 1].time, "BermudanProblem: exercise dates must be increasing");
        zetas_.push_back(model_.zeta(rights_[i].time));
        value_.push_back(terms(model_, zetas_.back(), rights_[i].flows));
        annuity_.push_back(terms(model_, zetas_.back(), rights_[i].annuity));
    }
}

std::vector<BermudanProblem::Term> BermudanProblem::terms(const models::Lgm1F& model, double zeta, const std::vector<models::BondFlow>& flows) {
    std::vector<Term> out;
    for (const models::BondFlow& f : flows) { // flows are sorted by time
        const double h = model.H(f.time);
        const double c = f.amount * f.df * std::exp(-0.5 * h * h * zeta);
        if (!out.empty() && std::fabs(out.back().h - h) < 1e-14) out.back().c += c;
        else out.push_back({h, c});
    }
    return out;
}

namespace {
template <class Terms>
double sum(const Terms& terms, double x) noexcept {
    double g = 0.0;
    for (const auto& t : terms) g += t.c * std::exp(-t.h * x);
    return g;
}
} // namespace

double BermudanProblem::exerciseValue(std::size_t i, double x) const noexcept { return sum(value_[i], x); }

double BermudanProblem::swapRate(std::size_t i, double x) const noexcept {
    return rights_[i].strike + (rights_[i].payer ? 1.0 : -1.0) * sum(value_[i], x) / sum(annuity_[i], x);
}

std::vector<double> BermudanProblem::europeanValues() const {
    std::vector<double> out;
    for (const auto& r : rights_) out.push_back(models::lgmSwaptionPrice(model_, r));
    return out;
}

CalibratedBermudan calibrateBermudan(const models::UsdMarket& market, const instruments::Swap& swap, double firstCallYears, int noticeDays, double kappa) {
    std::vector<models::ExerciseRight> rights;
    for (core::Date d : models::coterminalExerciseDates(swap, firstCallYears, noticeDays, "USD-SIFMA"))
        rights.push_back(models::makeExerciseRight(swap, market.curves, d));
    models::LgmCalibration cal = models::calibrateLgm(kappa, models::makeStrip(rights, market.cube));
    BermudanProblem problem(cal.model, std::move(rights));
    return CalibratedBermudan{std::move(problem), std::move(cal)};
}

CalibratedBermudan calibrateBermudanAtStrikes(const models::UsdMarket& market, const instruments::Swap& swap, double firstCallYears, int noticeDays, double kappa,
                                              const std::function<instruments::Swap(double)>& withRate, const std::vector<double>& strikes) {
    const std::vector<core::Date> dates = models::coterminalExerciseDates(swap, firstCallYears, noticeDays, "USD-SIFMA");
    RL_REQUIRE(strikes.size() == dates.size(), "calibrateBermudanAtStrikes: need one calibration strike per exercise date");
    std::vector<models::ExerciseRight> rights, calibrationRights;
    for (std::size_t i = 0; i < dates.size(); ++i) {
        rights.push_back(models::makeExerciseRight(swap, market.curves, dates[i]));
        calibrationRights.push_back(models::makeExerciseRight(withRate(strikes[i]), market.curves, dates[i]));
        RL_REQUIRE(std::fabs(calibrationRights.back().time - rights.back().time) < 1e-12, "calibrateBermudanAtStrikes: withRate changed the swap's dates");
    }
    // makeStrip reads each right's own strike off the smile, so the strip is at the requested strikes.
    models::LgmCalibration cal = models::calibrateLgm(kappa, models::makeStrip(calibrationRights, market.cube));
    BermudanProblem problem(cal.model, std::move(rights));
    return CalibratedBermudan{std::move(problem), std::move(cal)};
}

} // namespace rl::pricers
