#include "rl/bootstrap/detail/residuals.ipp"

#include "rl/core/error.hpp"
#include "rl/linear/annuity.hpp"

RL_BOOTSTRAP_INSTANTIATE_RESIDUALS(double)

namespace rl::bootstrap {

// Rate-form residuals are already in rate units, so the scale to basis points of quote is 1e4 for
// every instrument (a future's residual is a price difference / 100, i.e. also a rate). With PV-form
// residuals this is where 1e4 / annuity would go.
double residualScale(const BootstrapInstrument&, const market::CurveSet<double>&, const market::FixingStore&) { return 1e4; }

} // namespace rl::bootstrap
