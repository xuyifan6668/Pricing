// T = rl::aad::Real instantiations of the templated bootstrap residuals: checks the analytic
// Jacobian against AAD through residual<Real>.
#include "rl/aad/real.hpp"
#include "rl/bootstrap/detail/residuals.ipp"

RL_BOOTSTRAP_INSTANTIATE_RESIDUALS(rl::aad::Real)
