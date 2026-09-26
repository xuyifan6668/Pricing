// T = rl::aad::Real instantiations of the market templates (the module links rl::math_real).
#include "rl/aad/real.hpp"
#include "rl/market/detail/snapshot.ipp"
#include "rl/market/detail/swaption_cube.ipp"
#include "rl/market/detail/yield_curve.ipp"

RL_MARKET_INSTANTIATE_CURVES(rl::aad::Real)
RL_MARKET_INSTANTIATE_SNAPSHOT(rl::aad::Real)
RL_MARKET_INSTANTIATE_CUBE(rl::aad::Real)
