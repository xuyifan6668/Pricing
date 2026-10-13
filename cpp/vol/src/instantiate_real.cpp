// T = rl::aad::Real instantiations of the rl::vol templates (the module links rl::math_real). Vega by cube node
//  differentiates through sabrNormalVol, cubeParams and alphaFromAtmNormalVol under Real.
#include "rl/aad/real.hpp"
#include "rl/vol/detail/bachelier.ipp"
#include "rl/vol/detail/black.ipp"
#include "rl/vol/detail/cube.ipp"
#include "rl/vol/detail/sabr.ipp"

RL_VOL_INSTANTIATE_BACHELIER(rl::aad::Real)
RL_VOL_INSTANTIATE_BLACK(rl::aad::Real)
RL_VOL_INSTANTIATE_SABR(rl::aad::Real)
RL_VOL_INSTANTIATE_CUBE(rl::aad::Real)
