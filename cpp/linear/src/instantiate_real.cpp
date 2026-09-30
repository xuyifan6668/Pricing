// T = rl::aad::Real instantiations of rl::linear (the module links rl::math_real, so the tape is
// available; pure-double consumers still pay nothing at run time — Real is only instantiated, not used).
#include "rl/aad/real.hpp"
#include "rl/linear/detail/linear.ipp"

RL_LINEAR_INSTANTIATE(rl::aad::Real)
