// Explicit instantiations of every rl::math template for T = rl::aad::Real (target rl_math_real).
// Lives in the aad module so that rl_math stays free of the tape and
// downstream templated modules link rl::math_real to get both scalars.
//
// When a new template is added to math: add its .ipp include and its RL_MATH_INSTANTIATE_* line
// here, and the double line in the corresponding cpp/math/src/*.cpp. A missing line shows up as an
// unresolved external for the Real symbol — that is the intended failure mode.

#include "rl/aad/real.hpp"
#include "rl/math/detail/interpolation.ipp"
#include "rl/math/detail/normal.ipp"
#include "rl/math/detail/tridiagonal.ipp"

RL_MATH_INSTANTIATE_INTERPOLATION(rl::aad::Real)
RL_MATH_INSTANTIATE_NORMAL(rl::aad::Real)
RL_MATH_INSTANTIATE_TRIDIAGONAL(rl::aad::Real)
