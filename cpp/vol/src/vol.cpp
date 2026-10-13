// T = double instantiations of the rl::vol templates.
#include "rl/vol/detail/bachelier.ipp"
#include "rl/vol/detail/black.ipp"
#include "rl/vol/detail/cube.ipp"
#include "rl/vol/detail/sabr.ipp"

RL_VOL_INSTANTIATE_BACHELIER(double)
RL_VOL_INSTANTIATE_BLACK(double)
RL_VOL_INSTANTIATE_SABR(double)
RL_VOL_INSTANTIATE_CUBE(double)
