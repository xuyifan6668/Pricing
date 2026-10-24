# rl_pricers: the Bermudan swaption priced three ways under LGM (PDE, least-squares Monte Carlo lower
# bound, Andersen–Broadie upper bound), swap exposure profiles, and bump / adjoint risk.
rl_add_module(
  NAME pricers
  DEPS rl::core rl::math rl::math_real rl::market rl::instruments rl::linear rl::vol rl::bootstrap rl::models
  SOURCES
    src/bermudan.cpp
    src/bermudan_pde.cpp
    src/bermudan_mc.cpp
    src/exposure.cpp
    src/risk.cpp
  TEST_SOURCES
    test/test_bermudan_pde.cpp
    test/test_bermudan_mc.cpp
    test/test_exposure.cpp
    test/test_risk.cpp
)
