# rl_linear: leg PV, OIS compounding kernel, annuity, par rate, forwards, DV01.
# Templated on the scalar; instantiated for double (src/linear.cpp) and Real (src/instantiate_real.cpp).
rl_add_module(
  NAME linear
  DEPS rl::core rl::math rl::math_real rl::market rl::instruments
  SOURCES
    src/linear.cpp
    src/instantiate_real.cpp
  TEST_SOURCES
    test/test_leg_pv.cpp
    test/test_ois_kernel.cpp
    test/test_par_rate.cpp
    test/test_forward_dv01.cpp
)
