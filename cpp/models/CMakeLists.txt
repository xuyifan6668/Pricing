# rl_models: the LGM one-factor and G2++ two-factor models, European swaption pricing under them,
# and their calibration to swaptions. Double only: a calibration is a solve, and risk on model prices is
# taken by bumping the market and re-calibrating.
rl_add_module(
  NAME models
  DEPS rl::core rl::math rl::market rl::instruments rl::linear rl::vol rl::bootstrap
  SOURCES
    src/lgm.cpp
    src/flows.cpp
    src/lgm_swaption.cpp
    src/lgm_calibration.cpp
    src/usd_market.cpp
    src/g2pp.cpp
    src/g2pp_swaption.cpp
    src/g2pp_calibration.cpp
  TEST_SOURCES
    test/test_lgm.cpp
    test/test_flows.cpp
    test/test_lgm_swaption.cpp
    test/test_lgm_calibration.cpp
    test/test_g2pp.cpp
    test/test_g2pp_calibration.cpp
)
