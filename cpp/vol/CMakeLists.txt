# rl_vol: Bachelier / shifted Black and implied vols, Hagan SABR and its gradient, cube queries, the
# arbitrage diagnostic, vol conversion and the SABR calibration.
# Templated on the scalar; instantiated for double (src/vol.cpp) and Real (src/instantiate_real.cpp).
# DAG (docs/DESIGN.md section 1): math -> vol, market -> vol. vol does NOT link linear.
rl_add_module(
  NAME vol
  DEPS rl::core rl::math rl::math_real rl::market
  SOURCES
    src/vol.cpp
    src/instantiate_real.cpp
    src/implied_vol.cpp
    src/conversion.cpp
    src/sabr_fit.cpp
    src/arbitrage.cpp
  TEST_SOURCES
    test/test_bachelier.cpp
    test/test_black.cpp
    test/test_conversion.cpp
    test/test_sabr.cpp
    test/test_sabr_fit.cpp
    test/test_cube.cpp
    test/test_arbitrage.cpp
)
