# rl_market: quotes, fixings, hashing, the curve / curve-set / snapshot templates, and the vol
# quotes and the swaption cube container. Templates are
# declared in include/, defined in include/rl/market/detail/*.ipp, instantiated for double in
# src/*.cpp and for rl::aad::Real in src/instantiate_real.cpp.
rl_add_module(
  NAME market
  DEPS rl::core rl::math rl::math_real nlohmann_json::nlohmann_json
  SOURCES
    src/curve_interpolation.cpp
    src/quote.cpp
    src/fixings.cpp
    src/hash.cpp
    src/yield_curve.cpp
    src/snapshot.cpp
    src/swaption_cube.cpp
    src/vol_quotes.cpp
    src/instantiate_real.cpp
  TEST_SOURCES
    test/test_curve_interpolation.cpp
    test/test_quote.cpp
    test/test_fixings.cpp
    test/test_hash.cpp
    test/test_yield_curve.cpp
    test/test_curve_set.cpp
    test/test_snapshot.cpp
    test/test_swaption_cube.cpp
    test/test_vol_quotes.cpp
)
target_compile_definitions(rl_market PRIVATE RL_DATA_DIR="${RL_DATA_DIR}")
