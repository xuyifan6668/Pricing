# rl_math: double instantiations only. The Real instantiations live in rl_math_real (cpp/aad), so
# that pure-double consumers never link the tape (docs/DESIGN.md §1, "math | aad").
rl_add_module(
  NAME math
  DEPS rl::core Eigen3::Eigen
  SOURCES
    src/interpolation.cpp
    src/tridiagonal.cpp
    src/normal.cpp
    src/rootfind.cpp
    src/lm.cpp
    src/quadrature.cpp
    src/rng.cpp
    src/sobol.cpp
    src/brownian_bridge.cpp
  TEST_SOURCES
    test/test_scalar.cpp
    test/test_interpolation.cpp
    test/test_tridiagonal.cpp
    test/test_normal.cpp
    test/test_rootfind.cpp
    test/test_lm.cpp
    test/test_dual.cpp
    test/test_quadrature.cpp
    test/test_rng.cpp
    test/test_sobol.cpp
    test/test_brownian_bridge.cpp
  TEST_DEPS rl::math_real
)
target_compile_definitions(rl_math PRIVATE RL_DATA_DIR="${RL_DATA_DIR}")
