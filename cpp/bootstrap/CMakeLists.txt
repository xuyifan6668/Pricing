# rl_bootstrap: instrument set, pillars, problem, report, FD helpers; residuals,
# Newton, adjoint Jacobian, par-risk transform.
rl_add_module(
  NAME bootstrap
  DEPS rl::core rl::math rl::math_real rl::market rl::instruments rl::linear nlohmann_json::nlohmann_json
  SOURCES
    src/instrument_set.cpp
    src/pillars.cpp
    src/residuals.cpp
    src/instantiate_real.cpp
    src/global_bootstrap.cpp
    src/jacobian.cpp
    src/report.cpp
    src/bootstrap_problem.cpp
  TEST_SOURCES
    test/test_instrument_set.cpp
    test/test_pillars.cpp
    test/test_bootstrap.cpp
    test/test_jacobian.cpp
)
target_compile_definitions(rl_bootstrap PRIVATE RL_DATA_DIR="${RL_DATA_DIR}")
