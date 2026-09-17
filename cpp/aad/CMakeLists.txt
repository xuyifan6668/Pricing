# rl_aad: the Real type and its tape. Links core (for the error types) and Eigen (traits header) only.
rl_add_module(
  NAME aad
  DEPS rl::core Eigen3::Eigen
  SOURCES
    src/tape.cpp
    src/real.cpp
  TEST_SOURCES
    test/test_tape.cpp
    test/test_real.cpp
    test/test_eigen_traits.cpp
  TEST_DEPS rl::math rl::math_real
)

# rl_math_real: the T = Real instantiations of rl::math. Consumers that price under both scalars
# link this; pure-double consumers link rl::math alone.
add_library(rl_math_real STATIC src/instantiate_math.cpp)
add_library(rl::math_real ALIAS rl_math_real)
target_link_libraries(rl_math_real PUBLIC rl::math rl::aad)
set_target_properties(rl_math_real PROPERTIES FOLDER "modules")
