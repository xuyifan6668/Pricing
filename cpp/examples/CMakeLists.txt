# Demo executables the notebooks call through subprocess. Each prints CSV
# to stdout and takes its parameters on the command line; nothing here is library code.
foreach(demo rl_demo_schedule rl_demo_interp rl_demo_sobol)
  add_executable(${demo} ${demo}.cpp)
  target_link_libraries(${demo} PRIVATE rl::core rl::math)
  set_target_properties(${demo} PROPERTIES FOLDER "examples")
endforeach()

# demo (notebooks 04 and 05).
add_executable(rl_demo_bootstrap rl_demo_bootstrap.cpp)
target_link_libraries(rl_demo_bootstrap PRIVATE rl::core rl::math rl::market rl::instruments rl::linear rl::bootstrap)
set_target_properties(rl_demo_bootstrap PROPERTIES FOLDER "examples")

# demo (notebooks 06 and 07).
add_executable(rl_demo_vol rl_demo_vol.cpp)
target_link_libraries(rl_demo_vol PRIVATE rl::core rl::math rl::market rl::vol)
set_target_properties(rl_demo_vol PROPERTIES FOLDER "examples")

# Every number and plot in docs/RESULTS.md comes from this one executable (tools/make_results.py).
add_executable(rl_analysis rl_analysis.cpp)
target_link_libraries(rl_analysis PRIVATE rl::pricers)
set_target_properties(rl_analysis PROPERTIES FOLDER "examples")
