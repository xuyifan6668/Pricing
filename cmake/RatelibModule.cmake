# rl_add_module(NAME <name> [DEPS <targets...>] [SOURCES <files...>] [TEST_SOURCES <files...>]
#               [TEST_DEPS <targets...>] [HEADER_ONLY])
#
# Creates:
#   rl_<name>        STATIC library (or INTERFACE if HEADER_ONLY), alias rl::<name>,
#                    public include dir cpp/<name>/include, linked to rl::flags and DEPS.
#   rl_<name>_tests  Catch2 executable from cpp/<name>/test/*.cpp, registered with CTest.
#
# The module dependency order (docs/DESIGN.md) is encoded here and nowhere else.

function(rl_add_module)
  set(options HEADER_ONLY)
  set(oneValueArgs NAME)
  set(multiValueArgs DEPS SOURCES TEST_SOURCES TEST_DEPS)
  cmake_parse_arguments(RL "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  set(target rl_${RL_NAME})
  set(inc "${CMAKE_CURRENT_SOURCE_DIR}/include")

  if(RL_HEADER_ONLY)
    add_library(${target} INTERFACE)
    target_include_directories(${target} INTERFACE ${inc})
    target_link_libraries(${target} INTERFACE rl::flags ${RL_DEPS})
  else()
    add_library(${target} STATIC ${RL_SOURCES})
    target_include_directories(${target} PUBLIC ${inc})
    target_link_libraries(${target} PUBLIC rl::flags ${RL_DEPS})
    set_target_properties(${target} PROPERTIES FOLDER "modules")
  endif()
  add_library(rl::${RL_NAME} ALIAS ${target})

  if(RL_BUILD_TESTS AND RL_TEST_SOURCES)
    set(test_target ${target}_tests)
    add_executable(${test_target} ${RL_TEST_SOURCES})
    target_link_libraries(${test_target} PRIVATE ${target} ${RL_TEST_DEPS} Catch2::Catch2WithMain)
    target_compile_definitions(${test_target} PRIVATE RL_TEST_DATA_DIR="${RL_DATA_DIR}")
    set_target_properties(${test_target} PROPERTIES FOLDER "tests")
    catch_discover_tests(${test_target}
      TEST_PREFIX "${RL_NAME}."
      PROPERTIES ENVIRONMENT "RL_DATA_DIR=${RL_DATA_DIR}")
  endif()
endfunction()
