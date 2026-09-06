rl_add_module(
  NAME core
  SOURCES
    src/error.cpp
    src/date.cpp
    src/tenor.cpp
    src/daycount.cpp
    src/calendar.cpp
    src/currency.cpp
    src/index.cpp
    src/schedule.cpp
  TEST_SOURCES
    test/test_error.cpp
    test/test_date.cpp
    test/test_tenor.cpp
    test/test_daycount.cpp
    test/test_calendar.cpp
    test/test_index.cpp
    test/test_schedule.cpp
)
target_compile_definitions(rl_core PRIVATE RL_DATA_DIR="${RL_DATA_DIR}")
