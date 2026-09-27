# rl_instruments: conventions, IMM helpers, the trade types and the builders. Not templated on the
# scalar: a trade holds no market value (docs/DESIGN.md, "Trade representation").
rl_add_module(
  NAME instruments
  DEPS rl::core
  SOURCES
    src/conventions.cpp
    src/trades.cpp
    src/builder.cpp
  TEST_SOURCES
    test/test_conventions.cpp
    test/test_cashflow.cpp
    test/test_builder.cpp
)
