# rl_xva: counterparty credit and funding adjustments on top of the pricing models. A credit curve
# from CDS spreads, a Monte Carlo exposure engine for a netting set with collateral, CVA / DVA / FVA,
# the exposure of a Bermudan swaption, CVA sensitivities, and wrong-way risk.
rl_add_module(
  NAME xva
  DEPS rl::core rl::math rl::math_real rl::market rl::instruments rl::linear rl::vol rl::bootstrap rl::models rl::pricers
  SOURCES
    src/credit.cpp
    src/exposure_engine.cpp
    src/xva.cpp
    src/bermudan_exposure.cpp
    src/cva_risk.cpp
    src/wrong_way.cpp
  TEST_SOURCES
    test/test_credit.cpp
    test/test_exposure_engine.cpp
    test/test_xva.cpp
    test/test_bermudan_exposure.cpp
    test/test_cva_risk.cpp
    test/test_wrong_way.cpp
)
