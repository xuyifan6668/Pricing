cmake_minimum_required(VERSION 3.25)

project(ratelib
  VERSION 0.1.0
  DESCRIPTION "Exotic rates pricing library: multi-curve, SABR, LGM, Bermudan (PDE/LSMC/AB), AAD"
  LANGUAGES CXX)

# ---------------------------------------------------------------------------------------------------
# Global settings
# ---------------------------------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

option(RL_BUILD_TESTS    "Build the Catch2 test executables"        ON)
option(RL_BUILD_EXAMPLES "Build the demo executables used by notebooks" ON)
option(RL_SANITIZE       "Enable ASan (+UBSan on clang/gcc)"         OFF)

# Data directory baked in at configure time; runtime override via env RL_DATA_DIR (see core/data_dir).
set(RL_DATA_DIR "${CMAKE_CURRENT_SOURCE_DIR}/data" CACHE PATH "Directory holding holidays/ and sobol/")

include(CompilerFlags)
include(RatelibModule)

find_package(Eigen3 CONFIG REQUIRED)
find_package(nlohmann_json CONFIG REQUIRED)
if(RL_BUILD_TESTS)
  find_package(Catch2 3 CONFIG REQUIRED)
  include(CTest)
  include(Catch)
endif()

# ---------------------------------------------------------------------------------------------------
# Modules — the dependency DAG of docs/DESIGN.md §1 is encoded ONLY through the DEPS lists.
# A cyclic include fails to link instead of silently accumulating.
# ---------------------------------------------------------------------------------------------------
add_subdirectory(cpp/core)
add_subdirectory(cpp/math)
add_subdirectory(cpp/aad)
add_subdirectory(cpp/market)        # quotes, fixings, hash, curves, snapshot
add_subdirectory(cpp/instruments)   # conventions, trade types, builders
add_subdirectory(cpp/linear)        # leg PV, OIS compounding kernel, annuity, par rate, DV01
add_subdirectory(cpp/bootstrap)     # global multi-curve solve, adjoint Jacobian, par risk
add_subdirectory(cpp/vol)           # Bachelier/Black, SABR and its fit, cube queries, arbitrage check
add_subdirectory(cpp/models)        # LGM 1F, swaption pricing under it, strip calibration, the shared USD market
add_subdirectory(cpp/pricers)       # Bermudan by PDE, LSMC and Andersen-Broadie; swap exposure; bump and adjoint risk
add_subdirectory(cpp/xva)           # credit curve, exposure engine, CVA/DVA/FVA, Bermudan exposure, CVA risk

if(RL_BUILD_EXAMPLES)
  add_subdirectory(cpp/examples)
endif()
