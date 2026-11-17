cmake_minimum_required(VERSION 3.14)

project(options_pricing_engine
        VERSION 0.1.0
        DESCRIPTION "European option pricing & Greeks: Black-Scholes closed form + Monte Carlo"
        LANGUAGES CXX)

# ---- Global settings -------------------------------------------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(OPE_BUILD_TESTS "Build the unit test suite" ON)

# Warnings scoped to our own targets (via this interface library) so third-party
# code fetched by CMake — e.g. GoogleTest — isn't held to it.
add_library(ope_warnings INTERFACE)
target_compile_options(ope_warnings INTERFACE -Wall -Wextra -Wpedantic)

# ---- Core library ----------------------------------------------------------
# Sources are added milestone by milestone. M1: Black-Scholes closed form.
# M2: plain Monte Carlo. M6: implied volatility.
# M7: implied-vol surface. M8: SVI smile fit.
add_library(ope
    src/Option.cpp
    src/BlackScholes.cpp
    src/MonteCarlo.cpp
    src/ImpliedVol.cpp
    src/VolSurface.cpp
    src/Svi.cpp
    src/Ssvi.cpp
)
target_include_directories(ope PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(ope PRIVATE ope_warnings)

# ---- CLI: price one option, print price + Greeks ---------------------------
add_executable(pricer src/main.cpp)
target_link_libraries(pricer PRIVATE ope ope_warnings)

# ---- Validation tables & benchmarks (build Release for meaningful numbers) --
add_executable(ope_bench src/bench_main.cpp)
target_link_libraries(ope_bench PRIVATE ope ope_warnings)

# ---- CSV export for the plots in docs/img (plotting itself is scripts/plot.py) --
add_executable(ope_export src/export_main.cpp)
target_link_libraries(ope_export PRIVATE ope ope_warnings)

# ---- Tests -----------------------------------------------------------------
if(OPE_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()
