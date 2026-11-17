include(FetchContent)

FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG v1.15.2
)
# Honor the parent project's runtime library selection on Windows.
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

include(GoogleTest)

add_executable(ope_tests
    test_black_scholes.cpp
    test_monte_carlo.cpp
    test_mc_greeks.cpp
    test_implied_vol.cpp
    test_vol_surface.cpp
    test_svi.cpp
    test_ssvi.cpp
)
target_link_libraries(ope_tests PRIVATE ope ope_warnings GTest::gtest_main)

gtest_discover_tests(ope_tests)
