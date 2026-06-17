# Filter test harness executable
add_executable(quotient_filter_tests
        FilterTestHarness.cpp)

target_link_libraries(quotient_filter_tests PRIVATE
        svf_numeric
        Catch2::Catch2)

target_compile_features(quotient_filter_tests PRIVATE cxx_std_23)

# Default chart-data output directory baked in at compile time
target_compile_definitions(quotient_filter_tests PRIVATE
        QUOTIENT_DEFAULT_CHART_DIR="${QUOTIENT_CHART_DATA_DIR}")

# Register the individual Catch2 test cases with CTest
catch_discover_tests(quotient_filter_tests)

if(Python3_Interpreter_FOUND)
    add_custom_target(charts
            COMMAND quotient_filter_tests --export-charts --chart-dir "${QUOTIENT_CHART_DATA_DIR}"
            COMMAND "${Python3_EXECUTABLE}" "${QUOTIENT_ROOT}/tests/charts/plot_svf.py"
                    --data "${QUOTIENT_CHART_DATA_DIR}"
                    --out  "${QUOTIENT_ROOT}/tests/charts"
            DEPENDS quotient_filter_tests
            WORKING_DIRECTORY "${QUOTIENT_ROOT}/tests/charts"
            COMMENT "Running SVF tests, exporting chart data, and rendering matplotlib charts"
            USES_TERMINAL
            VERBATIM)
else()
    message(STATUS "Python3 interpreter not found; 'charts' target unavailable.")
endif()
