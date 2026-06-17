# ============================================================================
#  Quotient filter tests
# ============================================================================
cmake_minimum_required(VERSION 3.23.1)

# Standalone configuration provides its own project()
if(NOT DEFINED PROJECT_NAME)
    project(QuotientFilterTests CXX)
endif()

get_filename_component(QUOTIENT_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/.." ABSOLUTE)

if(NOT DEFINED CMAKE_CXX_STANDARD)
    set(CMAKE_CXX_STANDARD 23)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    set(CMAKE_CXX_EXTENSIONS OFF)
endif()

# Python is optional
find_package(Python3 COMPONENTS Interpreter)

set(QUOTIENT_CHART_DATA_DIR "${QUOTIENT_ROOT}/tests/charts/data"
        CACHE PATH "Directory for generated SVF chart-data CSVs")

enable_testing()

# Vendored Catch2 v3. gcem's cmake_minimum_required floor is also relaxed
add_subdirectory("${QUOTIENT_ROOT}/libs/Catch2" "${CMAKE_CURRENT_BINARY_DIR}/Catch2")

# Catch2 appends its CMake helpers
include(Catch)

add_subdirectory(numeric)
add_subdirectory(harnesses)
