# SVF numeric correctness tests + chart-data exporter
add_library(svf_numeric OBJECT
        SvfNumericTests.cpp
        SvfReference.h
        SvfChartData.h)

target_include_directories(svf_numeric PUBLIC
        "${CMAKE_CURRENT_SOURCE_DIR}"             # SvfReference.h / SvfChartData.h
        "${QUOTIENT_ROOT}/source"                 # dsp/filters/StateVariable.h
        "${QUOTIENT_ROOT}/libs/gcem/include")     # gcem.hpp

target_compile_features(svf_numeric PUBLIC cxx_std_23)
target_link_libraries(svf_numeric PUBLIC Catch2::Catch2)
