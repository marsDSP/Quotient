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

# FFT numeric correctness tests + chart-data exporter
add_library(fft_numeric OBJECT
        FftNumericTests.cpp
        FftSimdTests.cpp
        FftTestSupport.h
        FftChartData.h)

target_include_directories(fft_numeric PUBLIC
        "${CMAKE_CURRENT_SOURCE_DIR}"             # FftChartData.h
        "${QUOTIENT_ROOT}/source"                 # math/FFT.h
        "${QUOTIENT_ROOT}/libs/xsimd/include")    # xsimd SIMD batches (enables haveXSIMD)

target_compile_features(fft_numeric PUBLIC cxx_std_23)
target_link_libraries(fft_numeric PUBLIC Catch2::Catch2)
