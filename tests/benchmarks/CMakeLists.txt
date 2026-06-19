# ============================================================================
#  FFT micro-benchmark (Catch2 benchmarking)
#  Standalone executable; intentionally NOT registered with CTest so the timing
#  runs do not slow the normal ctest cycle. Build/run in Release:
#     cmake -S tests -B build-bench -DCMAKE_BUILD_TYPE=Release
#     cmake --build build-bench --target quotient_fft_bench
#     ./build-bench/benchmarks/quotient_fft_bench "[bench]"
# ============================================================================
add_executable(quotient_fft_bench FftBenchmark.cpp)

target_include_directories(quotient_fft_bench PRIVATE
        "${QUOTIENT_ROOT}/source")          # math/FFT.h

target_compile_features(quotient_fft_bench PRIVATE cxx_std_23)
target_link_libraries(quotient_fft_bench PRIVATE Catch2::Catch2WithMain)
