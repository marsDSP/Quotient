#pragma once

#include <string>

namespace quotient_fft_charts {

/**
 * Write FFT chart CSVs to outDir.
 *
 * Produces:
 *   fft_spectrum.csv    – magnitude spectrum of a pure tone (N=512, bin 40)
 *   fft_accuracy.csv    – max|FFT(x) - DFT(x)| and round-trip error vs N
 *   fft_parseval.csv    – relative Parseval energy error vs N
 *
 * Returns true on success, false if any file could not be written.
 */
bool exportChartData(const std::string &outDir);

/**
 * Write SIMD chart CSVs to outDir.
 *
 * Produces:
 *   fft_simd_accuracy.csv   – SoA SIMD error vs DFT and vs scalar FFT, vs N
 *   fft_phase_coherence.csv – per-bin phase/magnitude error of a swept tone
 *   fft_throughput.csv      – SoA SIMD vs scalar throughput (Msamples/s), vs N
 *
 * Returns true on success, false if any file could not be written.
 */
bool exportSimdChartData(const std::string &outDir);

} // namespace quotient_fft_charts
