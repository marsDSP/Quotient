// ============================================================================
//  FftSimdTests.cpp
//  Correctness, phase-coherence, and throughput tests for the
//  SoA/SIMD engine, plus the RealFFT wrapper that sits on top of it.
//  Also defines the SIMD chart-data exporter.
// ============================================================================
#include <catch2/catch_test_macros.hpp>

#include "FftTestSupport.h"
#include "FftChartData.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <numbers>
#include <string>
#include <vector>

using namespace ffttest;

namespace
{
    constexpr double kPi = std::numbers::pi_v<double>;

    // Wrap an angle (radians) into (-pi, pi].
    double wrapPi(double a)
    {
        while (a > kPi) a -= kTwoPi;
        while (a <= -kPi) a += kTwoPi;
        return a;
    }

    // Sizes that guarantee the SIMD kernels run on common ISAs (NEON W=2,
    // AVX W=4, AVX-512 W=8): once a stage has Ns >= W the batch path is taken.
    const std::vector<std::size_t> kSizes = {
        8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536
    };
}

// ----------------------------------------------------------------------------
//  Guard: the tests must actually be exercising SIMD.
// ----------------------------------------------------------------------------
TEST_CASE("SoA engine is compiled with SIMD enabled", "[fft][simd]")
{
    REQUIRE(MarsDSP::MathOps::haveXSIMD);
}

// ----------------------------------------------------------------------------
//  Correctness: SIMD output vs ground-truth DFT and vs the scalar engine.
// ----------------------------------------------------------------------------
TEST_CASE("SoA SIMD forward matches the reference DFT", "[fft][simd][correctness]")
{
    for (const std::size_t N : {2uz, 4uz, 8uz, 16uz, 32uz, 64uz, 128uz, 256uz, 512uz})
    {
        const auto x = makeRandom(N, 1234u + N);
        const double err = maxAbsDiff(soaForward(x), dftRef(x));
        CAPTURE(N, err);
        REQUIRE(err < 1e-9 * static_cast<double>(N));
    }
}

TEST_CASE("SoA SIMD matches the scalar recursive FFT across sizes", "[fft][simd][correctness]")
{
    for (const std::size_t N : kSizes)
    {
        const auto x = makeRandom(N, 7u + N);
        const double fe = maxAbsDiff(soaForward(x), scalarForward(x));
        CAPTURE(N, fe);
        REQUIRE(fe < 1e-9 * static_cast<double>(N));
    }
}

TEST_CASE("SoA SIMD inverse round-trips to the input", "[fft][simd][roundtrip]")
{
    for (const std::size_t N : kSizes)
    {
        const auto x = makeRandom(N, 99u + N);
        auto r = soaInverse(soaForward(x));
        for (auto &v : r) v /= static_cast<double>(N);
        const double err = maxAbsDiff(r, x);
        CAPTURE(N, err);
        REQUIRE(err < 1e-9);
    }
}

TEST_CASE("SoA SIMD unit impulse yields a flat spectrum", "[fft][simd]")
{
    for (const std::size_t N : {8uz, 64uz, 256uz, 1024uz, 4096uz})
    {
        std::vector<Complex> impulse(N);
        impulse[0] = Complex{1.0, 0.0};
        const auto X = soaForward(impulse);
        for (std::size_t k = 0; k < N; ++k)
        {
            CAPTURE(N, k);
            REQUIRE(std::abs(X[k] - Complex{1.0, 0.0}) < 1e-9);
        }
    }
}

// ----------------------------------------------------------------------------
//  Phase coherence / no detuning.
// ----------------------------------------------------------------------------
TEST_CASE("SoA SIMD places tones on the exact bin with the exact phase", "[fft][simd][phase]")
{
    const std::size_t N = 4096;
    for (const double phase : {0.0, 0.3, 1.1, -2.0, 3.0})
    {
        for (const std::size_t k0 : {1uz, 5uz, 64uz, 333uz, 2047uz})
        {
            const auto x = makeRealTone(N, static_cast<double>(k0), 1.0, phase);
            const auto X = soaForward(x);

            std::size_t peak = 1;
            double peakMag = 0.0;
            for (std::size_t k = 1; k < N / 2; ++k)
            {
                const double m = std::abs(X[k]);
                if (m > peakMag) { peakMag = m; peak = k; }
            }
            CAPTURE(N, k0, phase, peak);
            REQUIRE(peak == k0);                                          // no bin shift
            REQUIRE(std::abs(std::abs(X[k0]) - 0.5 * N) < 1e-6 * N);      // correct magnitude
            REQUIRE(std::abs(wrapPi(std::arg(X[k0]) - phase)) < 1e-6);    // correct phase
        }
    }
}

TEST_CASE("SoA SIMD introduces no detuning vs scalar (off-bin tones + noise)", "[fft][simd][phase]")
{
    const std::size_t N = 8192;
    for (const double f : {10.5, 100.25, 1000.75, 3333.3})
    {
        const auto x = makeRealTone(N, f, 1.0, 0.7);
        const double e = maxAbsDiff(soaForward(x), scalarForward(x));
        CAPTURE(N, f, e);
        REQUIRE(e < 1e-9 * static_cast<double>(N));
    }
    for (std::uint64_t s = 0; s < 8; ++s)
    {
        const auto x = makeRandom(N, 4000u + s);
        const double e = maxAbsDiff(soaForward(x), scalarForward(x));
        CAPTURE(N, s, e);
        REQUIRE(e < 1e-9 * static_cast<double>(N));
    }
}

TEST_CASE("SoA SIMD preserves every partial's phase in a multitone", "[fft][simd][phase]")
{
    const std::size_t N = 2048;
    struct Tone { std::size_t k; double amp; double phase; };
    const std::vector<Tone> tones = {
        {16, 1.0, 0.2}, {37, 0.7, -1.0}, {200, 0.5, 2.5}, {512, 0.9, 0.0}
    };

    std::vector<Complex> x(N);
    for (const auto &t : tones)
    {
        const auto s = makeRealTone(N, static_cast<double>(t.k), t.amp, t.phase);
        for (std::size_t i = 0; i < N; ++i) x[i] += s[i];
    }

    const auto X = soaForward(x);
    for (const auto &t : tones)
    {
        CAPTURE(N, t.k, t.amp, t.phase);
        REQUIRE(std::abs(std::abs(X[t.k]) - 0.5 * t.amp * N) < 1e-6 * N);
        REQUIRE(std::abs(wrapPi(std::arg(X[t.k]) - t.phase)) < 1e-6);
    }
    REQUIRE(maxAbsDiff(X, scalarForward(x)) < 1e-9 * static_cast<double>(N));
}

// ----------------------------------------------------------------------------
//  RealFFT (real-input wrapper over the SIMD engine).
// ----------------------------------------------------------------------------
TEST_CASE("RealFFT forward matches the DFT and round-trips", "[fft][real][simd]")
{
    for (const std::size_t n : {8uz, 16uz, 64uz, 256uz, 1024uz})
    {
        std::mt19937_64 rng(1234u + n);
        std::uniform_real_distribution<double> dist(-1.0, 1.0);
        std::vector<double> in(n), out(n);
        for (auto &v : in) v = dist(rng);

        MarsDSP::MathOps::RealFFT<double> rfft(n);
        const std::size_t half = n / 2;
        std::vector<Complex> spec(half);
        rfft.forward(in.data(), spec.data());

        // Compare unique bins against the full complex DFT of the real signal.
        std::vector<Complex> xc(n);
        for (std::size_t i = 0; i < n; ++i) xc[i] = Complex{in[i], 0.0};
        const auto full = dftRef(xc);
        for (std::size_t k = 1; k < half; ++k)
        {
            CAPTURE(n, k);
            REQUIRE(std::abs(spec[k] - full[k]) < 1e-9 * static_cast<double>(n));
        }
        // out[0] packs DC (real) and Nyquist (imag).
        REQUIRE(std::abs(spec[0].real() - full[0].real()) < 1e-9 * static_cast<double>(n));
        REQUIRE(std::abs(spec[0].imag() - full[half].real()) < 1e-9 * static_cast<double>(n));

        // The inverse is unnormalised, matching the complex engine convention
        // (inverse(forward(x)) == n*x), so divide by n to recover the input.
        rfft.inverse(spec.data(), out.data());
        double err = 0.0;
        for (std::size_t i = 0; i < n; ++i)
            err = std::max(err, std::abs(out[i] / static_cast<double>(n) - in[i]));
        CAPTURE(n, err);
        REQUIRE(err < 1e-9);
    }
}

// ----------------------------------------------------------------------------
//  Throughput is measurable (no hard perf threshold — CI timing is noisy).
// ----------------------------------------------------------------------------
TEST_CASE("SoA SIMD throughput is measurable", "[fft][simd][throughput]")
{
    const std::size_t N = 4096;
    const auto x = makeRandom(N, 5u);
    SoaFFT e;
    e.setSize(N);
    std::vector<Complex> out(N);

    const double tns = medianTimeNs([&]
    {
        e.fwd(x.data(), out.data());
        return out[0].real() + out[N / 2].imag();
    });
    const double msps = (static_cast<double>(N) / (tns * 1e-9)) / 1e6;
    CAPTURE(N, tns, msps);
    REQUIRE(tns > 0.0);
    REQUIRE(msps > 0.0);
}

// ============================================================================
//  SIMD chart-data export
// ============================================================================
bool quotient_fft_charts::exportSimdChartData(const std::string &outDir)
{
    namespace fs = std::filesystem;
    std::error_code ec;
    fs::create_directories(outDir, ec);
    if (ec) return false;

    const fs::path dir(outDir);
    bool ok = true;

    // (a) SIMD accuracy vs N: SoA SIMD compared to DFT and to the scalar engine.
    {
        const std::vector<std::size_t> sizes = {16, 32, 64, 128, 256, 512, 1024, 2048, 4096};
        std::vector<double> Nv, vsDft, vsScalar;
        for (const std::size_t N : sizes)
        {
            const auto x = makeRandom(N, 0xA11CEu + N);
            const auto Xs = soaForward(x);
            vsDft.push_back(maxAbsDiff(Xs, dftRef(x)));
            vsScalar.push_back(maxAbsDiff(Xs, scalarForward(x)));
            Nv.push_back(static_cast<double>(N));
        }
        ok &= writeCsv(dir / "fft_simd_accuracy.csv",
                       {"N", "simd_vs_dft", "simd_vs_scalar"}, {Nv, vsDft, vsScalar});
    }

    // (b) Phase coherence: sweep a fixed-phase tone across bins; record the
    //     phase error (deg) and magnitude error (dB). Flat ~0 => coherent.
    {
        const std::size_t N = 2048;
        const double phase = 0.3;
        const double ref = 0.5 * static_cast<double>(N);
        std::vector<double> bin, phaseErrDeg, magErrDb;
        for (std::size_t k = 1; k < N / 2; k += 4)
        {
            const auto x = makeRealTone(N, static_cast<double>(k), 1.0, phase);
            const auto X = soaForward(x);
            bin.push_back(static_cast<double>(k));
            phaseErrDeg.push_back(wrapPi(std::arg(X[k]) - phase) * 180.0 / kPi);
            magErrDb.push_back(20.0 * std::log10(std::max(std::abs(X[k]) / ref, 1e-12)));
        }
        ok &= writeCsv(dir / "fft_phase_coherence.csv",
                       {"bin", "phase_err_deg", "mag_err_db"}, {bin, phaseErrDeg, magErrDb});
    }

    // (c) Throughput vs N: SoA SIMD vs scalar recursive engine (Msamples/s).
    {
        const std::vector<std::size_t> sizes = {
            64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536
        };
        std::vector<double> Nv, simdMsps, scalarMsps, speedup;
        for (const std::size_t N : sizes)
        {
            const auto x = makeRandom(N, 0xBEEFu + N);
            SoaFFT e;
            e.setSize(N);
            std::vector<Complex> out(N);
            ScalarFFT f(N);

            const double tSimd = medianTimeNs([&]
            {
                e.fwd(x.data(), out.data());
                return out[0].real();
            });
            const double tScalar = medianTimeNs([&]
            {
                const auto r = f.fwd(x);
                return r[0].real();
            });

            Nv.push_back(static_cast<double>(N));
            simdMsps.push_back((static_cast<double>(N) / (tSimd * 1e-9)) / 1e6);
            scalarMsps.push_back((static_cast<double>(N) / (tScalar * 1e-9)) / 1e6);
            speedup.push_back(tScalar / tSimd);
        }
        ok &= writeCsv(dir / "fft_throughput.csv",
                       {"N", "simd_msps", "scalar_msps", "speedup"},
                       {Nv, simdMsps, scalarMsps, speedup});
    }

    return ok;
}
