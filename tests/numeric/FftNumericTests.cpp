// ============================================================================
//  FftNumericTests.cpp
//  Correctness tests + chart-data exporter for MarsDSP::MathOps::FFT
// ============================================================================
#include <catch2/catch_test_macros.hpp>

#include "FftChartData.h"
#include "math/FFT.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <numbers>
#include <random>
#include <vector>

namespace
{
    using Complex = std::complex<double>;
    using FFTd    = MarsDSP::MathOps::FFT<double>;

    constexpr double twoPI = 2.0 * std::numbers::pi_v<double>;

    // -------------------------------------------------------------------------
    //  Reference O(N²) DFT — used only for small N in accuracy checks
    // -------------------------------------------------------------------------
    std::vector<Complex> dftRef(const std::vector<Complex> &x)
    {
        const std::size_t N = x.size();
        std::vector<Complex> X(N);
        for (std::size_t k = 0; k < N; ++k)
            for (std::size_t n = 0; n < N; ++n)
            {
                const double phase =
                    -twoPI * static_cast<double>(k * n) / static_cast<double>(N);
                X[k] += x[n] * std::polar(1.0, phase);
            }
        return X;
    }

    double maxAbsDiff(const std::vector<Complex> &a, const std::vector<Complex> &b)
    {
        double m = 0.0;
        for (std::size_t i = 0; i < a.size(); ++i)
            m = std::max(m, std::abs(a[i] - b[i]));
        return m;
    }

    std::vector<Complex> makeRandom(std::size_t N, std::uint64_t seed = 0xDEADBEEFULL)
    {
        std::mt19937_64 rng(seed);
        std::uniform_real_distribution<double> dist(-1.0, 1.0);
        std::vector<Complex> v(N);
        for (auto &z : v) z = Complex{dist(rng), dist(rng)};
        return v;
    }

    // Pure complex exponential e^{j2π k0 n/N}: all energy concentrates at bin k0
    std::vector<Complex> makeComplexTone(std::size_t N, std::size_t k0)
    {
        std::vector<Complex> x(N);
        for (std::size_t n = 0; n < N; ++n)
        {
            const double phase =
                twoPI * static_cast<double>(k0 * n) / static_cast<double>(N);
            x[n] = Complex{std::cos(phase), std::sin(phase)};
        }
        return x;
    }

    // Real cosine: cos(2π k0 n/N)
    std::vector<Complex> makeCosineTone(std::size_t N, std::size_t k0)
    {
        std::vector<Complex> x(N);
        for (std::size_t n = 0; n < N; ++n)
        {
            const double phase =
                twoPI * static_cast<double>(k0 * n) / static_cast<double>(N);
            x[n] = Complex{std::cos(phase), 0.0};
        }
        return x;
    }

    bool writeCsv(const std::filesystem::path        &path,
                  const std::vector<std::string>      &header,
                  const std::vector<std::vector<double>> &cols)
    {
        std::ofstream os(path);
        if (!os) return false;
        for (std::size_t c = 0; c < header.size(); ++c)
            os << header[c] << (c + 1 < header.size() ? "," : "\n");
        const std::size_t rows = cols.empty() ? 0uz : cols.front().size();
        os << std::setprecision(10);
        for (std::size_t r = 0; r < rows; ++r)
            for (std::size_t c = 0; c < cols.size(); ++c)
                os << cols[c][r] << (c + 1 < cols.size() ? "," : "\n");
        return static_cast<bool>(os);
    }
} // namespace


// =============================================================================
//  Test cases
// =============================================================================

TEST_CASE("FFT forward matches reference O(N²) DFT", "[fft][correctness]")
{
    for (const std::size_t N : {2uz, 4uz, 8uz, 16uz, 32uz, 64uz, 128uz, 256uz})
    {
        const auto x    = makeRandom(N, 42 + N);
        FFTd       fft(N);
        const auto X    = fft.fwd(x);
        const auto Xref = dftRef(x);
        const double err = maxAbsDiff(X, Xref);
        CAPTURE(N, err);
        // Rounding error grows at most O(log2 N) * eps * sqrt(N) for radix-2 FFT
        REQUIRE(err < 1e-9 * static_cast<double>(N));
    }
}

TEST_CASE("FFT round-trip identity: inv(fwd(x)) / N == x", "[fft][roundtrip]")
{
    for (const std::size_t N : {4uz, 8uz, 64uz, 256uz, 1024uz, 4096uz})
    {
        const auto x = makeRandom(N, 99 + N);
        FFTd fft(N);
        auto recovered = fft.inv(fft.fwd(x));
        for (auto &v : recovered) v /= static_cast<double>(N);
        const double err = maxAbsDiff(recovered, x);
        CAPTURE(N, err);
        REQUIRE(err < 1e-9);
    }
}

TEST_CASE("FFT linearity: fwd(a·x + b·y) == a·fwd(x) + b·fwd(y)", "[fft][linearity]")
{
    constexpr std::size_t N = 128;
    const auto x = makeRandom(N, 111);
    const auto y = makeRandom(N, 222);
    constexpr double a = 2.5, b = -1.3;

    std::vector<Complex> axy(N);
    for (std::size_t i = 0; i < N; ++i) axy[i] = a * x[i] + b * y[i];

    FFTd fft(N);
    const auto Faxy = fft.fwd(axy);
    const auto Fx   = fft.fwd(x);
    const auto Fy   = fft.fwd(y);

    std::vector<Complex> aFbG(N);
    for (std::size_t i = 0; i < N; ++i) aFbG[i] = a * Fx[i] + b * Fy[i];

    const double err = maxAbsDiff(Faxy, aFbG);
    CAPTURE(err);
    REQUIRE(err < 1e-10);
}

TEST_CASE("FFT unit impulse yields flat (all-ones) spectrum", "[fft][impulse]")
{
    for (const std::size_t N : {8uz, 64uz, 256uz, 1024uz})
    {
        std::vector<Complex> impulse(N);
        impulse[0] = Complex{1.0, 0.0};
        FFTd fft(N);
        const auto X = fft.fwd(impulse);
        for (std::size_t k = 0; k < N; ++k)
        {
            CAPTURE(N, k);
            REQUIRE(std::abs(X[k] - Complex{1.0, 0.0}) < 1e-10);
        }
    }
}

TEST_CASE("FFT DC signal concentrates all energy at bin 0", "[fft][dc]")
{
    for (const std::size_t N : {8uz, 32uz, 128uz, 512uz})
    {
        const std::vector<Complex> dc(N, Complex{1.0, 0.0});
        FFTd fft(N);
        const auto X = fft.fwd(dc);

        REQUIRE(std::abs(X[0] - Complex{static_cast<double>(N), 0.0}) < 1e-8);
        for (std::size_t k = 1; k < N; ++k)
        {
            CAPTURE(N, k);
            REQUIRE(std::abs(X[k]) < 1e-8);
        }
    }
}

TEST_CASE("FFT complex exponential concentrates energy at one bin", "[fft][tone]")
{
    constexpr std::size_t N = 256;
    FFTd fft(N);
    for (const std::size_t k0 : {1uz, 7uz, 32uz, 64uz})
    {
        const auto X = fft.fwd(makeComplexTone(N, k0));
        // X[k0] should equal N; all other bins should be negligible
        REQUIRE(std::abs(X[k0] - Complex{static_cast<double>(N), 0.0}) < 1e-6);
        for (std::size_t k = 0; k < N; ++k)
        {
            if (k == k0) continue;
            CAPTURE(N, k0, k);
            REQUIRE(std::abs(X[k]) < 1e-6);
        }
    }
}

TEST_CASE("FFT Parseval's theorem: sum|x|² == sum|X|² / N", "[fft][parseval]")
{
    for (const std::size_t N : {8uz, 64uz, 256uz, 1024uz, 4096uz})
    {
        const auto x = makeRandom(N);
        FFTd fft(N);
        const auto X = fft.fwd(x);

        double Et = 0.0;
        for (const auto &s : x) Et += std::norm(s);

        double Ef = 0.0;
        for (const auto &S : X) Ef += std::norm(S);
        Ef /= static_cast<double>(N);

        CAPTURE(N, Et, Ef);
        REQUIRE(std::abs(Et - Ef) < 1e-8 * Et);
    }
}

TEST_CASE("FFT Hermitian symmetry for real-valued inputs", "[fft][symmetry]")
{
    for (const std::size_t N : {16uz, 64uz, 256uz})
    {
        // Imaginary part set to zero → X[N-k] must equal conj(X[k])
        const auto r = makeRandom(N, 777);
        std::vector<Complex> x(N);
        for (std::size_t i = 0; i < N; ++i) x[i] = Complex{r[i].real(), 0.0};

        FFTd fft(N);
        const auto X = fft.fwd(x);
        for (std::size_t k = 1; k < N / 2; ++k)
        {
            CAPTURE(N, k);
            REQUIRE(std::abs(X[N - k] - std::conj(X[k])) < 1e-10);
        }
    }
}


// =============================================================================
//  Chart-data export
// =============================================================================

bool quotient_fft_charts::exportChartData(const std::string &outDir)
{
    namespace fs = std::filesystem;
    std::error_code ec;
    fs::create_directories(outDir, ec);
    if (ec) return false;

    const fs::path dir(outDir);
    bool ok = true;

    // -------------------------------------------------------------------------
    // (a) Magnitude spectrum of a real cosine tone
    //     N=512, tone at bin 40; normalised so peak = 0 dB.
    // -------------------------------------------------------------------------
    {
        constexpr std::size_t N  = 512;
        constexpr std::size_t k0 = 40;
        FFTd fft(N);
        const auto X = fft.fwd(makeCosineTone(N, k0));

        // Cosine gives |X[k0]| = |X[N-k0]| = N/2.  Divide by N/2 → 0 dB peak.
        const double norm = 0.5 * static_cast<double>(N);
        std::vector<double> bins(N), magDb(N);
        for (std::size_t k = 0; k < N; ++k)
        {
            bins[k]  = static_cast<double>(k);
            magDb[k] = 20.0 * std::log10(std::max(std::abs(X[k]) / norm, 1e-12));
        }
        ok &= writeCsv(dir / "fft_spectrum.csv", {"bin", "magnitude_db"}, {bins, magDb});
    }

    // -------------------------------------------------------------------------
    // (b) Forward accuracy and round-trip error vs transform size.
    //     DFT reference is used for all sizes; acceptable since DFT for
    //     N ≤ 4096 completes in well under a second.
    // -------------------------------------------------------------------------
    {
        const std::vector<std::size_t> sizes = {
            16, 32, 64, 128, 256, 512, 1024, 2048, 4096
        };
        std::vector<double> Nvec, fwdErr, rtErr;
        Nvec.reserve(sizes.size());
        fwdErr.reserve(sizes.size());
        rtErr.reserve(sizes.size());

        for (const std::size_t N : sizes)
        {
            const auto x = makeRandom(N, 0xABCDEFULL);
            FFTd fft(N);
            const auto X = fft.fwd(x);

            fwdErr.push_back(maxAbsDiff(X, dftRef(x)));

            auto rec = fft.inv(X);
            for (auto &v : rec) v /= static_cast<double>(N);
            rtErr.push_back(maxAbsDiff(rec, x));

            Nvec.push_back(static_cast<double>(N));
        }
        ok &= writeCsv(dir / "fft_accuracy.csv",
                       {"N", "fwd_max_error", "roundtrip_max_error"},
                       {Nvec, fwdErr, rtErr});
    }

    // -------------------------------------------------------------------------
    // (c) Parseval relative energy error vs transform size
    // -------------------------------------------------------------------------
    {
        const std::vector<std::size_t> sizes = {
            8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096
        };
        std::vector<double> Nvec, parsErr;
        Nvec.reserve(sizes.size());
        parsErr.reserve(sizes.size());

        for (const std::size_t N : sizes)
        {
            const auto x = makeRandom(N, 0xFADEFACEULL);
            FFTd fft(N);
            const auto X = fft.fwd(x);

            double Et = 0.0;
            for (const auto &s : x) Et += std::norm(s);
            double Ef = 0.0;
            for (const auto &S : X) Ef += std::norm(S);
            Ef /= static_cast<double>(N);

            parsErr.push_back(Et > 0.0 ? std::abs(Et - Ef) / Et : 0.0);
            Nvec.push_back(static_cast<double>(N));
        }
        ok &= writeCsv(dir / "fft_parseval.csv",
                       {"N", "parseval_relative_error"},
                       {Nvec, parsErr});
    }

    return ok;
}
