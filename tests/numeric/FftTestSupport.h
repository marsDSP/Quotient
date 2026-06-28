#pragma once
// ============================================================================
//  FftTestSupport.h
//  Inline helpers shared between the FFT numeric and SIMD test translation
//  units: a reference DFT, signal generators, wrappers around the SoA/SIMD
//  engine and the scalar recursive FFT, a median-timing helper, and a CSV
//  writer.
// ============================================================================
#include "math/FFT.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <numbers>
#include <random>
#include <string>
#include <vector>

namespace ffttest
{
    using Complex = std::complex<double>;
    using ScalarFFT = MarsDSP::MathOps::FFT<double>;
    using SoaFFT = MarsDSP::MathOps::detail::SOAPow2FFT<double>;

    inline constexpr double kTwoPi = 2.0 * std::numbers::pi_v<double>;

    // ------------------------------------------------------------------------
    //  Reference O(N^2) DFT — ground truth for correctness checks.
    // ------------------------------------------------------------------------
    inline std::vector<Complex> dftRef(const std::vector<Complex> &x)
    {
        const std::size_t N = x.size();
        std::vector<Complex> X(N);
        for (std::size_t k = 0; k < N; ++k)
        {
            Complex acc{0.0, 0.0};
            for (std::size_t n = 0; n < N; ++n)
            {
                const double phase =
                    -kTwoPi * static_cast<double>(k * n) / static_cast<double>(N);
                acc += x[n] * std::polar(1.0, phase);
            }
            X[k] = acc;
        }
        return X;
    }

    inline double maxAbsDiff(const std::vector<Complex> &a, const std::vector<Complex> &b)
    {
        double m = 0.0;
        const std::size_t n = std::min(a.size(), b.size());
        for (std::size_t i = 0; i < n; ++i)
            m = std::max(m, std::abs(a[i] - b[i]));
        return m;
    }

    // ------------------------------------------------------------------------
    //  Signal generators
    // ------------------------------------------------------------------------
    inline std::vector<Complex> makeRandom(std::size_t N, std::uint64_t seed = 0xDEADBEEFULL)
    {
        std::mt19937_64 rng(seed);
        std::uniform_real_distribution<double> dist(-1.0, 1.0);
        std::vector<Complex> v(N);
        for (auto &z : v) z = Complex{dist(rng), dist(rng)};
        return v;
    }

    // Pure complex exponential e^{j 2pi k0 n / N}; all energy at bin k0.
    inline std::vector<Complex> makeComplexTone(std::size_t N, std::size_t k0)
    {
        std::vector<Complex> x(N);
        for (std::size_t n = 0; n < N; ++n)
        {
            const double ph = kTwoPi * static_cast<double>(k0 * n) / static_cast<double>(N);
            x[n] = Complex{std::cos(ph), std::sin(ph)};
        }
        return x;
    }

    // Real cosine A*cos(2pi k0 n / N + phase). For 0 < k0 < N/2 the forward FFT
    // gives X[k0] = (A*N/2) e^{j phase}.
    inline std::vector<Complex> makeRealTone(std::size_t N, double k0, double amp = 1.0, double phase = 0.0)
    {
        std::vector<Complex> x(N);
        for (std::size_t n = 0; n < N; ++n)
        {
            const double ph = kTwoPi * k0 * static_cast<double>(n) / static_cast<double>(N) + phase;
            x[n] = Complex{amp * std::cos(ph), 0.0};
        }
        return x;
    }

    inline std::vector<Complex> makeCosineTone(std::size_t N, std::size_t k0)
    {
        return makeRealTone(N, static_cast<double>(k0), 1.0, 0.0);
    }

    // ------------------------------------------------------------------------
    //  Engine wrappers (std::vector in/out)
    // ------------------------------------------------------------------------
    inline std::vector<Complex> soaForward(const std::vector<Complex> &x)
    {
        SoaFFT e;
        e.setSize(x.size());
        std::vector<Complex> out(x.size());
        e.fwd(x.data(), out.data());
        return out;
    }

    inline std::vector<Complex> soaInverse(const std::vector<Complex> &X)
    {
        SoaFFT e;
        e.setSize(X.size());
        std::vector<Complex> out(X.size());
        e.inv(X.data(), out.data());
        return out;
    }

    inline std::vector<Complex> scalarForward(const std::vector<Complex> &x)
    {
        ScalarFFT f(x.size());
        return f.fwd(x);
    }

    inline std::vector<Complex> scalarInverse(const std::vector<Complex> &X)
    {
        ScalarFFT f(X.size());
        return f.inv(X);
    }

    // ------------------------------------------------------------------------
    //  Median wall-clock time (nanoseconds) of fn() over `iters` runs after
    //  `warmup` untimed runs. fn must return a double derived from its output
    //  so the work cannot be optimised away.
    // ------------------------------------------------------------------------
    template <class Fn>
    inline double medianTimeNs(Fn &&fn, int iters = 9, int warmup = 2)
    {
        using clock = std::chrono::steady_clock;
        volatile double sink = 0.0;
        for (int i = 0; i < warmup; ++i)
            sink += fn();
        std::vector<double> times;
        times.reserve(static_cast<std::size_t>(iters));
        for (int i = 0; i < iters; ++i)
        {
            const auto t0 = clock::now();
            const double v = fn();
            const auto t1 = clock::now();
            sink += v;
            times.push_back(std::chrono::duration<double, std::nano>(t1 - t0).count());
        }
        (void) sink;
        std::sort(times.begin(), times.end());
        return times[times.size() / 2];
    }

    // ------------------------------------------------------------------------
    //  CSV writer (column-major input).
    // ------------------------------------------------------------------------
    inline bool writeCsv(const std::filesystem::path &path,
                         const std::vector<std::string> &header,
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
} // namespace ffttest
