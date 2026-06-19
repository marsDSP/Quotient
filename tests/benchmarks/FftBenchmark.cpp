// A/B micro-benchmark for the FFT butterfly helpers.
//
// "baseline" below is a verbatim copy of the FFT recursion as it existed
// before complexMul was integrated: it uses std::complex's operator* in the
// butterfly and the conjugate-trick inverse. The "optimized" side is the real
// production class, which now uses detail::complexMul.

#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>

#include "math/FFT.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <numbers>
#include <random>
#include <string>
#include <vector>

namespace
{
    using img = std::complex<double>;

    std::vector<img> makeInput(std::size_t N)
    {
        std::mt19937_64 rng(0xC0FFEEULL);
        std::uniform_real_distribution<double> dist(-1.0, 1.0);
        std::vector<img> v(N);
        for (auto &z : v)
            z = img{dist(rng), dist(rng)};
        return v;
    }

    double maxAbsDiff(const std::vector<img> &a, const std::vector<img> &b)
    {
        double m = 0.0;
        for (std::size_t i = 0; i < a.size(); ++i)
            m = std::max(m, std::abs(a[i] - b[i]));
        return m;
    }
}

// ---- baseline: pre-integration implementation (std::complex operator*) -------
namespace baseline
{
    using img = std::complex<double>;

    inline std::vector<img> transform(std::vector<img> x)
    {
        const std::size_t N = x.size();
        if (N == 1) return x;
        std::vector<img> even(N / 2);
        std::vector<img> odd(N / 2);
        for (std::size_t k = 0; k < N / 2; ++k)
        {
            even[k] = x[2 * k];
            odd[k] = x[2 * k + 1];
        }
        const std::vector<img> E = transform(even);
        const std::vector<img> O = transform(odd);
        std::vector<img> result(N);
        constexpr double twoPI = 2 * std::numbers::pi_v<double>;
        for (std::size_t k = 0; k < N / 2; ++k)
        {
            const double angle = -twoPI * static_cast<double>(k) / static_cast<double>(N);
            const img W = std::polar(1.0, angle);
            const img t = W * O[k];                 // <-- std::complex operator*
            result[k] = E[k] + t;
            result[k + N / 2] = E[k] - t;
        }
        return result;
    }

    inline std::vector<img> fwd(const std::vector<img> &x) { return transform(x); }

    inline std::vector<img> inv(const std::vector<img> &X)
    {
        const std::size_t N = X.size();
        std::vector<img> xc(N);
        for (std::size_t k = 0; k < N; ++k)
            xc[k] = std::conj(X[k]);
        std::vector<img> y = transform(std::move(xc));
        for (std::size_t k = 0; k < N; ++k)
            y[k] = std::conj(y[k]) / static_cast<double>(N);
        return y;
    }
}

TEST_CASE("FFT helpers match the baseline and round-trip", "[fft]")
{
    for (const std::size_t N : {std::size_t{2}, std::size_t{8}, std::size_t{64}, std::size_t{1024}})
    {
        const auto in = makeInput(N);
        MarsDSP::MathOps::FFT fft(N);

        const auto optF = fft.fwd(in);
        const auto baseF = baseline::fwd(in);
        REQUIRE(maxAbsDiff(optF, baseF) < 1e-9);          // forward matches old code

        const auto optR = fft.inv(optF);
        REQUIRE(maxAbsDiff(optR, in) < 1e-9);             // round-trip identity
        REQUIRE(maxAbsDiff(optR, baseline::inv(baseF)) < 1e-9); // inverse matches old code
    }
}

TEST_CASE("FFT forward: complexMul vs std::complex operator*", "[fft][bench]")
{
    for (const std::size_t N : {std::size_t{1} << 10, std::size_t{1} << 12, std::size_t{1} << 14})
    {
        const auto in = makeInput(N);
        MarsDSP::MathOps::FFT fft(N);
        BENCHMARK("baseline  fwd  N=" + std::to_string(N)) { return baseline::fwd(in); };
        BENCHMARK("complexMul fwd N=" + std::to_string(N)) { return fft.fwd(in); };
    }
}

TEST_CASE("FFT inverse: complexMul<true> vs conjugate-trick", "[fft][bench]")
{
    for (const std::size_t N : {std::size_t{1} << 10, std::size_t{1} << 12, std::size_t{1} << 14})
    {
        const auto in = makeInput(N);
        MarsDSP::MathOps::FFT fft(N);
        const auto spec = fft.fwd(in);
        BENCHMARK("baseline  inv  N=" + std::to_string(N)) { return baseline::inv(spec); };
        BENCHMARK("complexMul inv N=" + std::to_string(N)) { return fft.inv(spec); };
    }
}
