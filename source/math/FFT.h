#pragma once

#ifndef QUOTIENT_FFT_H
#define QUOTIENT_FFT_H

#include <bit>
#include <complex>
#include <vector>
#include <cstddef>
#include <cassert>
#include <numbers>

namespace MarsDSP::MathOps
{
    // recursive radix-2 cooley-tukey FFT primitive
    class FFT
    {
    public:
        using img = std::complex<double>;
        explicit FFT(const std::size_t size) : sz(size)
        {
            assert(std::has_single_bit(sz) && "FFT size must be a power of two!");
        }

        [[nodiscard]] std::size_t size() const
        {
            return sz;
        }

        [[nodiscard]] std::vector<img> fwd(const std::vector<img> &x) const
        {
            assert(x.size() == sz && "input size must match FFT size!");
            return transform(x);
        }

        [[nodiscard]] std::vector<img> inv(const std::vector<img> &X) const
        {
            assert(X.size() == sz && "input size must match FFT size!");
            const std::size_t N = sz;
            std::vector<img> xc(N);
            for (auto k {0uz}; k < N; ++k)
                xc[k] = std::conj(X[k]);
            std::vector<img> y = transform(std::move(xc));
            for (auto k {0uz}; k < N; ++k)
                y[k] = std::conj(y[k]) / static_cast<double>(N);
            return y;
        }

    private:
        std::size_t sz;
        static constexpr double twoPI = 2 * std::numbers::pi_v<double>;

        // split even/odd, recurse on each half, recombine with butterflies
        static std::vector<img> transform(std::vector<img> x)
        {
            const std::size_t N = x.size();
            if (N == 1) return x;                       // 1-point DFT is identity
            std::vector<img> even(N / 2);
            std::vector<img> odd(N / 2);
            for (auto k {0uz}; k < N / 2; ++k)
            {
                even[k] = x[2 * k];
                odd[k] = x[2 * k + 1];
            }
            const std::vector<img> E = transform(even); // DFT of even samples
            const std::vector<img> O = transform(odd);  // DFT of odd samples
            std::vector<img> result(N);
            for (auto k {0uz}; k < N / 2; ++k)
            {
                double angle = -twoPI * static_cast<double>(k) / static_cast<double>(N);
                img W = std::polar(1.0, angle); // twiddle W_N^k
                img t = W * O[k];               // shared product feeds both outputs
                result[k] = E[k] + t;           // a + W*b
                result[k + N / 2] = E[k] - t;   // a - W*b (W_N^{k+N/2} = -W_N^k)
            }
            return result;
        }
    };
}
#endif