#pragma once

#ifndef QUOTIENT_FFT_H
#define QUOTIENT_FFT_H

#include <bit>
#include <cmath>
#include <array>
#include <memory>
#include <complex>
#include <vector>
#include <cstddef>
#include <cassert>
#include <numbers>

namespace MarsDSP::MathOps {
    namespace detail {
        // Branch-free complex multiply that sidesteps std::complex's slow,
        // Annex-G-compliant operator* (no inf/NaN fix-up; fully inlinable and
        // vectorizable). conjugateSecond == true computes a * conj(b), which the
        // inverse transform uses to flip the twiddle sign for free.
        template <bool conjugateSecond, typename V>
        inline std::complex<V> complexMul(const std::complex<V> &a, const std::complex<V> &b)
        {
            if constexpr (conjugateSecond)
                return std::complex<V>{
                    b.real() * a.real() + b.imag() * a.imag(),
                    b.real() * a.imag() - b.imag() * a.real()
                };
            else
                return std::complex<V>{
                    a.real() * b.real() - a.imag() * b.imag(),
                    a.real() * b.imag() + a.imag() * b.real()
                };
        }

        // a +/- i*b with no multiplies (radix-4 / split-radix butterfly
        // primitive). flipped == false -> a + i*b; flipped == true -> a - i*b.
        // Unused by the radix-2 butterfly below; provided for a future radix-4 path.
        template <bool flipped, typename V>
        inline std::complex<V> complexAddI(const std::complex<V> &a, const std::complex<V> &b)
        {
            if constexpr (flipped)
                return std::complex<V>{a.real() + b.imag(), a.imag() - b.real()};
            else
                return std::complex<V>{a.real() - b.imag(), a.imag() + b.real()};
        }
    }

    class FFT {
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
            return transform<false>(x);
        }

        [[nodiscard]] std::vector<img> inv(const std::vector<img> &X) const
        {
            assert(X.size() == sz && "input size must match FFT size!");
            // Inverse reuses the same recursion with conjugated twiddles
            // (transform<true>), so the two extra std::conj passes the previous
            // conjugate-trick required are gone; just normalise by N.
            const auto N = static_cast<double>(sz);
            std::vector<img> y = transform<true>(X);
            for (auto &v : y)
                v /= N;
            return y;
        }

    private:
        std::size_t sz;
        static constexpr double twoPI = 2 * std::numbers::pi_v<double>;

        // split even/odd, recurse on each half, recombine with butterflies.
        // inverse == true flips the twiddle sign via complexMul's conjugate path.
        template <bool inverse>
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
            const std::vector<img> E = transform<inverse>(even); // DFT of even samples
            const std::vector<img> O = transform<inverse>(odd);  // DFT of odd samples
            std::vector<img> result(N);
            for (auto k {0uz}; k < N / 2; ++k)
            {
                const double angle = -twoPI * static_cast<double>(k) / static_cast<double>(N);
                const img W = std::polar(1.0, angle);   // forward twiddle W_N^k = e^{-i2(pi)k/N}
                // forward: W * O[k]; inverse: conj(W) * O[k] == O[k] * conj(W)
                const img t = inverse ? detail::complexMul<true>(O[k], W)
                                      : detail::complexMul<false>(W, O[k]);
                result[k] = E[k] + t;                   // a + W*b
                result[k + N / 2] = E[k] - t;           // a - W*b (W_N^{k+N/2} = -W_N^k)
            }
            return result;
        }
    };
}
#endif