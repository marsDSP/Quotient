#pragma once

#ifndef QUOTIENT_FFT_H
#define QUOTIENT_FFT_H

#include <bit>
#include <complex>
#include <vector>
#include <cstddef>
#include <cassert>
#include <numbers>

namespace MarsDSP::MathOps {
    namespace detail {
        template <bool conjugateSecond, typename V>
        inline std::complex<V> complexMul(const std::complex<V> &a, const std::complex<V> &b)
        {
            if constexpr (conjugateSecond)
                return std::complex<V> {
                    b.real() * a.real() + b.imag() * a.imag(),
                    b.real() * a.imag() - b.imag() * a.real()
                };
            else
                return std::complex<V> {
                    a.real() * b.real() - a.imag() * b.imag(),
                    a.real() * b.imag() + a.imag() * b.real()
                };
        }

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
            return transform<true>(X);
        }

    private:
        std::size_t sz;
        static constexpr double twoPI = 2 * std::numbers::pi_v<double>;

        template <bool inverse>
        static std::vector<img> transform(std::vector<img> x)
        {
            const std::size_t N = x.size();
            if (N == 1) return x;                       // 1-point DFT is identity
            if (N == 2)                                 // radix-2 butterfly base case
            {
                const img a = x[0];
                const img b = x[1];
                return {a + b, a - b};                  // W_2^0 = 1, W_2^1 = -1
            }

            const std::size_t M = N / 4;                // sub-transform length
            std::vector<img> x0(M), x1(M), x2(M), x3(M);
            for (auto k {0uz}; k < M; ++k)
            {
                x0[k] = x[4 * k];
                x1[k] = x[4 * k + 1];
                x2[k] = x[4 * k + 2];
                x3[k] = x[4 * k + 3];
            }
            const std::vector<img> X0 = transform<inverse>(x0); // DFT of n%4==0
            const std::vector<img> X1 = transform<inverse>(x1); // DFT of n%4==1
            const std::vector<img> X2 = transform<inverse>(x2); // DFT of n%4==2
            const std::vector<img> X3 = transform<inverse>(x3); // DFT of n%4==3

            std::vector<img> result(N);
            for (auto k {0uz}; k < M; ++k)
            {
                const double angle = -twoPI * static_cast<double>(k) / static_cast<double>(N);
                const img W1 = std::polar(1.0, angle);
                const img W2 = detail::complexMul<false>(W1, W1);
                const img W3 = detail::complexMul<false>(W2, W1);

                const img t0 = X0[k];
                const img t1 = inverse ? detail::complexMul<true>(X1[k], W1)
                                       : detail::complexMul<false>(W1, X1[k]);
                const img t2 = inverse ? detail::complexMul<true>(X2[k], W2)
                                       : detail::complexMul<false>(W2, X2[k]);
                const img t3 = inverse ? detail::complexMul<true>(X3[k], W3)
                                       : detail::complexMul<false>(W3, X3[k]);

                const img sum02 = t0 + t2;
                const img dif02 = t0 - t2;
                const img sum13 = t1 + t3;
                const img dif13 = t1 - t3;

                result[k]         = sum02 + sum13;
                result[k + M]     = detail::complexAddI<!inverse>(dif02, dif13);
                result[k + 2 * M] = sum02 - sum13;
                result[k + 3 * M] = detail::complexAddI<inverse>(dif02, dif13);
            }
            return result;
        }
    };
}
#endif