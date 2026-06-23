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

    template <typename ST = double>
    class FFT {
    public:
        using Complex = std::complex<double>;
        explicit FFT(const std::size_t size) : fftSize(size)
        {
            assert(std::has_single_bit(fftSize) && "FFT size must be a power of two!");
        }

        [[nodiscard]] std::size_t size() const
        {
            return fftSize;
        }

        [[nodiscard]] std::vector<Complex> fwd(const std::vector<Complex> &x) const
        {
            assert(x.size() == fftSize && "input size must match FFT size!");
            return transform<false>(x);
        }

        [[nodiscard]] std::vector<Complex> inv(const std::vector<Complex> &X) const
        {
            assert(X.size() == fftSize && "input size must match FFT size!");
            return transform<true>(X);
        }

    private:
        enum class StageKind
        {
            mixed, radix2, radix3, radix4
        };

        struct Stage
        {
            StageKind kind;
            std::size_t radix;
            std::size_t offset;
            std::size_t innerCount;
            std::size_t outerCount;
            std::size_t twiddleOffset;
        };

        struct ReorderEntry
        {
            std::size_t outputIndex;
            std::size_t inputIndex;
        };

        std::size_t fftSize;
        std::size_t planLength = 0;
        std::vector<std::size_t> radixFactors;
        std::vector<Complex> twiddleTable;
        std::vector<ReorderEntry> reorderTable;
        std::vector<Stage> stages;

        static constexpr double twoPi = 2 * std::numbers::pi_v<double>;

        void buildPlan()
        {
            radixFactors.clear();
            {
                std::size_t remaining = planLength;
                std::size_t candidate = 2;
                while (remaining > 1)
                {
                    if (remaining % candidate == 0) { radixFactors.push_back(candidate); remaining /= candidate; }
                    else if (candidate * candidate > remaining) { candidate = remaining; }
                    else { ++candidate; }
                }
            }

            stages.clear();
            twiddleTable.clear();
            if (planLength > 1) appendStages(0, 0, planLength, 1);
            reorderTable.clear();
            reorderTable.push_back(ReorderEntry{0, 0});
            std::size_t lowFactorIdx = 0;
            std::size_t highFactorIdx = radixFactors.size();
            std::size_t inStrideLow = planLength;
            std::size_t outStrideLow = 1;
            std::size_t inStrideHigh = 1;
            std::size_t outStrideHigh = planLength;

            while (outStrideLow * inStrideHigh < planLength)
            {
                std::size_t radix;
                std::size_t inStride;
                std::size_t outStride;
                if (outStrideLow <= inStrideHigh)
                {
                    radix = radixFactors[lowFactorIdx++];
                    inStride = (inStrideLow /= radix);
                    outStride = outStrideLow;
                    outStrideLow *= radix;
                }
                else
                {
                    radix = radixFactors[--highFactorIdx];
                    inStride = inStrideHigh;
                    inStrideHigh *= radix;
                    outStride = (outStrideHigh /= radix);
                }
                const std::size_t prevCount = reorderTable.size();
                for (auto i {1uz}; i < radix; ++i)
                    for (auto j {0uz}; j < prevCount; ++j)
                    {
                        ReorderEntry entry = reorderTable[j];
                        entry.outputIndex += i * inStride;
                        entry.inputIndex += i * outStride;
                        reorderTable.push_back(entry);
                    }
            }
        }

        void appendStages(std::size_t factorIndex,
                          const std::size_t offset,
                          const std::size_t length,
                          const std::size_t repeatCount)
        {
            if (factorIndex >= radixFactors.size()) return;
            std::size_t radix = radixFactors[factorIndex];
            if (factorIndex + 1 < radixFactors.size()
                && radixFactors[factorIndex] == 2
                && radixFactors[factorIndex + 1] == 2)
            {
                ++factorIndex;
                radix = 4;
            }

            const std::size_t subLength = length / radix;
            Stage stage { StageKind::mixed, radix, offset, subLength, repeatCount, twiddleTable.size() };
            if (radix == 2) stage.kind = StageKind::radix2;
            else if (radix == 3) stage.kind = StageKind::radix3;
            else if (radix == 4) stage.kind = StageKind::radix4;

            bool reusedTwiddles = false;
            for (const Stage &existing : stages)
                if (existing.radix == stage.radix && existing.innerCount == stage.innerCount)
                {
                    stage.twiddleOffset = existing.twiddleOffset;
                    reusedTwiddles = true;
                    break;
                }
            if (!reusedTwiddles)
                for (auto i {0uz}; i < subLength; ++i)
                    for (auto r {0uz}; r < radix; ++r)
                    {
                        const double phase = twoPi * static_cast<double>(i) *
                                                     static_cast<double>(r) /
                                                     static_cast<double>(length);
                        twiddleTable.push_back(Complex(static_cast<ST>(std::cos(phase)),
                                                       static_cast<ST>(-std::sin(phase))));
                    }
            appendStages(factorIndex + 1, offset, subLength, repeatCount * radix);
            stages.push_back(stage);
        }

        template <bool inverse>
        void execute(const Complex *in, Complex *out) const
        {
            for (const ReorderEntry &entry : reorderTable)
                out[entry.outputIndex] = in[entry.inputIndex];
            for (const Stage &stage : stages)
                switch (stage.kind)
                {
                    case StageKind::radix2: butterflyRadix2<inverse>(out + stage.offset, stage); break;
                    case StageKind::radix3: butterflyRadix3<inverse>(out + stage.offset, stage); break;
                    case StageKind::radix4: butterflyRadix4<inverse>(out + stage.offset, stage); break;
                    case StageKind::mixed: butterflyMixed<inverse>(out + stage.offset, stage); break;
                    default:;
                }
        }

        template <bool inverse>
        static std::vector<Complex> transform(std::vector<Complex> x)
        {
            const std::size_t N = x.size();
            if (N == 1) return x;                       // 1-point DFT is identity
            if (N == 2)                                 // radix-2 butterfly base case
            {
                const Complex s0 = x[0];
                const Complex s1 = x[1];
                return {s0 + s1, s0 - s1};              // W_2^0 = 1, W_2^1 = -1
            }

            const std::size_t quarter = N / 4;          // sub-transform length
            std::vector<Complex> sub0(quarter);
            std::vector<Complex> sub1(quarter);
            std::vector<Complex> sub2(quarter);
            std::vector<Complex> sub3(quarter);
            for (auto k{0uz}; k < quarter; ++k)
            {
                sub0[k] = x[4 * k];
                sub1[k] = x[4 * k + 1];
                sub2[k] = x[4 * k + 2];
                sub3[k] = x[4 * k + 3];
            }
            const std::vector<Complex> dft0 = transform<inverse>(sub0); // DFT of n%4==0
            const std::vector<Complex> dft1 = transform<inverse>(sub1); // DFT of n%4==1
            const std::vector<Complex> dft2 = transform<inverse>(sub2); // DFT of n%4==2
            const std::vector<Complex> dft3 = transform<inverse>(sub3); // DFT of n%4==3

            std::vector<Complex> result(N);
            for (auto k {0uz}; k < quarter; ++k)
            {
                const double angle = -twoPi * static_cast<double>(k) / static_cast<double>(N);
                const Complex tw1 = std::polar(1.0, angle);
                const Complex tw2 = detail::complexMul<false>(tw1, tw1);
                const Complex tw3 = detail::complexMul<false>(tw2, tw1);

                const Complex term0 = dft0[k];
                const Complex term1 = inverse ? detail::complexMul<true>(dft1[k], tw1)
                                              : detail::complexMul<false>(tw1, dft1[k]);
                const Complex term2 = inverse ? detail::complexMul<true>(dft2[k], tw2)
                                              : detail::complexMul<false>(tw2, dft2[k]);
                const Complex term3 = inverse ? detail::complexMul<true>(dft3[k], tw3)
                                              : detail::complexMul<false>(tw3, dft3[k]);

                const Complex sum02 = term0 + term2;
                const Complex diff02 = term0 - term2;
                const Complex sum13 = term1 + term3;
                const Complex diff13 = term1 - term3;

                result[k]               = sum02 + sum13;
                result[k + quarter]     = detail::complexAddI<!inverse>(diff02, diff13);
                result[k + 2 * quarter] = sum02 - sum13;
                result[k + 3 * quarter] = detail::complexAddI<inverse>(diff02, diff13);
            }
            return result;
        }
    };
}
#endif