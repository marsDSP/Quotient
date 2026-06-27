#pragma once

#include <bit>
#include <complex>
#include <vector>
#include <cstddef>
#include <cassert>
#include <numbers>

#if !defined(QUOTIENT_FFT_DISABLE_XSIMD) && defined(__has_include) && __has_include(<xsimd/xsimd.hpp>)
#  include <xsimd/xsimd.hpp>
#  define QUOTIENT_FFT_XSIMD_INCLUDED_
#endif

namespace MarsDSP::MathOps {
#ifdef QUOTIENT_FFT_XSIMD_INCLUDED_
    inline constexpr bool haveXSIMD = true;
#else
    inline constexpr bool haveXSIMD = false;
#endif

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

        static constexpr double twoPI = 2 * std::numbers::pi_v<double>;

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
                    outStrideHigh /= radix;
                    outStride = outStrideHigh;
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

        void appendStages(std::size_t factorIndex, const std::size_t offset, const std::size_t length, const std::size_t repeatCount)
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
                        const double phase = twoPI * static_cast<double>(i) * static_cast<double>(r) / static_cast<double>(length);
                        twiddleTable.push_back(Complex(static_cast<ST>(std::cos(phase)), static_cast<ST>(-std::sin(phase))));
                    }
            appendStages(factorIndex + 1, offset, subLength, repeatCount * radix);
            stages.push_back(stage);
        }

        template <bool inverse>
        void execute(const Complex *in, Complex *out) const
        {
            for (const ReorderEntry &entry : reorderTable) out[entry.outputIndex] = in[entry.inputIndex];
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
            if (N == 1) return x;
            if (N == 2)
            {
                const Complex s0 = x[0];
                const Complex s1 = x[1];
                return {s0 + s1, s0 - s1};
            }

            const std::size_t quarter = N / 4;
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
            const std::vector<Complex> dft0 = transform<inverse>(sub0);
            const std::vector<Complex> dft1 = transform<inverse>(sub1);
            const std::vector<Complex> dft2 = transform<inverse>(sub2);
            const std::vector<Complex> dft3 = transform<inverse>(sub3);

            std::vector<Complex> result(N);
            for (auto k {0uz}; k < quarter; ++k)
            {
                const double angle = -twoPI * static_cast<double>(k) / static_cast<double>(N);
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

                result[k] = sum02 + sum13;
                result[k + quarter] = detail::complexAddI<!inverse>(diff02, diff13);
                result[k + 2 * quarter] = sum02 - sum13;
                result[k + 3 * quarter] = detail::complexAddI<inverse>(diff02, diff13);
            }
            return result;
        }
    };

    template <typename ST = double>
    class RealFFT {
    public:
        using Complex = std::complex<ST>;
        RealFFT() = default;
        explicit RealFFT(const std::size_t n)
        {
            setSize(n);
        }

        [[nodiscard]] std::size_t size() const
        {
            return complexEngine.size() * 2;
        }

        std::size_t setSize(const std::size_t n)
        {
            assert(n % 2 == 0 && "RealFFT length must be even!");
            const std::size_t half = n / 2;
            packed.resize(half);
            spectrum.resize(half);

            const std::size_t binCount = half / 2 + 1;
            twiddles.resize(binCount);
            for (auto k {0uz}; k < binCount; ++k)
            {
                const double theta = -twoPI * static_cast<double>(k) / static_cast<double>(n);
                twiddles[k] = Complex(static_cast<ST>(std::cos(theta)), static_cast<ST>(std::sin(theta)));
            }
            complexEngine.setSize(half);
            return n;
        }

        void forward(const ST *in, Complex *out)
        {
            const std::size_t half = complexEngine.size();

            for (auto k {0uz}; k < half; ++k)
                packed[k] = Complex(in[2 * k], in[2 * k + 1]);

            complexEngine.fwd(packed.data(), spectrum.data());

            const ST re = spectrum[0].real();
            const ST im = spectrum[0].imag();
            out[0] = Complex(re + im, re - im);

            for (auto k {1uz}; k <= half / 2; ++k)
            {
                const std::size_t mirror = half - k;
                const Complex lo = spectrum[k];
                const Complex hi = std::conj(spectrum[mirror]);
                const Complex evenBin = (lo + hi) * static_cast<ST>(0.5);
                const Complex halfDiff = (lo - hi) * static_cast<ST>(0.5);
                const Complex oddBin = Complex(halfDiff.imag(), -halfDiff.real());
                const Complex rotated = detail::complexMul<false>(twiddles[k], oddBin);
                out[k] = evenBin + rotated;
                out[mirror] = std::conj(evenBin - rotated);
            }
        }

        void inverse(const Complex *in, ST *out)
        {
            const std::size_t half = complexEngine.size();
            spectrum[0] = Complex(in[0].real() + in[0].imag(), in[0].real() - in[0].imag());
            for (auto k {1uz}; k <= half / 2; ++k)
            {
                const std::size_t mirror = half - k;
                const Complex sum = in[k] + std::conj(in[mirror]);
                const Complex diff = in[k] - std::conj(in[mirror]);
                const Complex twoOdd = detail::complexMul<true>(diff, twiddles[k]);
                const Complex twoIOdd = Complex(-twoOdd.imag(), twoOdd.real());
                spectrum[k] = sum + twoIOdd;
                spectrum[mirror] = std::conj(sum - twoIOdd);
            }
            complexEngine.inverse(spectrum.data(), packed.data());
            for (auto k {0uz}; k < half; ++k)
            {
                out[2*k] = packed[k].real();
                out[2*k+1] = packed[k].imag();
            }
        }

    private:
        static constexpr double twoPI = 2.0 * std::numbers::pi_v<double>;
        FFT<ST> complexEngine;
        std::vector<Complex> packed;
        std::vector<Complex> spectrum;
        std::vector<Complex> twiddles;
    };

    // CTAD
    RealFFT(std::size_t) -> RealFFT<>;
}
