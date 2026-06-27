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
            const V ar = a.real();
            const V ai = a.imag();
            const V br = b.real();
            const V bi = conjugateSecond ? -b.imag() : b.imag();
            return std::complex<V>{ar * br - ai * bi, ar * bi + ai * br};
        }

        template <bool flipped, typename V>
        inline std::complex<V> complexAddI(const std::complex<V> &a, const std::complex<V> &b)
        {
            if constexpr (flipped)
                return std::complex<V>{a.real() + b.imag(), a.imag() - b.real()};
            else
                return std::complex<V>{a.real() - b.imag(), a.imag() + b.real()};
        }

        template <typename ST>
        class SOAPow2FFT {
        public:
            using Complex = std::complex<ST>;
            [[nodiscard]] std::size_t size() const { return N; }

            void setSize(const std::size_t n)
            {
                N = n;
                aRe.assign(n, ST(0)); aIm.assign(n, ST(0));
                bRe.assign(n, ST(0)); bIm.assign(n, ST(0));
                schedule.clear();
                twRe.clear(); twIm.clear();
                tw4Re.clear(); tw4Im.clear();

                std::size_t Ns = 1;
                while (Ns*4 <= n)
                {
                    const std::size_t off = tw4Re.size();
                    for (auto r {1uz}; r <= 3; ++r)
                    {
                        for (auto k {0uz}; k < Ns; ++k)
                        {
                            const double th = -twoPI * static_cast<double>(r) *
                                                       static_cast<double>(k) /
                                                       static_cast<double>(4*Ns);

                            tw4Re.push_back(static_cast<ST>(std::cos(th)));
                            tw4Im.push_back(static_cast<ST>(std::sin(th)));
                        }
                    }
                    schedule.push_back({4, Ns, off});
                    Ns *= 4;
                }

                if (Ns < n)
                {
                    const std::size_t off = twRe.size();
                    for (auto k {0uz}; k < Ns; ++k)
                    {
                        const double th = -twoPI * static_cast<double>(k) /
                                                   static_cast<double>(2*Ns);

                        twRe.push_back(static_cast<ST>(std::cos(th)));
                        twIm.push_back(static_cast<ST>(std::sin(th)));
                    }
                    schedule.push_back({2, Ns, off});
                }
            }

            void fwd(const Complex *in, Complex *out) { run<false>(in, out); }
            void inv(const Complex *in, Complex *out) { run<true>(in, out); }

        private:
            std::size_t N = 0;
            static constexpr double twoPI = 2.0 * std::numbers::pi_v<double>;

            struct StageInfo
            {
                std::size_t radix;
                std::size_t Ns;
                std::size_t twOff;
            };

            std::vector<StageInfo> schedule;

#ifdef QUOTIENT_FFT_XSIMD_INCLUDED_
            template <class U> using Alloc = xsimd::aligned_allocator<U>;
#else
            template <class U> using Alloc = std::allocator<U>;
#endif
            using AVec = std::vector<ST, Alloc<ST>>;

            AVec twRe;
            AVec twIm;  // radix-2 cleanup twiddles (single stage, Ns = N/2)
            AVec tw4Re;
            AVec tw4Im; // radix-4 twiddles: per stage [tw1|tw2|tw3] over k
            AVec aRe;
            AVec aIm;
            AVec bRe;
            AVec bIm;

            template <bool inverse>
            void run(const Complex *in, Complex *out)
            {
                for (auto i {0uz}; i < N; ++i)
                {
                    aRe[i] = in[i].real();
                    aIm[i] = in[i].imag();
                }
                ST *sr = aRe.data();
                ST *si = aIm.data();
                ST *dr = bRe.data();
                ST *di = bIm.data();

                for (const StageInfo &st : schedule)
                {
                    if (st.radix == 4) stage4<inverse>(sr, si, dr, di, st.Ns, st.twOff);
                    else stage<inverse>(sr, si, dr, di, st.Ns, st.twOff);
                    std::swap(sr, dr);
                    std::swap(si, di);
                }
                for (auto i {0uz}; i < N; ++i) out[i] = Complex(sr[i], si[i]);
            }

            template <bool inverse>
            void stage(const ST *__restrict sr, const ST *__restrict si, ST *__restrict dr, ST *__restrict di,
                       const std::size_t Ns,    const std::size_t off)
            {
                const std::size_t half = N >> 1;
                const std::size_t nblk = half / Ns;
                const ST *twR = twRe.data();
                const ST *twI = twIm.data();
#ifdef QUOTIENT_FFT_XSIMD_INCLUDED_
                using B = xsimd::batch<ST>;
                constexpr std::size_t W = B::size;
                if (Ns >= W)
                {
                    for (auto blk {0uz}; blk < nblk; ++blk)
                    {
                        const std::size_t jb = blk * Ns;
                        const std::size_t ob = blk * 2 * Ns;
                        for (auto k {0uz}; k < Ns; k += W)
                        {
                            const std::size_t j = jb + k;
                            const B v0r = B::load_aligned(sr + j);
                            const B v0i = B::load_aligned(si + j);
                            const B ar = B::load_aligned(sr + j + half);
                            const B ai = B::load_aligned(si + j + half);
                            const B wr = B::load_unaligned(twR + off + k);
                            B wi = B::load_unaligned(twI + off + k);
                            if constexpr (inverse) wi = -wi;
                            const B v1r = xsimd::fnma(ai, wi, ar * wr);
                            const B v1i = xsimd::fma(ai, wr, ar * wi);
                            (v0r + v1r).store_aligned(dr + ob + k);
                            (v0i + v1i).store_aligned(di + ob + k);
                            (v0r - v1r).store_aligned(dr + ob + k + Ns);
                            (v0i - v1i).store_aligned(di + ob + k + Ns);
                        }
                    }
                    return;
                }
#endif
                // scalar
                for (auto blk {0uz}; blk < nblk; ++blk)
                {
                    const std::size_t jb = blk * Ns;
                    const std::size_t ob = blk * 2 * Ns;
                    for (auto k {0uz}; k < Ns; ++k)
                    {
                        const std::size_t j = jb + k;
                        const ST v0r = sr[j];
                        const ST v0i = si[j];
                        const ST ar = sr[j + half];
                        const ST ai = si[j + half];
                        const ST wr = twR[off + k];
                        const ST wi = inverse ? -twI[off + k] : twI[off + k];
                        const ST v1r = ar * wr - ai * wi;
                        const ST v1i = ar * wi + ai * wr;
                        dr[ob + k] = v0r + v1r; di[ob + k] = v0i + v1i;
                        dr[ob + k + Ns] = v0r - v1r; di[ob + k + Ns] = v0i - v1i;
                    }
                }
            }

            template <bool inverse>
            void stage4(const ST *__restrict sr, const ST *__restrict si, ST *__restrict dr, ST *__restrict di,
                        const std::size_t Ns,    const std::size_t off)
            {
                const std::size_t q = N >> 2;                // N/4 read stride
                const std::size_t nblk = q / Ns;
                const ST *w1R = tw4Re.data() + off;          const ST *w1I = tw4Im.data() + off;
                const ST *w2R = tw4Re.data() + off + Ns;     const ST *w2I = tw4Im.data() + off + Ns;
                const ST *w3R = tw4Re.data() + off + 2 * Ns; const ST *w3I = tw4Im.data() + off + 2 * Ns;
#ifdef QUOTIENT_FFT_XSIMD_INCLUDED_
                using B = xsimd::batch<ST>;
                constexpr std::size_t W = B::size;
                if (Ns >= W)
                {
                    for (auto blk {0uz}; blk < nblk; ++blk)
                    {
                        const std::size_t jb = blk * Ns;
                        const std::size_t ob = blk * 4 * Ns;
                        for (auto k {0uz}; k < Ns; k += W)
                        {
                            const std::size_t j = jb + k;
                            const B a0r = B::load_aligned(sr + j),         a0i = B::load_aligned(si + j);
                            const B a1r = B::load_aligned(sr + j + q),     a1i = B::load_aligned(si + j + q);
                            const B a2r = B::load_aligned(sr + j + 2 * q), a2i = B::load_aligned(si + j + 2 * q);
                            const B a3r = B::load_aligned(sr + j + 3 * q), a3i = B::load_aligned(si + j + 3 * q);
                            B w1r = B::load_unaligned(w1R + k), w1i = B::load_unaligned(w1I + k);
                            B w2r = B::load_unaligned(w2R + k), w2i = B::load_unaligned(w2I + k);
                            B w3r = B::load_unaligned(w3R + k), w3i = B::load_unaligned(w3I + k);
                            if constexpr (inverse) { w1i = -w1i; w2i = -w2i; w3i = -w3i; }
                            const B v1r = xsimd::fnma(a1i, w1i, a1r * w1r), v1i = xsimd::fma(a1i, w1r, a1r * w1i);
                            const B v2r = xsimd::fnma(a2i, w2i, a2r * w2r), v2i = xsimd::fma(a2i, w2r, a2r * w2i);
                            const B v3r = xsimd::fnma(a3i, w3i, a3r * w3r), v3i = xsimd::fma(a3i, w3r, a3r * w3i);
                            const B s02r = a0r + v2r, s02i = a0i + v2i;
                            const B d02r = a0r - v2r, d02i = a0i - v2i;
                            const B s13r = v1r + v3r, s13i = v1i + v3i;
                            const B d13r = v1r - v3r, d13i = v1i - v3i;
                            (s02r + s13r).store_aligned(dr + ob + k);
                            (s02i + s13i).store_aligned(di + ob + k);
                            (s02r - s13r).store_aligned(dr + ob + k + 2 * Ns);
                            (s02i - s13i).store_aligned(di + ob + k + 2 * Ns);
                            if constexpr (!inverse)
                            {
                                (d02r + d13i).store_aligned(dr + ob + k + Ns);
                                (d02i - d13r).store_aligned(di + ob + k + Ns);
                                (d02r - d13i).store_aligned(dr + ob + k + 3 * Ns);
                                (d02i + d13r).store_aligned(di + ob + k + 3 * Ns);
                            }
                            else
                            {
                                (d02r - d13i).store_aligned(dr + ob + k + Ns);
                                (d02i + d13r).store_aligned(di + ob + k + Ns);
                                (d02r + d13i).store_aligned(dr + ob + k + 3 * Ns);
                                (d02i - d13r).store_aligned(di + ob + k + 3 * Ns);
                            }
                        }
                    }
                    return;
                }
#endif
                // scalar
                for (auto blk {0uz}; blk < nblk; ++blk)
                {
                    const std::size_t jb = blk * Ns;
                    const std::size_t ob = blk * 4 * Ns;
                    for (auto k {0uz}; k < Ns; ++k)
                    {
                        const std::size_t j = jb + k;
                        const ST a0r = sr[j];
                        const ST a0i = si[j];
                        const ST a1r = sr[j + q];
                        const ST a1i = si[j + q];
                        const ST a2r = sr[j + 2 * q];
                        const ST a2i = si[j + 2 * q];
                        const ST a3r = sr[j + 3 * q];
                        const ST a3i = si[j + 3 * q];
                        const ST w1r = w1R[k];
                        const ST w1i = inverse ? -w1I[k] : w1I[k];
                        const ST w2r = w2R[k];
                        const ST w2i = inverse ? -w2I[k] : w2I[k];
                        const ST w3r = w3R[k];
                        const ST w3i = inverse ? -w3I[k] : w3I[k];
                        const ST v1r = a1r * w1r - a1i * w1i;
                        const ST v1i = a1r * w1i + a1i * w1r;
                        const ST v2r = a2r * w2r - a2i * w2i;
                        const ST v2i = a2r * w2i + a2i * w2r;
                        const ST v3r = a3r * w3r - a3i * w3i;
                        const ST v3i = a3r * w3i + a3i * w3r;
                        const ST s02r = a0r + v2r;
                        const ST s02i = a0i + v2i;
                        const ST d02r = a0r - v2r;
                        const ST d02i = a0i - v2i;
                        const ST s13r = v1r + v3r;
                        const ST s13i = v1i + v3i;
                        const ST d13r = v1r - v3r;
                        const ST d13i = v1i - v3i;
                        dr[ob + k] = s02r + s13r;
                        di[ob + k] = s02i + s13i;
                        dr[ob + k + 2 * Ns] = s02r - s13r;
                        di[ob + k + 2 * Ns] = s02i - s13i;
                        if constexpr (!inverse)
                        {
                            dr[ob + k + Ns] = d02r + d13i;
                            di[ob + k + Ns] = d02i - d13r;
                            dr[ob + k + 3 * Ns] = d02r - d13i;
                            di[ob + k + 3 * Ns] = d02i + d13r;
                        }
                        else
                        {
                            dr[ob + k + Ns] = d02r - d13i;
                            di[ob + k + Ns] = d02i + d13r;
                            dr[ob + k + 3 * Ns] = d02r + d13i;
                            di[ob + k + 3 * Ns] = d02i - d13r;
                        }
                    }
                }
            }
        };
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
        detail::SOAPow2FFT<ST> soa;              // power-of-two SoA/SIMD fast path
        bool usePow2 = false;

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
            for (auto k {0uz}; k < quarter; ++k)
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
    FFT(std::size_t) -> FFT<>;

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
            assert(std::has_single_bit(half) && "RealFFT half-size (n/2) must be a power of two!");
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

            for (auto k {0uz}; k < half; ++k) packed[k] = Complex(in[2 * k], in[2 * k + 1]);

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
            complexEngine.inv(spectrum.data(), packed.data());
            for (auto k {0uz}; k < half; ++k)
            {
                out[2*k] = packed[k].real();
                out[2*k+1] = packed[k].imag();
            }
        }

    private:
        static constexpr double twoPI = 2.0 * std::numbers::pi_v<double>;
        detail::SOAPow2FFT<ST> complexEngine;
        std::vector<Complex> packed;
        std::vector<Complex> spectrum;
        std::vector<Complex> twiddles;
    };
    RealFFT(std::size_t) -> RealFFT<>;
}
