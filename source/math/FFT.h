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
        using cpx = std::complex<double>;
        explicit FFT(const std::size_t size) : sz(size)
        {
            assert(std::has_single_bit(sz) && "FFT size must be a power of two!");
        }

        [[nodiscard]] std::size_t size() const
        {
            return sz;
        }

        [[nodiscard]] std::vector<cpx> fwd(const std::vector<cpx> &x) const
        {
            assert(x.size() == sz && "input size must match FFT size!");
            return transform<false>(x);
        }

        [[nodiscard]] std::vector<cpx> inv(const std::vector<cpx> &X) const
        {
            assert(X.size() == sz && "input size must match FFT size!");
            return transform<true>(X);
        }

    private:
        enum class StepType
        {
            generic, step2, step3, step4
        };

        struct Step
        {
            StepType type;
            std::size_t factor;
            std::size_t start;
            std::size_t innerReps;
            std::size_t outerReps;
            std::size_t twiddleIdx;
        };

        struct PermutationPair
        {
            std::size_t from;
            std::size_t to;
        };

        std::size_t sz;
        std::size_t N = 0;
        std::vector<std::size_t> factors;
        std::vector<cpx> twiddleVector;
        std::vector<PermutationPair> permutation;
        std::vector<Step> plan;

        static constexpr double twoPI = 2 * std::numbers::pi_v<double>;

        void setPlan()
        {
            factors.clear();
            {
                std::size_t s = N;
                std::size_t f = 2;
                while (s > 1)
                {
                    if (s % f == 0) { factors.push_back(f); s /= f; }
                    else if (f * f > s) { f = s; }
                    else { ++f; }
                }
            }

            plan.clear();
            twiddleVector.clear();
            if (N > 1) addPlanSteps(0, 0, N, 1);
            permutation.clear();
            permutation.push_back(PermutationPair{0, 0});
            std::size_t idxLow = 0;
            std::size_t idxHigh = factors.size();
            std::size_t inLow = N;
            std::size_t outLow = 1;
            std::size_t inHigh = 1;
            std::size_t outHigh = N;

            while (outLow * inHigh < N)
            {
                std::size_t f;
                std::size_t inStep;
                std::size_t outStep;
                if (outLow <= inHigh)
                {
                    f = factors[idxLow++];
                    inStep = (inLow /= f);
                    outStep = outLow;
                    outLow *= f;
                }
                else
                {
                    f = factors[--idxHigh];
                    inStep = inHigh;
                    inHigh *= f;
                    outStep = (outHigh /= f);
                }
                const std::size_t old = permutation.size();
                for (auto i {1uz}; i < f; ++i)
                    for (auto j {0uz}; j < old; ++j)
                    {
                        PermutationPair p = permutation[j];
                        p.from += i * inStep;
                        p.to += i * outStep;
                        permutation.push_back(p);
                    }
            }
        }

        void addPlanSteps(std::size_t factorIndex,
                          const std::size_t start,
                          const std::size_t length,
                          const std::size_t repeats)
        {
            if (factorIndex >= factors.size()) return;
            std::size_t factor = factors[factorIndex];
            if (factorIndex + 1 < factors.size()
                && factors[factorIndex] == 2
                && factors[factorIndex + 1] == 2)
            {
                ++factorIndex;
                factor = 4;
            }

            const std::size_t subLength = length / factor;
            Step step { StepType::generic, factor, start, subLength, repeats, twiddleVector.size() };
            if (factor == 2) step.type = StepType::step2;
            else if (factor == 3) step.type = StepType::step3;
            else if (factor == 4) step.type = StepType::step4;

            bool found = false;
            for (const Step &e : plan)
                if (e.factor == step.factor && e.innerReps == step.innerReps)
                {
                    step.twiddleIdx = e.twiddleIdx;
                    found = true;
                    break;
                }
            if (!found)
                for (auto i {0uz}; i < subLength; ++i)
                    for (auto f {0uz}; f < factor; ++f)
                    {
                        const double phase = twoPI * static_cast<double>(i) *
                                                     static_cast<double>(f) /
                                                     static_cast<double>(length);
                        twiddleVector.push_back(cpx(static_cast<ST>(std::cos(phase)),
                                                    static_cast<ST>(-std::sin(phase))));
                    }
            addPlanSteps(factorIndex + 1, start, subLength, repeats * factor);
            plan.push_back(step);
        }

        template <bool inverse>
        void run(const cpx *in, cpx *out) const
        {
            for (const PermutationPair &p : permutation)
                out[p.from] = in[p.to];
            for (const Step &s : plan)
                switch (s.type)
                {
                    case StepType::step2: fftStep2<inverse>(out + s.start, s); break;
                    case StepType::step3: fftStep3<inverse>(out + s.start, s); break;
                    case StepType::step4: fftStep4<inverse>(out + s.start, s); break;
                    case StepType::generic: fftStepGeneric<inverse>(out + s.start, s); break;
                    default:;
                }
        }

        template <bool inverse>
        static std::vector<cpx> transform(std::vector<cpx> x)
        {
            const std::size_t N = x.size();
            if (N == 1) return x;                       // 1-point DFT is identity
            if (N == 2)                                 // radix-2 butterfly base case
            {
                const cpx a = x[0];
                const cpx b = x[1];
                return {a + b, a - b};                  // W_2^0 = 1, W_2^1 = -1
            }

            const std::size_t M = N / 4;                // sub-transform length
            std::vector<cpx> x0(M);
            std::vector<cpx> x1(M);
            std::vector<cpx> x2(M);
            std::vector<cpx> x3(M);
            for (auto k{0uz}; k < M; ++k)
            {
                x0[k] = x[4 * k];
                x1[k] = x[4 * k + 1];
                x2[k] = x[4 * k + 2];
                x3[k] = x[4 * k + 3];
            }
            const std::vector<cpx> X0 = transform<inverse>(x0); // DFT of n%4==0
            const std::vector<cpx> X1 = transform<inverse>(x1); // DFT of n%4==1
            const std::vector<cpx> X2 = transform<inverse>(x2); // DFT of n%4==2
            const std::vector<cpx> X3 = transform<inverse>(x3); // DFT of n%4==3

            std::vector<cpx> result(N);
            for (auto k {0uz}; k < M; ++k)
            {
                const double angle = -twoPI * static_cast<double>(k) / static_cast<double>(N);
                const cpx W1 = std::polar(1.0, angle);
                const cpx W2 = detail::complexMul<false>(W1, W1);
                const cpx W3 = detail::complexMul<false>(W2, W1);

                const cpx t0 = X0[k];
                const cpx t1 = inverse ? detail::complexMul<true>(X1[k], W1)
                                       : detail::complexMul<false>(W1, X1[k]);
                const cpx t2 = inverse ? detail::complexMul<true>(X2[k], W2)
                                       : detail::complexMul<false>(W2, X2[k]);
                const cpx t3 = inverse ? detail::complexMul<true>(X3[k], W3)
                                       : detail::complexMul<false>(W3, X3[k]);

                const cpx sum02 = t0 + t2;
                const cpx dif02 = t0 - t2;
                const cpx sum13 = t1 + t3;
                const cpx dif13 = t1 - t3;

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