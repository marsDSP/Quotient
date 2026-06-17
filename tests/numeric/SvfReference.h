#pragma once

// ============================================================================
//  SvfReference.h
//
//  An INDEPENDENT reference implementation of the state-variable filter
//  described in:
//
//      Andrew Simper (Cytomic), "Solving the continuous SVF equations using
//      trapezoidal integration and equivalent currents", 2013 (rev. 2016).
//      docs/papers/SvfLinearTrapOptimised2.pdf
//
//  This header deliberately does NOT include the production filter
//  (source/dsp/filters/StateVariable.h) and uses only <cmath>/<complex> so the
//  tests have a truly independent oracle to compare against. Everything here is
//  transcribed directly from the paper:
//
//    * coefficient table .............. "Algorithm for every response" (pp.31-32)
//    * per-sample update (tick) ....... "Algorithm ... bounded terms" (pp.5-6)
//    * discrete transfer functions .... "Transfer functions ..." (p.12)
//
//  NOTE on bandpass: the paper's generic `band` response uses m1 = 1, which
//  yields a band-pass whose peak gain equals Q. The production code normalises
//  this to a 0 dB peak by using m1 = 1/Q. We keep the *paper* value here so the
//  difference can be asserted explicitly in the tests.
// ============================================================================

#include <cmath>
#include <complex>
#include <functional>

namespace svfref
{
    inline constexpr double kPi = 3.14159265358979323846;

    enum class Type
    {
        LowPass,
        HighPass,
        BandPass,
        Notch,
        AllPass,
        Bell,
        LowShelf,
        HighShelf
    };

    struct Coeffs
    {
        double g  = 0.0; // tan(pi * fc / fs), possibly shelf-adjusted
        double k  = 1.0; // damping (1/Q, or 1/(Q*A) for bell)
        double m0 = 0.0; // input mix
        double m1 = 0.0; // band-pass (v1) mix
        double m2 = 0.0; // low-pass  (v2) mix
    };

    // Gain term A = 10^(gainDB/40), exactly as Power[10, gaindB/40] in the paper.
    inline double gainToA(double gainDB) { return std::pow(10.0, gainDB / 40.0); }

    // Coefficient table transcribed verbatim from the paper (pp.31-32).
    inline Coeffs paperCoeffs(Type t, double fs, double fc, double Q, double gainDB)
    {
        const double A  = gainToA(gainDB);
        const double gt = std::tan(kPi * fc / fs);
        const double kk = 1.0 / Q;

        Coeffs c;
        c.g = gt;
        c.k = kk;

        switch (t)
        {
            case Type::LowPass:                                   // case low
                c.m0 = 0.0; c.m1 = 0.0;            c.m2 = 1.0;            break;
            case Type::HighPass:                                  // case high
                c.m0 = 1.0; c.m1 = -kk;            c.m2 = -1.0;           break;
            case Type::BandPass:                                  // case band
                c.m0 = 0.0; c.m1 = 1.0;            c.m2 = 0.0;            break;
            case Type::Notch:                                     // case notch
                c.m0 = 1.0; c.m1 = -kk;            c.m2 = 0.0;            break;
            case Type::AllPass:                                   // case all
                c.m0 = 1.0; c.m1 = -2.0 * kk;      c.m2 = 0.0;            break;
            case Type::Bell:                                      // case bell
                c.k  = 1.0 / (Q * A);
                c.m0 = 1.0; c.m1 = c.k * (A * A - 1.0); c.m2 = 0.0;       break;
            case Type::LowShelf:                                  // case low shelf
                c.g  = gt / std::sqrt(A);
                c.m0 = 1.0; c.m1 = kk * (A - 1.0); c.m2 = (A * A - 1.0);  break;
            case Type::HighShelf:                                 // case high shelf
                c.g  = gt * std::sqrt(A);
                c.m0 = A * A; c.m1 = kk * (1.0 - A) * A; c.m2 = (1.0 - A * A); break;
        }
        return c;
    }

    // Independent transcription of the bounded trapezoidal algorithm (pp.5-6).
    struct PaperSvf
    {
        Coeffs c;
        double a1 = 1.0, a2 = 0.0, a3 = 0.0;
        double ic1eq = 0.0, ic2eq = 0.0;

        void setParams(Type t, double fs, double fc, double Q, double gainDB)
        {
            c  = paperCoeffs(t, fs, fc, Q, gainDB);
            a1 = 1.0 / (1.0 + c.g * (c.g + c.k));
            a2 = c.g * a1;
            a3 = c.g * a2;
        }

        void reset() { ic1eq = 0.0; ic2eq = 0.0; }

        double tick(double v0)
        {
            const double v3 = v0 - ic2eq;
            const double v1 = a1 * ic1eq + a2 * v3;
            const double v2 = ic2eq + a2 * ic1eq + a3 * v3;
            ic1eq = 2.0 * v1 - ic1eq;
            ic2eq = 2.0 * v2 - ic2eq;
            return c.m0 * v0 + c.m1 * v1 + c.m2 * v2;
        }
    };

    // Exact discrete transfer function H(e^{jw}).
    //
    // The Cytomic SVF is the bilinear transform (with cutoff pre-warp
    // g = tan(pi*fc/fs)) of the analogue prototype
    //
    //      H_a(s_n) = [ m0 s_n^2 + (m0 k + m1) s_n + (m0 + m2) ]
    //                 / [ s_n^2 + k s_n + 1 ]
    //
    // On the unit circle z = e^{jw} the bilinear map sends s_n -> j*Omega with
    //
    //      Omega = tan(pi*f/fs) / g.
    //
    // Substituting reproduces the paper's z-domain transfer functions (p.12)
    // exactly, so this is an analytic oracle (no truncation error).
    inline std::complex<double> paperResponse(const Coeffs& c, double f, double fs)
    {
        const double Omega = std::tan(kPi * f / fs) / c.g;
        const std::complex<double> s(0.0, Omega);
        const std::complex<double> num =
            c.m0 * s * s + (c.m0 * c.k + c.m1) * s + (c.m0 + c.m2);
        const std::complex<double> den = s * s + c.k * s + 1.0;
        return num / den;
    }

    inline double paperMagnitude(const Coeffs& c, double f, double fs)
    {
        return std::abs(paperResponse(c, f, fs));
    }

    inline double paperMagnitudeDb(const Coeffs& c, double f, double fs)
    {
        return 20.0 * std::log10(paperMagnitude(c, f, fs));
    }

    // ------------------------------------------------------------------------
    // Black-box frequency-response measurement.
    //
    // Drives an arbitrary per-sample callable with a unit impulse and computes
    // the single-bin DFT of the (truncated) impulse response at frequency f.
    // This measures the *time-domain* behaviour of whatever filter is supplied,
    // independent of any analytic magnitude helper, so it cross-checks that the
    // recurrence actually realises the transfer function above.
    // ------------------------------------------------------------------------
    template <class Step>
    std::complex<double> measureResponse(Step&& step, double f, double fs, int numSamples)
    {
        const double w = 2.0 * kPi * f / fs;
        std::complex<double> acc(0.0, 0.0);
        for (int n = 0; n < numSamples; ++n)
        {
            const double x = (n == 0) ? 1.0 : 0.0;
            const double y = step(x);
            acc += y * std::exp(std::complex<double>(0.0, -w * static_cast<double>(n)));
        }
        return acc;
    }
}
