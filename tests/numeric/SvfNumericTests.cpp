#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "SvfReference.h"
#include "SvfChartData.h"

#include <dsp/filters/StateVariable.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>

namespace
{
    using Svf = MarsDSP::Filters::TwoPoleSVF;
    using OnePole = MarsDSP::Filters::OnePoleTPT;
    using RefType = svfref::Type;

    Svf::SVFType toProd(RefType t)
    {
        switch (t)
        {
            case RefType::LowPass: return Svf::SVFType::LowPass;
            case RefType::HighPass: return Svf::SVFType::HighPass;
            case RefType::BandPass: return Svf::SVFType::BandPass;
            case RefType::Notch: return Svf::SVFType::Notch;
            case RefType::AllPass: return Svf::SVFType::AllPass;
            case RefType::Bell: return Svf::SVFType::Bell;
            case RefType::LowShelf: return Svf::SVFType::LowShelf;
            case RefType::HighShelf: return Svf::SVFType::HighShelf;
        }
        return Svf::SVFType::LowPass;
    }

    const char *typeName(RefType t)
    {
        switch (t)
        {
            case RefType::LowPass: return "lowpass";
            case RefType::HighPass: return "highpass";
            case RefType::BandPass: return "bandpass";
            case RefType::Notch: return "notch";
            case RefType::AllPass: return "allpass";
            case RefType::Bell: return "bell";
            case RefType::LowShelf: return "lowshelf";
            case RefType::HighShelf: return "highshelf";
        }
        return "?";
    }

    Svf makeProd(RefType t, double fs, double fc, double Q, double gainDB)
    {
        Svf f;
        f.setParams(toProd(t), fs, fc, Q, gainDB);
        return f;
    }

    bool close(double a, double b, double absTol, double relTol)
    {
        const double diff = std::abs(a - b);
        return diff <= absTol
               || diff <= relTol * std::max(std::abs(a), std::abs(b));
    }

    std::vector<double> logspace(double lo, double hi, int n)
    {
        std::vector<double> v(static_cast<std::size_t>(n));
        const double ratio = std::log(hi / lo);
        for (int i = 0; i < n; ++i)
            v[static_cast<std::size_t>(i)] =
                    lo * std::exp(ratio * (static_cast<double>(i) / static_cast<double>(n - 1)));
        return v;
    }

    std::complex<double> measureProd(RefType t, double fs, double fc, double Q,
                                     double gainDB, double f, int n)
    {
        Svf filt = makeProd(t, fs, fc, Q, gainDB);
        auto step = [&filt](double x) { return filt.processSample<double>(x); };
        return svfref::measureResponse(step, f, fs, n);
    }
}

TEST_CASE (
"SVF magnitude() matches the paper transfer function"
,
"[svf][paper][magnitude]"
)
{
    const double fs = 48000.0;
    const std::vector<double> freqs = logspace(20.0, 0.45 * fs, 200);

    struct Case { RefType type; double fc; double Q; double gainDB; };
    const std::vector<Case> cases = {
        { RefType::LowPass,   1000.0, 0.7071, 0.0 },
        { RefType::LowPass,    220.0, 4.0,    0.0 },
        { RefType::HighPass,  1000.0, 0.7071, 0.0 },
        { RefType::HighPass,  3000.0, 2.0,    0.0 },
        { RefType::Notch,     1000.0, 1.0,    0.0 },
        { RefType::AllPass,   1000.0, 0.7071, 0.0 },
        { RefType::AllPass,    500.0, 3.0,    0.0 },
        { RefType::Bell,      1000.0, 2.0,   12.0 },
        { RefType::Bell,      1000.0, 2.0,  -12.0 },
        { RefType::LowShelf,   800.0, 0.7071, 9.0 },
        { RefType::LowShelf,   800.0, 0.7071,-9.0 },
        { RefType::HighShelf, 1200.0, 0.7071, 9.0 },
        { RefType::HighShelf, 1200.0, 0.7071,-9.0 },
    };

    for (const auto& c : cases)
    {
        const Svf prod = makeProd(c.type, fs, c.fc, c.Q, c.gainDB);
        const svfref::Coeffs ref = svfref::paperCoeffs(c.type, fs, c.fc, c.Q, c.gainDB);

        for (double f : freqs)
        {
            const double mc = prod.magnitude(f, fs);
            const double mp = svfref::paperMagnitude(ref, f, fs);
            CAPTURE(typeName(c.type), c.fc, c.Q, c.gainDB, f, mc, mp);
            REQUIRE(close(mc, mp, 1e-9, 1e-7));
        }
    }
}

TEST_CASE (
"SVF processSample() matches the paper bounded algorithm"
,
"[svf][paper][tick]"
)
{
    const double fs = 48000.0;

    struct Case { RefType type; double fc; double Q; double gainDB; };
    const std::vector<Case> cases = {
        { RefType::LowPass,   1000.0, 0.7071, 0.0 },
        { RefType::HighPass,   600.0, 2.0,    0.0 },
        { RefType::Notch,     1500.0, 1.5,    0.0 },
        { RefType::AllPass,    900.0, 3.0,    0.0 },
        { RefType::Bell,       750.0, 2.0,   10.0 },
        { RefType::LowShelf,   400.0, 0.7071, 8.0 },
        { RefType::HighShelf, 2000.0, 0.7071,-8.0 },
    };

    for (const auto& c : cases)
    {
        Svf prod = makeProd(c.type, fs, c.fc, c.Q, c.gainDB);
        svfref::PaperSvf ref;
        ref.setParams(c.type, fs, c.fc, c.Q, c.gainDB);

        double phase = 0.0;
        for (int n = 0; n < 4096; ++n)
        {
            const double saw = 2.0 * (phase - std::floor(phase + 0.5));
            const double x = 0.6 * saw + 0.3 * std::sin(0.01 * n) + ((n % 97 == 0) ? 1.0 : 0.0);
            phase += 437.0 / fs;

            const double yc = prod.processSample<double>(x);
            const double yr = ref.tick(x);
            CAPTURE(typeName(c.type), n, x, yc, yr);
            REQUIRE(close(yc, yr, 1e-9, 1e-9));
        }
    }
}

TEST_CASE (
"BandPass is normalised to 0 dB peak; paper band peaks at Q"
,
"[svf][paper][bandpass]"
)
{
    const double fs = 48000.0;
    const double fc = 1000.0;

    for (double Q : { 0.5, 1.0, 2.0, 4.0, 8.0 })
    {
        const Svf prod = makeProd(RefType::BandPass, fs, fc, Q, 0.0);
        const svfref::Coeffs ref = svfref::paperCoeffs(RefType::BandPass, fs, fc, Q, 0.0);

        const double peak = prod.magnitude(fc, fs);
        CAPTURE(Q, peak);
        REQUIRE(close(peak, 1.0, 1e-6, 1e-6));

        REQUIRE(close(svfref::paperMagnitude(ref, fc, fs), Q, 1e-6, 1e-6));

        for (double f : logspace(20.0, 0.45 * fs, 100))
        {
            const double mc = prod.magnitude(f, fs);
            const double mp = svfref::paperMagnitude(ref, f, fs);
            CAPTURE(Q, f, mc, mp);
            REQUIRE(close(mc * Q, mp, 1e-7, 1e-6));
        }
    }
}

TEST_CASE (
"SVF measured impulse response matches the analytic response"
,
"[svf][measure]"
)
{
    const double fs = 48000.0;
    const int    N  = 1 << 16;

    struct Case { RefType type; double fc; double Q; double gainDB; };
    const std::vector<Case> cases = {
        { RefType::LowPass,   1000.0, 0.7071, 0.0 },
        { RefType::HighPass,  1000.0, 1.5,    0.0 },
        { RefType::BandPass,   800.0, 2.0,    0.0 },
        { RefType::AllPass,    900.0, 1.0,    0.0 },
        { RefType::Bell,      1000.0, 2.0,    9.0 },
        { RefType::LowShelf,   600.0, 0.7071, 9.0 },
        { RefType::HighShelf, 1500.0, 0.7071, 9.0 },
    };

    for (const auto& c : cases)
    {
        const Svf prod = makeProd(c.type, fs, c.fc, c.Q, c.gainDB);
        for (double f : { 50.0, 200.0, 1000.0, 4000.0, 12000.0 })
        {
            const double measured = std::abs(measureProd(c.type, fs, c.fc, c.Q, c.gainDB, f, N));
            const double analytic = prod.magnitude(f, fs);
            CAPTURE(typeName(c.type), c.fc, c.Q, c.gainDB, f, measured, analytic);
            REQUIRE(close(measured, analytic, 5e-3, 5e-3));
        }
    }
}

TEST_CASE (
"SVF characteristic-frequency behaviour"
,
"[svf][points]"
)
{
    const double fs = 48000.0;
    const double fc = 1000.0;
    const double lowF = 20.0;
    const double hiF  = 0.45 * fs;

    SECTION("low-pass passes DC, rejects highs")
    {
        const Svf f = makeProd(RefType::LowPass, fs, fc, 0.7071, 0.0);
        REQUIRE(close(f.magnitude(lowF, fs), 1.0, 1e-3, 1e-3));
        REQUIRE(f.magnitude(hiF, fs) < 0.05);
    }
    SECTION("high-pass rejects DC, passes highs")
    {
        const Svf f = makeProd(RefType::HighPass, fs, fc, 0.7071, 0.0);
        REQUIRE(f.magnitude(lowF, fs) < 0.05);
        REQUIRE(close(f.magnitude(hiF, fs), 1.0, 1e-2, 1e-2));
    }
    SECTION("notch nulls at the cutoff")
    {
        const Svf f = makeProd(RefType::Notch, fs, fc, 1.0, 0.0);
        REQUIRE(20.0 * std::log10(f.magnitude(fc, fs)) < -60.0);
        REQUIRE(close(f.magnitude(lowF, fs), 1.0, 1e-2, 1e-2));
        REQUIRE(close(f.magnitude(hiF, fs), 1.0, 1e-2, 1e-2));
    }
    SECTION("all-pass is unity magnitude everywhere")
    {
        const Svf f = makeProd(RefType::AllPass, fs, fc, 3.0, 0.0);
        for (double freq : logspace(lowF, hiF, 64))
        {
            CAPTURE(freq);
            REQUIRE(close(f.magnitude(freq, fs), 1.0, 1e-9, 1e-9));
        }
    }
    SECTION("bell hits its gain at cutoff and is flat far away")
    {
        for (double g : { 12.0, 6.0, -6.0, -12.0 })
        {
            const Svf f = makeProd(RefType::Bell, fs, fc, 3.0, g);
            CAPTURE(g);
            REQUIRE(close(20.0 * std::log10(f.magnitude(fc, fs)), g, 1e-2, 1e-3));
            REQUIRE(close(f.magnitude(lowF, fs), 1.0, 5e-2, 5e-2));
            REQUIRE(close(f.magnitude(hiF, fs), 1.0, 5e-2, 5e-2));
        }
    }
    SECTION("low shelf lifts DC, leaves highs alone")
    {
        const double g = 12.0;
        const Svf f = makeProd(RefType::LowShelf, fs, fc, 0.7071, g);
        REQUIRE(close(20.0 * std::log10(f.magnitude(lowF, fs)), g, 2e-1, 1e-2));
        REQUIRE(close(20.0 * std::log10(f.magnitude(hiF, fs)), 0.0, 2e-1, 1e-2));
    }
    SECTION("high shelf lifts highs, leaves DC alone")
    {
        const double g = 12.0;
        const Svf f = makeProd(RefType::HighShelf, fs, fc, 0.7071, g);
        REQUIRE(close(20.0 * std::log10(f.magnitude(lowF, fs)), 0.0, 2e-1, 1e-2));
        REQUIRE(close(20.0 * std::log10(f.magnitude(hiF, fs)), g, 3e-1, 1e-2));
    }
}

TEST_CASE (
"SVF stays stable for high-Q, broadband input"
,
"[svf][stability]"
)
{
    const double fs = 48000.0;
    const RefType types[] = { RefType::LowPass, RefType::HighPass, RefType::BandPass,
                              RefType::Notch, RefType::AllPass, RefType::Bell,
                              RefType::LowShelf, RefType::HighShelf };

    for (RefType t : types)
    {
        Svf f = makeProd(t, fs, 500.0, 20.0, 12.0); // intentionally high Q
        unsigned int seed = 12345u;
        for (int n = 0; n < 200000; ++n)
        {
            seed = seed * 1664525u + 1013904223u;
            const double noise = (static_cast<double>(seed) / 4294967295.0) * 2.0 - 1.0;
            const double y = f.processSample<double>(noise);
            CAPTURE(typeName(t), n, y);
            REQUIRE(std::isfinite(y));
            REQUIRE(std::abs(y) < 1.0e3);
        }

        Svf g = makeProd(t, fs, 500.0, 20.0, 12.0);
        double tail = 0.0;
        for (int n = 0; n < 200000; ++n)
        {
            const double y = g.processSample<double>(n == 0 ? 1.0 : 0.0);
            if (n > 190000) tail += std::abs(y);
        }
        CAPTURE(typeName(t), tail);
        REQUIRE(tail < 1e-6);
    }
}

TEST_CASE (
"OnePoleTPT basic low/high-pass behaviour"
,
"[onepole][sanity]"
)
{
    const double fs = 48000.0;
    const double fc = 1000.0;

    OnePole lp;
    lp.setParams(OnePole::Type::LowPass, fs, fc);
    REQUIRE(close(lp.magnitude(20.0, fs), 1.0, 1e-3, 1e-3));
    REQUIRE(close(lp.magnitude(fc, fs), 1.0 / std::sqrt(2.0), 1e-9, 1e-9)); // -3 dB at cutoff
    REQUIRE(lp.magnitude(0.45 * fs, fs) < 0.1);

    OnePole hp;
    hp.setParams(OnePole::Type::HighPass, fs, fc);
    REQUIRE(hp.magnitude(20.0, fs) < 0.1);
    REQUIRE(close(hp.magnitude(fc, fs), 1.0 / std::sqrt(2.0), 1e-9, 1e-9));
    REQUIRE(close(hp.magnitude(0.45 * fs, fs), 1.0, 1e-2, 1e-2));

    for (double f : logspace(20.0, 0.45 * fs, 32))
    {
        const double a = lp.magnitude(f, fs);
        const double b = hp.magnitude(f, fs);
        CAPTURE(f, a, b);
        REQUIRE(close(a * a + b * b, 1.0, 1e-9, 1e-9));
    }
}

TEST_CASE (
"TiltShelf tilts the spectrum symmetrically"
,
"[svf][tilt][sanity]"
)
{
    const double fs = 48000.0;
    const double fc = 1000.0;
    const double gainDB = 12.0;
    const double A = svfref::gainToA(gainDB);

    Svf f;
    f.setParams(Svf::SVFType::TiltShelf, fs, fc, 0.7071, gainDB);

    const double lowGain = f.magnitude(20.0, fs);
    const double hiGain  = f.magnitude(0.45 * fs, fs);
    CAPTURE(A, lowGain, hiGain);
    REQUIRE(close(lowGain, 1.0 / A, 5e-2, 5e-2));
    REQUIRE(close(hiGain, A, 5e-2, 5e-2));
    REQUIRE(hiGain > lowGain);
}

namespace
{
    bool writeCsv(const std::filesystem::path &path,
                  const std::vector<std::string> &header,
                  const std::vector<std::vector<double> > &cols)
    {
        std::ofstream os(path);
        if (!os) return false;

        for (std::size_t c = 0; c < header.size(); ++c)
            os << header[c] << (c + 1 < header.size() ? "," : "\n");

        const std::size_t rows = cols.empty() ? 0 : cols.front().size();
        os << std::setprecision(10);
        for (std::size_t r = 0; r < rows; ++r)
        {
            for (std::size_t c = 0; c < cols.size(); ++c)
                os << cols[c][r] << (c + 1 < cols.size() ? "," : "\n");
        }
        return static_cast<bool>(os);
    }

    std::vector<double> magnitudeDbSweep(RefType t, double fs, double fc, double Q,
                                         double gainDB, const std::vector<double> &freqs)
    {
        const Svf f = makeProd(t, fs, fc, Q, gainDB);
        std::vector<double> out;
        out.reserve(freqs.size());
        for (double freq: freqs)
            out.push_back(20.0 * std::log10(f.magnitude(freq, fs)));
        return out;
    }
}

bool quotient_charts::exportChartData(const std::string &outDir)
{
    namespace fs_ns = std::filesystem;
    std::error_code ec;
    fs_ns::create_directories(outDir, ec);
    if (ec) return false;

    const fs_ns::path dir(outDir);
    const double fs = 48000.0;
    const double fc = 1000.0;
    const std::vector<double> freqs = logspace(20.0, 0.49 * fs, 512);
    bool ok = true;

    // (a) Core responses overlaid (lp/hp/bp/notch/allpass), Q = 1/sqrt(2)
    ok &= writeCsv(dir / "mag_core.csv",
                   {"freq", "lowpass", "highpass", "bandpass", "notch", "allpass"},
                   {
                       freqs,
                       magnitudeDbSweep(RefType::LowPass, fs, fc, 0.7071, 0.0, freqs),
                       magnitudeDbSweep(RefType::HighPass, fs, fc, 0.7071, 0.0, freqs),
                       magnitudeDbSweep(RefType::BandPass, fs, fc, 0.7071, 0.0, freqs),
                       magnitudeDbSweep(RefType::Notch, fs, fc, 0.7071, 0.0, freqs),
                       magnitudeDbSweep(RefType::AllPass, fs, fc, 0.7071, 0.0, freqs)
                   });

    // (b) Low-pass resonance for a family of Q values
    ok &= writeCsv(dir / "mag_lowpass_resonance.csv",
                   {"freq", "Q=0.5", "Q=1", "Q=2", "Q=4", "Q=8"},
                   {
                       freqs,
                       magnitudeDbSweep(RefType::LowPass, fs, fc, 0.5, 0.0, freqs),
                       magnitudeDbSweep(RefType::LowPass, fs, fc, 1.0, 0.0, freqs),
                       magnitudeDbSweep(RefType::LowPass, fs, fc, 2.0, 0.0, freqs),
                       magnitudeDbSweep(RefType::LowPass, fs, fc, 4.0, 0.0, freqs),
                       magnitudeDbSweep(RefType::LowPass, fs, fc, 8.0, 0.0, freqs)
                   });

    // (c) Bell boost/cut
    ok &= writeCsv(dir / "mag_bell.csv",
                   {"freq", "+12dB", "+6dB", "-6dB", "-12dB"},
                   {
                       freqs,
                       magnitudeDbSweep(RefType::Bell, fs, fc, 2.0, 12.0, freqs),
                       magnitudeDbSweep(RefType::Bell, fs, fc, 2.0, 6.0, freqs),
                       magnitudeDbSweep(RefType::Bell, fs, fc, 2.0, -6.0, freqs),
                       magnitudeDbSweep(RefType::Bell, fs, fc, 2.0, -12.0, freqs)
                   });

    // (d) Shelving filters, +/-12 dB
    ok &= writeCsv(dir / "mag_shelf.csv",
                   {"freq", "lowshelf+12", "lowshelf-12", "highshelf+12", "highshelf-12"},
                   {
                       freqs,
                       magnitudeDbSweep(RefType::LowShelf, fs, fc, 0.7071, 12.0, freqs),
                       magnitudeDbSweep(RefType::LowShelf, fs, fc, 0.7071, -12.0, freqs),
                       magnitudeDbSweep(RefType::HighShelf, fs, fc, 0.7071, 12.0, freqs),
                       magnitudeDbSweep(RefType::HighShelf, fs, fc, 0.7071, -12.0, freqs)
                   });

    // (e) Verification overlay code magnitude() vs paper transfer function
    {
        const svfref::Coeffs lpRef = svfref::paperCoeffs(RefType::LowPass, fs, fc, 2.0, 0.0);
        const Svf lp = makeProd(RefType::LowPass, fs, fc, 2.0, 0.0);
        std::vector<double> codeDb, paperDb;
        codeDb.reserve(freqs.size());
        paperDb.reserve(freqs.size());
        for (double f: freqs)
        {
            codeDb.push_back(20.0 * std::log10(lp.magnitude(f, fs)));
            paperDb.push_back(svfref::paperMagnitudeDb(lpRef, f, fs));
        }
        ok &= writeCsv(dir / "verify_lowpass.csv",
                       {"freq", "code_db", "paper_db"}, {freqs, codeDb, paperDb});
    }

    // (f) Verification overlay for bandpass; shows the 1/Q normalisation
    {
        const double Q = 4.0;
        const svfref::Coeffs bpRef = svfref::paperCoeffs(RefType::BandPass, fs, fc, Q, 0.0);
        const Svf bp = makeProd(RefType::BandPass, fs, fc, Q, 0.0);
        std::vector<double> codeDb, paperDb;
        codeDb.reserve(freqs.size());
        paperDb.reserve(freqs.size());
        for (double f: freqs)
        {
            codeDb.push_back(20.0 * std::log10(bp.magnitude(f, fs)));
            paperDb.push_back(svfref::paperMagnitudeDb(bpRef, f, fs));
        }
        ok &= writeCsv(dir / "verify_bandpass.csv",
                       {"freq", "code_normalised_db", "paper_band_db"}, {freqs, codeDb, paperDb});
    }

    // (g) All-pass phase: measured (production IR DFT) vs paper analytic
    {
        const svfref::Coeffs apRef = svfref::paperCoeffs(RefType::AllPass, fs, fc, 0.7071, 0.0);
        const int N = 1 << 16;
        std::vector<double> measDeg, paperDeg;
        measDeg.reserve(freqs.size());
        paperDeg.reserve(freqs.size());
        for (double f: freqs)
        {
            const auto meas = measureProd(RefType::AllPass, fs, fc, 0.7071, 0.0, f, N);
            measDeg.push_back(std::arg(meas) * 180.0 / svfref::kPi);
            paperDeg.push_back(std::arg(svfref::paperResponse(apRef, f, fs)) * 180.0 / svfref::kPi);
        }
        ok &= writeCsv(dir / "phase_allpass.csv",
                       {"freq", "measured_deg", "paper_deg"}, {freqs, measDeg, paperDeg});
    }

    // (h) Time-domain saw test (mirrors the paper's pp.9-10 figures)
    {
        const double fsTime = 44100.0;
        const double t1 = 0.005;
        const double sawHz = 500.0;
        const double fcTime = 1000.0;
        const int n = static_cast<int>(t1 * fsTime) + 1;

        std::vector<double> t(static_cast<std::size_t>(n)), in(static_cast<std::size_t>(n));
        const RefType outTypes[] = {
            RefType::LowPass, RefType::HighPass, RefType::BandPass,
            RefType::Notch, RefType::AllPass, RefType::Bell,
            RefType::LowShelf, RefType::HighShelf
        };
        std::vector<Svf> filters;
        for (RefType ty: outTypes)
            filters.push_back(makeProd(ty, fsTime, fcTime, 0.7071,
                                       (ty == RefType::Bell || ty == RefType::LowShelf
                                        || ty == RefType::HighShelf)
                                           ? 12.0
                                           : 0.0));
        std::vector<std::vector<double> > outs(std::size(outTypes),
                                               std::vector<double>(static_cast<std::size_t>(n)));

        for (int i = 0; i < n; ++i)
        {
            const double time = static_cast<double>(i) / fsTime;
            const double ph = sawHz * time;
            const double saw = 2.0 * (ph - std::floor(ph + 0.5)); // band-unlimited saw
            t[static_cast<std::size_t>(i)] = time;
            in[static_cast<std::size_t>(i)] = saw;
            for (std::size_t f = 0; f < filters.size(); ++f)
                outs[f][static_cast<std::size_t>(i)] = filters[f].processSample<double>(saw);
        }

        std::vector<std::string> header = {"t", "input"};
        std::vector<std::vector<double> > cols = {t, in};
        for (std::size_t f = 0; f < std::size(outTypes); ++f)
        {
            header.push_back(typeName(outTypes[f]));
            cols.push_back(outs[f]);
        }
        ok &= writeCsv(dir / "time_saw.csv", header, cols);
    }

    return ok;
}
