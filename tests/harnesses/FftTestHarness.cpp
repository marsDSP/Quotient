#include <catch2/catch_session.hpp>

#include "FftChartData.h"

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#ifndef QUOTIENT_DEFAULT_CHART_DIR
#define QUOTIENT_DEFAULT_CHART_DIR "charts_data"
#endif

namespace
{
    std::string resolveChartDir(const std::string &cliDir)
    {
        if (!cliDir.empty())
            return cliDir;
        if (const char *env = std::getenv("QUOTIENT_CHART_DATA_DIR"); env && *env)
            return env;
        return QUOTIENT_DEFAULT_CHART_DIR;
    }
} // namespace

int main(int argc, char *argv[])
{
    bool        exportCharts = false;
    std::string chartDir;

    std::vector<char *> forwarded;
    forwarded.reserve(static_cast<std::size_t>(argc));
    forwarded.push_back(argv[0]);

    for (int i = 1; i < argc; ++i)
    {
        const std::string arg = argv[i];
        if (arg == "--export-charts")
        {
            exportCharts = true;
        }
        else if (arg == "--chart-dir")
        {
            if (i + 1 < argc)
                chartDir = argv[++i];
            else
            {
                std::cerr << "error: --chart-dir requires a path argument\n";
                return 2;
            }
        }
        else if (arg.rfind("--chart-dir=", 0) == 0)
        {
            chartDir = arg.substr(std::strlen("--chart-dir="));
        }
        else
        {
            forwarded.push_back(argv[i]);
        }
    }

    Catch::Session session;
    const int cliResult =
        session.applyCommandLine(static_cast<int>(forwarded.size()), forwarded.data());
    if (cliResult != 0)
        return cliResult;

    const int testResult = session.run();

    if (exportCharts)
    {
        if (testResult != 0)
        {
            std::cerr << "skipping chart export: tests failed (" << testResult << ")\n";
            return testResult;
        }

        const std::string dir = resolveChartDir(chartDir);
        std::cout << "exporting FFT chart data to: " << dir << '\n';
        if (!quotient_fft_charts::exportChartData(dir))
        {
            std::cerr << "error: failed to write chart data to " << dir << '\n';
            return 3;
        }
        if (!quotient_fft_charts::exportSimdChartData(dir))
        {
            std::cerr << "error: failed to write SIMD chart data to " << dir << '\n';
            return 3;
        }
        std::cout << "chart data written.\n";
    }
    return testResult;
}
