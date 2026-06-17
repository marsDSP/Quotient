#pragma once

#include <string>

// Chart-data export entry point.
namespace quotient_charts
{
    // Returns true on success. Never throws; on I/O failure it returns false.
    bool exportChartData(const std::string& outDir);
}
