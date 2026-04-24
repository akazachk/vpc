#pragma once

#include <ctime>
#include <string>
#include <vector>

#include "analysis.hpp"
#include "CglVPC.hpp"
#include "VPCParameters.hpp"

namespace main_startup {

int startUp(
    int argc,
    char** argv,
    VPCParametersNamespace::VPCParameters& params,
    CglVPC::ExitReason& exitReason,
    SummaryBoundInfo& boundInfo,
    std::time_t& start_time_t,
    char* start_time_string,
    size_t start_time_string_size,
    char* hostname,
    int hostname_max,
    std::string& cpu_model,
    int& cpu_id,
    std::string& dir,
    std::string& filename_stub,
    std::string& instname,
    std::string& in_file_ext,
    const std::vector<std::string>& overallTimeStatsName);

} // namespace main_startup
