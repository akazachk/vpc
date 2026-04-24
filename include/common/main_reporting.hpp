#pragma once

#include <ctime>
#include <string>
#include <vector>

#include <OsiSolverInterface.hpp>

#include "analysis.hpp"
#include "CglVPC.hpp"
#include "TimeStats.hpp"
#include "VPCParameters.hpp"

namespace main_reporting {

int wrapUp(
    int retCode,
    int argc,
    char** argv,
    VPCParametersNamespace::VPCParameters& params,
    TimeStats& timer,
    const std::vector<std::string>& overallTimeStatsName,
    CglVPC::ExitReason exitReason,
    std::time_t start_time_t,
    const char* start_time_string,
    const char* hostname,
    const std::string& cpu_model,
    int cpu_id,
    const std::string& instname,
    const std::string& cut_output,
    const std::string& bb_output,
    const SummaryBoundInfo& boundInfo,
    const SummaryBBInfo& info_nocuts,
    const SummaryBBInfo& info_mycuts,
    const SummaryCutInfo& cutInfoGMICs,
    const SummaryCutInfo& cutInfo,
    const SummaryDisjunctionInfo& disjInfo,
    OsiSolverInterface*& solver,
    OsiSolverInterface*& origSolver,
    OsiSolverInterface*& GMICSolver,
    OsiSolverInterface*& VPCSolver,
    OsiSolverInterface*& allCutsSolver,
    const std::string& vpc_version_string,
    const std::string& cbc_version_string,
    const std::string& clp_version_string,
    const std::string& gurobi_version_string,
    const std::string& cplex_version_string);

} // namespace main_reporting
