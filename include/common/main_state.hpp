#pragma once

#include <ctime>
#include <string>
#include <vector>

#include <OsiSolverInterface.hpp>

#include "analysis.hpp"
#include "CglVPC.hpp"
#include "TimeStats.hpp"
#include "VPCParameters.hpp"

struct MainAppState {
  VPCParametersNamespace::VPCParameters& params;
  TimeStats& timer;
  const std::vector<std::string>& overallTimeStatsName;

  CglVPC::ExitReason& exitReason;
  SummaryBoundInfo& boundInfo;
  SummaryBBInfo& info_nocuts;
  SummaryBBInfo& info_mycuts;
  SummaryCutInfo& cutInfoGMICs;
  SummaryCutInfo& cutInfo;
  SummaryDisjunctionInfo& disjInfo;

  std::time_t& start_time_t;
  char* start_time_string;
  size_t start_time_string_size;

  char* hostname;
  int hostname_max;
  std::string& cpu_model;
  int& cpu_id;

  std::string& dir;
  std::string& filename_stub;
  std::string& instname;
  std::string& in_file_ext;

  std::string& cut_output;
  std::string& bb_output;

  OsiSolverInterface*& solver;
  OsiSolverInterface*& origSolver;
  OsiSolverInterface*& GMICSolver;
  OsiSolverInterface*& VPCSolver;
  OsiSolverInterface*& allCutsSolver;

  const std::string& vpc_version_string;
  const std::string& cbc_version_string;
  const std::string& clp_version_string;
  const std::string& gurobi_version_string;
  const std::string& cplex_version_string;
};
