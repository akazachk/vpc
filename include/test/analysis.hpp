// Name:     analysis.hpp
// Author:   A. M. Kazachkov
// Date:     2019-Mar-04
//-----------------------------------------------------------------------------
#pragma once

#include <string>
#include <vector>
#include <limits>

class OsiCuts;
struct VPCParameters;
struct SummaryBBInfo;

struct SummaryBoundInfo {
  double lp_obj = std::numeric_limits<double>::max();
  double best_disj_obj = std::numeric_limits<double>::lowest();
  double worst_disj_obj = std::numeric_limits<double>::lowest();
  double ip_obj = std::numeric_limits<double>::max();
  double gmic_obj = std::numeric_limits<double>::max();
  double lpc_obj = std::numeric_limits<double>::max();
  double vpc_obj = std::numeric_limits<double>::max();
  double gmic_vpc_obj = std::numeric_limits<double>::max();
  int num_gmic = 0, num_lpc = 0, num_vpc = 0;
}; /* SummaryBoundInfo */

void printHeader(const VPCParameters& params,
    const std::vector<std::string>& time_name,
    const char SEP = ',');
void printBoundAndGapInfo(const SummaryBoundInfo& boundInfo, FILE* logfile,
    const char SEP = ',');
void printBBInfo(const SummaryBBInfo& info_nocuts,
    const SummaryBBInfo& info_mycuts,
//    const SummaryBBInfo& info_allcuts,
    FILE *logfile, const int amountToPrint, const char SEP = ',');
void printBBInfo(const std::vector<SummaryBBInfo>& info, FILE* myfile, const bool print_blanks = false,
    const char SEP = ',');
//void printBBInfo(const SummaryBBInfo& info_mycuts,
//    const SummaryBBInfo& info_allcuts, FILE *myfile, const int amountToPrint,
//    const char SEP = ',');

void analyzeStrength(const VPCParameters& params,
    const SummaryBoundInfo& boundInfo, const OsiCuts* const vpc,
    std::string& output);
void analyzeBB(const VPCParameters& params, SummaryBBInfo& info_nocuts,
    SummaryBBInfo& info_mycuts, SummaryBBInfo& info_allcuts, std::string& output);
