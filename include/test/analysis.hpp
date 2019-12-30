// Name:     analysis.hpp
// Author:   A. M. Kazachkov
// Date:     2019-Mar-04
//-----------------------------------------------------------------------------
#pragma once

#include <string>
#include <vector>
#include <limits>

class OsiSolverInterface;
class OsiCuts;

#include "CglVPC.hpp"
namespace VPCParametersNamespace {
  struct VPCParameters;
}

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
  double all_cuts_obj = std::numeric_limits<double>::max();
  int num_gmic = 0, num_lpc = 0, num_vpc = 0;
}; /* SummaryBoundInfo */

struct SummaryDisjunctionInfo {
  double avg_num_terms = 0;
  double avg_density_prlp = 0.;
  double avg_num_rows_prlp = 0.;
  double avg_num_cols_prlp = 0.;
  double avg_num_points_prlp = 0.;
  double avg_num_rays_prlp = 0.;
  double avg_explored_nodes = 0.;
  double avg_pruned_nodes = 0.;
  double avg_min_depth = 0.;
  double avg_max_depth = 0.;
}; /* SummaryDisjunctionInfo */

struct SummaryCutInfo {
  int num_cuts = 0;
  int num_active = 0;
  int num_obj_tried = 0, num_failures = 0;
  int num_rounds = 0;
  int min_support = std::numeric_limits<int>::max();
  int max_support = 0;
  double avg_support = 0.;
  std::vector<CglVPC::CutType> cutType; // one entry per cut
  std::vector<CglVPC::ObjectiveType> objType; // one entry per cut

  std::vector<int> numCutsOfType;
  std::vector<int> numCutsFromHeur, numObjFromHeur, numFailsFromHeur, numActiveFromHeur;
  std::vector<int> numFails;
}; /* SummaryCutInfo */

void printHeader(const VPCParametersNamespace::VPCParameters& params,
    const std::vector<std::string>& time_name,
    const char SEP = ',');
void printBoundAndGapInfo(const SummaryBoundInfo& boundInfo, FILE* logfile,
    const char SEP = ',');
void printSummaryBBInfo(const std::vector<SummaryBBInfo>& info, FILE* myfile,
    const bool print_blanks = false, const char SEP = ',');
void printFullBBInfo(const std::vector<SummaryBBInfo>& info, FILE* myfile,
    const bool print_blanks = false, const char SEP = ',');
void printOrigProbInfo(const OsiSolverInterface* const solver, FILE* logfile,
    const char SEP = ',');
void printPostCutProbInfo(const OsiSolverInterface* const solver,
    const SummaryCutInfo& cutInfoGMICs, const SummaryCutInfo& cutInfoVPCs,
    FILE* logfile, const char SEP = ',');
void printDisjInfo(const SummaryDisjunctionInfo& disjInfo, FILE* logfile,
    const char SEP = ',');
void printCutInfo(const SummaryCutInfo& cutInfoGMICs,
    const SummaryCutInfo& cutInfo, FILE* logfile, const char SEP = ',');

void analyzeStrength(const VPCParametersNamespace::VPCParameters& params, const OsiSolverInterface* solver,
    SummaryCutInfo& cutInfoGMICs, SummaryCutInfo& cutInfoVPCs, 
    const OsiCuts* const gmics, const OsiCuts* const vpcs,
    const SummaryBoundInfo& boundInfo, std::string& output);
void analyzeBB(const VPCParametersNamespace::VPCParameters& params, SummaryBBInfo& info_nocuts,
    SummaryBBInfo& info_mycuts, SummaryBBInfo& info_allcuts, std::string& output);
double getNumGomoryRounds(const VPCParametersNamespace::VPCParameters& params,
    const OsiSolverInterface* const origSolver,
    const OsiSolverInterface* const postCutSolver);

void updateDisjInfo(SummaryDisjunctionInfo& disjInfo, const int num_disj, const CglVPC& gen);
void updateCutInfo(SummaryCutInfo& cutInfo, const CglVPC& gen);
void setCutInfo(SummaryCutInfo& cutInfo, const int num_rounds, const SummaryCutInfo* const oldCutInfos);

