/**
 * @file analysis.hpp
 * @author A. M. Kazachkov
 * @date 2019-03-04
 */
#pragma once

#include <string>
#include <vector>
#include <limits>

class OsiSolverInterface;
class OsiCuts;

#include "CglVPC.hpp" // CutType, ObjectiveType
namespace VPCParametersNamespace {
  struct VPCParameters;
}

struct SummaryBBInfo; // BBHelper.hpp

// Defined here
struct SummaryBoundInfo; // analysis.hpp
struct SummaryDisjunctionInfo; // analysis.hpp
struct SummaryCutInfo; // analysis.hpp

/// @brief Information about objective value at various points in the solution process
/// @details Gives objective for the LP and IP, and after adding GMICs, L&PCs, VPCs, and combinations of these cuts
/// and also keeps number of GMICs, L&PCs, and VPCs applied
struct SummaryBoundInfo {
  double lp_obj = std::numeric_limits<double>::max();
  double best_disj_obj = std::numeric_limits<double>::lowest();
  double worst_disj_obj = std::numeric_limits<double>::lowest();
  double root_obj = std::numeric_limits<double>::lowest();
  double ip_obj = std::numeric_limits<double>::max();
  double gmic_obj = std::numeric_limits<double>::max();
  double lpc_obj = std::numeric_limits<double>::max();
  double vpc_obj = std::numeric_limits<double>::max();
  double gmic_vpc_obj = std::numeric_limits<double>::max();
  double all_cuts_obj = std::numeric_limits<double>::max();
  int num_root_bounds_changed = 0, num_gmic = 0, num_lpc = 0, num_vpc = 0;
}; /* SummaryBoundInfo */

/// @brief Summary statistics for the disjunction generated
struct SummaryDisjunctionInfo {
  int num_disj = 0;
  int num_integer_sol = 0;
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

/// @brief Track statistics for a particular family of cuts
struct SummaryCutInfo {
  int num_cuts = 0;
  int num_active_gmic = 0, num_active_lpc = 0, num_active_vpc, num_active_all = 0;
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

/// @brief Write header for output to log in \p params.logfile
void printHeader(const VPCParametersNamespace::VPCParameters& params,
    const std::vector<std::string>& time_name,
    const char SEP = ',');
/// @brief Write to log objective values and gaps if IP value is known
void printBoundAndGapInfo(const SummaryBoundInfo& boundInfo, FILE* logfile,
    const char SEP = ',');
/// @brief Write to log summary of branch-and-bound time and number nodes
void printSummaryBBInfo(const std::vector<SummaryBBInfo>& info, FILE* logfile,
    const bool print_blanks = false, const char SEP = ',');
/// @brief Write to log full branch-and-bound information
void printFullBBInfo(const std::vector<SummaryBBInfo>& info, FILE* logfile,
    const bool print_blanks = false, const char SEP = ',');
/// @brief Write to log statistics about number of rows, cols, etc. for original problem
void printOrigProbInfo(const OsiSolverInterface* const solver, FILE* logfile,
    const char SEP = ',');
/// @brief Write to log statistics about number of rows, cols, etc. about problem after adding cuts, and information about active cuts
void printPostCutProbInfo(const OsiSolverInterface* const solver,
    const SummaryCutInfo& cutInfoGMICs, const SummaryCutInfo& cutInfoVPCs,
    FILE* logfile, const char SEP = ',');
/// @brief Write to log statistics about the disjunction stored in #SummaryDisjunctionInfo
void printDisjInfo(const SummaryDisjunctionInfo& disjInfo, FILE* logfile,
    const char SEP = ',');
/// @brief Write to log information about cuts generated and applied stored in #SummaryCutInfo
void printCutInfo(const SummaryCutInfo& cutInfoGMICs,
    const SummaryCutInfo& cutInfo, FILE* logfile, const char SEP = ',');

/// @brief Check cut density and update min/max support in \p cutInfo
int getCutSupport(
    const OsiRowCut* const cut,
    const double EPS = 1e-14);

/// @brief Check cut activity in solver and report cut density
bool checkCutActivity(
  const OsiSolverInterface* const solver,
  const OsiRowCut* const cut);

/// @brief Compute gap closed and active cuts
void analyzeStrength(const VPCParametersNamespace::VPCParameters& params, 
    const OsiSolverInterface* const solver_gmic,
    const OsiSolverInterface* const solver_vpc,
    const OsiSolverInterface* const solver_all,
    SummaryCutInfo& cutInfoGMICs, SummaryCutInfo& cutInfoVPCs, 
    const OsiCuts* const gmics, const OsiCuts* const vpcs,
    const SummaryBoundInfo& boundInfo, std::string& output);
/// @brief Prepare summary of #SummaryBBInfo
void analyzeBB(const VPCParametersNamespace::VPCParameters& params, SummaryBBInfo& info_nocuts,
    SummaryBBInfo& info_mycuts, SummaryBBInfo& info_allcuts, std::string& output);
/// @brief Get number rounds of SICs needed to meet bound from VPCs+SICs
double getNumGomoryRounds(const VPCParametersNamespace::VPCParameters& params,
    const OsiSolverInterface* const origSolver,
    const OsiSolverInterface* const postCutSolver);

/// @brief Update average number of terms, density, rows, cols, points, rays, and partial tree information if applicable
void updateDisjInfo(SummaryDisjunctionInfo& disjInfo, const int num_disj, const CglVPC& gen);
/// @brief Add to cut information after a round, such as number of cuts, objectives, failures, etc.
void updateGMICInfo(SummaryCutInfo& cutInfo, const OsiCuts* cuts, const double EPS = 1e-14);
/// @brief Add to cut information after a round, such as number of cuts, objectives, failures, etc.
void updateCutInfo(SummaryCutInfo& cutInfo, const CglVPC& gen, const OsiCuts* cuts = NULL, const double EPS = 1e-14);
/// @brief Use this to merge cut info from multiple rounds
void setCutInfo(SummaryCutInfo& cutInfo, const int num_rounds, const SummaryCutInfo* const oldCutInfos);

