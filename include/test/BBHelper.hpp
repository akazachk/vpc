/**
 * @file BBHelper.hpp
 * @author A. M. Kazachkov
 * @date 2018-11-19
 * @brief Helper functions for branch-and-bound
 */
#pragma once

#include <vector>
#include <string>

class PartialBBDisjunction;
class OsiSolverInterface;
class TimeStats;
class OsiCuts;

namespace VPCParametersNamespace {
  struct VPCParameters;
}

// Defined here
struct BBInfo;
struct SummaryBBInfo;

enum class BB_Strategy_Options {
  off = 0,
  cbc = 2,
  cplex = 4,
  gurobi = 8,
  user_cuts = 16,
  all_cuts_off = 32,
  all_cuts_on = 64,
  gmics_off = 128,
  gmics_on = 256,
  presolve_off = 512,
  presolve_on = 1024,
  heuristics_off = 2048,
  heuristics_on = 4096,
  use_best_bound = 8192,
  strong_branching_on = 16384
}; /* BB_Strategy_Options */

inline bool use_bb_option(const int strategy, const BB_Strategy_Options option) {
  return strategy & static_cast<int>(option);
}
inline int enable_bb_option(const int strategy, const BB_Strategy_Options option) {
  return strategy | static_cast<int>(option);
}
inline int disable_bb_option(const int strategy, const BB_Strategy_Options option) {
  return strategy & ~static_cast<int>(option);
}

struct BBInfo {
  double obj = 0.; // objective value of best IP-feasible solution
  double bound = 0.; // best dual bound found (best objective of any leaf node)
  long iters = 0; // # iters to solve the instance
  long nodes = 0; // # nodes to solve the instance
  long root_passes = 0; // # passes of cuts at the root node
  double first_cut_pass = 0.; // bound after one round of cuts at the root
  double last_cut_pass = 0.; // bound after last round of cuts at the root
  long root_iters = 0; // # iters spent at root node
  double root_time = 0.; // time spent at the root node
  double last_sol_time = 0.; // time that best IP-feasible solution was found
  double time = 0.; // total time to solve the instance
}; /* BBInfo */
enum BBInfoEnum {
  OBJ_BB_INFO_IND,
  BOUND_BB_INFO_IND,
  ITERS_BB_INFO_IND,
  NODES_BB_INFO_IND,
  ROOT_PASSES_BB_INFO_IND,
  FIRST_CUT_PASS_BB_INFO_IND,
  LAST_CUT_PASS_BB_INFO_IND,
  ROOT_ITERS_BB_INFO_IND,
  ROOT_TIME_BB_INFO_IND,
  LAST_SOL_TIME_BB_INFO_IND,
  TIME_BB_INFO_IND,
  NUM_BB_INFO
};
const std::vector<std::string> BB_INFO_CONTENTS = {
    "OBJ", "BOUND", "ITERS", "NODES", "ROOT_PASSES", "FIRST_CUT_PASS", "LAST_CUT_PASS", "ROOT_ITERS", "ROOT_TIME", "LAST_SOL_TIME", "TIME"
};

struct SummaryBBInfo {
  int num_cuts = 0;
  BBInfo first_bb_info, best_bb_info, avg_bb_info;
  std::vector<BBInfo> vec_bb_info;
};

void runBBTests(const VPCParametersNamespace::VPCParameters& base_params, SummaryBBInfo* const info_nocuts,
    SummaryBBInfo* const info_mycuts, SummaryBBInfo* const info_allcuts,
    const std::string fullfilename, const OsiSolverInterface* const solver,
    const double best_bound, const OsiCuts* vpcs, const OsiCuts* const gmics = NULL,
    const std::vector<double>& ip_solution = {});

inline void initializeBBInfo(BBInfo& info, double obj = 0.) {
  info.obj = obj;
  info.bound = obj;
  info.iters = 0;
  info.nodes = 0;
  info.root_passes = 0;
  info.first_cut_pass = 0.;
  info.last_cut_pass = 0.;
  info.root_iters = 0;
  info.root_time = 0.;
  info.last_sol_time = 0.;
  info.time = 0.;
}
void updateBestBBInfo(BBInfo& min_info, const BBInfo& curr_info, const bool first);
void averageBBInfo(BBInfo& avg_info, const std::vector<BBInfo>& info);

void createStringFromBBInfoVec(const std::vector<BBInfo>& vec_info,
    std::vector<std::string>& vec_str);

/**
 * Creates temporary file (in /tmp) so that it can be read by a different solver
 * It does not delete the file
 */
void createTmpFileCopy(const VPCParametersNamespace::VPCParameters& params,
    const OsiSolverInterface* const solver, std::string& f_name);

// COIN-OR
#ifdef USE_CBC
void doBranchAndBoundNoCuts(const VPCParametersNamespace::VPCParameters& params, const OsiSolverInterface* const solver, BBInfo& info);
void doBranchAndBoundYesCuts(const VPCParametersNamespace::VPCParameters& params, const OsiSolverInterface* const solver,
    BBInfo& info, const OsiCuts& structCuts, const bool doCutSelection,
    const int numCutsToAddPerRound, const int maxRounds,
    const std::string logstring);
#endif /* USE_CBC */
