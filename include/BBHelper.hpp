//============================================================================
// Name        : BBHelper.hpp
// Author      : akazachk
// Version     : 2018.11.19
// Description : Helper functions for branch-and-bound
//
// This code is licensed under the terms of the Eclipse Public License (EPL).
//============================================================================
#pragma once

// Project files
#include "VPCParameters.hpp"

class PartialBBDisjunction;

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

struct BBInfo {
  double obj;
  double bound;
  long iters;
  long nodes;
  long root_passes;
  double first_cut_pass;
  double last_cut_pass;
  double root_time;
  double last_sol_time;
  double time;
}; /* BBInfo */

enum BBInfoEnum {
  OBJ_BB_INFO_IND,
  BOUND_BB_INFO_IND,
  ITERS_BB_INFO_IND,
  NODES_BB_INFO_IND,
  ROOT_PASSES_BB_INFO_IND,
  FIRST_CUT_PASS_BB_INFO_IND,
  LAST_CUT_PASS_BB_INFO_IND,
  ROOT_TIME_BB_INFO_IND,
  LAST_SOL_TIME_BB_INFO_IND,
  TIME_BB_INFO_IND,
  NUM_BB_INFO
};

const std::vector<std::string> BB_INFO_CONTENTS = {
    "OBJ", "BOUND", "ITERS", "NODES", "ROOT_PASSES", "FIRST_CUT_PASS", "LAST_CUT_PASS", "ROOT_TIME", "LAST_SOL_TIME", "TIME"
};

inline void initializeBBInfo(BBInfo& info, double obj = 0.) {
  info.obj = obj;
  info.bound = obj;
  info.iters = 0;
  info.nodes = 0;
  info.root_passes = 0;
  info.first_cut_pass = 0.;
  info.last_cut_pass = 0.;
  info.root_time = 0.;
  info.last_sol_time = 0.;
  info.time = 0.;
}
void updateBestBBInfo(BBInfo& min_info, const BBInfo& curr_info, const bool first);
void averageBBInfo(BBInfo& avg_info, const std::vector<BBInfo>& info);
void printBBInfo(const BBInfo& info, FILE* myfile, const bool print_blanks =
    false, const char SEP = ',');
void printBBInfo(const BBInfo& info_mycuts, const BBInfo& info_allcuts,
    FILE* myfile, const bool print_blanks = false, const char SEP = ',');
void createStringFromBBInfoVec(const std::vector<BBInfo>& vec_info,
    std::vector<std::string>& vec_str);

// COIN-OR
#ifdef USE_CBC
#include "CbcModel.hpp"

/**
 * Sets message handler and special options when using solver as part of B&B
 * (in which we want to run full strong branching and enable the fixing of variables)
 */
void setupClpForCbc(OsiClpSolverInterface* const solver,
    const int hot_start_iter_limit = std::numeric_limits<int>::max()); 

/**
 * Set parameters for Cbc used for VPCs, as well as the custom branching decision
 */
void setCbcParametersForPartialBB(
    const VPCParameters& param,
    CbcModel* const cbc_model,
    CbcEventHandler* eventHandler = NULL, const int numStrong = 5,
    const int numBeforeTrusted = 10, const double max_time =
        std::numeric_limits<double>::max());

void generatePartialBBTree(PartialBBDisjunction* const owner, CbcModel* cbc_model,
    const OsiSolverInterface* const solver, const int max_nodes,
    const int num_strong, const int num_before_trusted);
#endif // USE_CBC
