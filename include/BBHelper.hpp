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

enum BB_Strategy_Options {
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
};

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

void generatePartialBBTree(const VPCParameters& param, CbcModel* cbc_model,
    const OsiSolverInterface* const solver, const int max_nodes,
    const int num_strong, const int num_before_trusted);
#endif // USE_CBC
