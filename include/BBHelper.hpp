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
