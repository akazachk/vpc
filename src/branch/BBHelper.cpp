//============================================================================
// Name        : BBHelper.cpp
// Author      : A. M. Kazachkov
// Version     : 2018-Dec-24
// Description : Helper functions for branch-and-bound
//
// This code is licensed under the terms of the Eclipse Public License (EPL).
//============================================================================

#include <cstdio> // for tmpnam

// Project files 
#include "BBHelper.hpp"
#include "VPCEventHandler.hpp"
#include "utility.hpp"

// COIN-OR
#include <CoinTime.hpp>

// Cbc
#include <CbcSolver.hpp>
#include <CbcTree.hpp>

// Variable selection
#include "CbcBranchDefaultDecision.hpp"
#include "CbcBranchStrongDecision.hpp"
#include "OsiChooseStrongCustom.hpp"
#include "CbcBranchDynamic.hpp"

// Node selection
#include "CbcCompareDefault.hpp"
#include "CbcCompareBFS.hpp"
#include "CbcCompareDepth.hpp"
#include "CbcCompareEstimate.hpp"
#include "CbcCompareObjective.hpp"

// General strategy
#include "CbcStrategy.hpp"

/**
 * Sets message handler and special options when using solver as part of B&B
 * (in which we want to run full strong branching and enable the fixing of variables)
 */
void setupClpForCbc(OsiClpSolverInterface* const solver,
    const int hot_start_iter_limit) {
  setLPSolverParameters(solver);
  solver->setHintParam(OsiDoPresolveInInitial, false);
  solver->setHintParam(OsiDoPresolveInResolve, false);
  solver->setIntParam(OsiMaxNumIterationHotStart, hot_start_iter_limit);
  solver->setSpecialOptions(16); // use standard strong branching rather than clp's
  // Do not switch from dual to primal, or something to this effect;
  // This allows infeasible branches to be fixed during strong branching
  solver->getModelPtr()->setMoreSpecialOptions(solver->getModelPtr()->moreSpecialOptions()+256);
} /* setupClpForCbc */

/**
 * Set parameters for Cbc used for VPCs, as well as the custom branching decision
 */
void setCbcParametersForPartialBB(
		const VPCParameters& params,
		CbcModel* const cbc_model,
    CbcEventHandler* eventHandler,
    const int numStrong,
    const int numBeforeTrusted,
    const double max_time) {
  setIPSolverParameters(cbc_model);
  cbc_model->solver()->setIntParam(OsiMaxNumIterationHotStart, 100);

  // What is the partial strategy?
	const int strategy = std::abs(params.get(intParam::PARTIAL_BB_STRATEGY));
  const int sign = (params.get(intParam::PARTIAL_BB_STRATEGY) < 0) ? -1 : 1;
  const int compare_strategy = strategy % 10; // ones digit
  const int branch_strategy = (strategy % 100 - compare_strategy) / 10; // tens digit
  const int choose_strategy = sign * (strategy % 1000 - compare_strategy - 10 * branch_strategy) / 100; // hundreds digit

  // Branching decision (tens digit)
  // Given a branching variable, which direction to choose?
  // (The choice of branching variable is through OsiChooseVariable)
  // 0: default, 1: dynamic, 2: strong, 3: none
  CbcBranchDecision* branch;
  if (branch_strategy == 1) {
    branch = new CbcBranchDynamicDecision();
  } else if (branch_strategy == 2) {
    branch = new CbcBranchStrongDecision();
  } else if (branch_strategy == 3) {
    branch = NULL;
  } else {
    branch = new CbcBranchDefaultDecision();
  }

  // Set comparison for nodes (ones digit)
  // Given a tree, which node to pick next?
  // 0: default: 1: bfs, 2: depth, 3: estimate, 4: objective, 5: objective_reverse
  CbcCompareBase* compare;
  if (compare_strategy == 1) {
    compare = new CbcCompareBFS();
  } else if (compare_strategy == 2) {
    compare = new CbcCompareDepth();
  } else if (compare_strategy == 3) {
    compare = new CbcCompareEstimate();
  } else if (compare_strategy == 4) {
    compare = new CbcCompareObjective();
  } else {
    compare = new CbcCompareDefault();
  }

  cbc_model->setTypePresolve(0);
  cbc_model->setMaximumSeconds(max_time);
  cbc_model->setMaximumCutPassesAtRoot(0);
  cbc_model->setMaximumCutPasses(0);
  cbc_model->setWhenCuts(0);
  if (numStrong >= 0) {
    // Maximum number of strong branching candidates to consider each time
    cbc_model->setNumberStrong(numStrong);
  }
  if (numBeforeTrusted >= 0) {
    // # before switching to pseudocosts, I think; 0 disables dynamic strong branching, doesn't work well
    cbc_model->setNumberBeforeTrust(numBeforeTrusted);
  }

  if (branch) {
    OsiChooseStrongCustom choose;
		if (numStrong >= 0) {
      choose.setNumberStrong(numStrong);
    }
    if (numBeforeTrusted >= 0) {
      choose.setNumberBeforeTrusted(numBeforeTrusted);
    }
    choose.setMethod(choose_strategy);
    branch->setChooseMethod(choose);
    // From CbcModel::convertToDynamic, we see that the branching decision may be ignored without a choose method
    if ((branch->whichMethod()&1) == 0 && !branch->chooseMethod()) {
      OsiChooseStrong choose;
      choose.setNumberStrong(5);
      choose.setNumberBeforeTrusted(0);
      branch->setChooseMethod(choose);
    }
    cbc_model->setBranchingMethod(*branch);
  } /* check that branch is not NULL */

  if (eventHandler) {
    cbc_model->passInEventHandler(eventHandler);
  }

  cbc_model->setNodeComparison(compare);

  if (branch) {
    delete branch;
  }
  if (compare) {
    delete compare;
  }
} /* setCbcParametersForPartialBB */

/************************************************************/
/**
 * Generate a partial branch-and-bound tree with at most max_leaf_nodes leaf nodes 
 */
void generatePartialBBTree(const VPCParameters& params, CbcModel* cbc_model,
    const OsiSolverInterface* const solver, const int max_leaf_nodes,
    const int num_strong, const int num_before_trusted) {
  //  double partial_timelimit = 100 * max_leaf_nodes * solver->getNumCols()
  //      * GlobalVariables::timeStats.get_time(GlobalConstants::INIT_SOLVE_TIME);
  const double partial_timelimit = params.get(PARTIAL_BB_TIMELIMIT); // will be checked manually by the eventHandler

  // Set up options
  VPCEventHandler* eventHandler = new VPCEventHandler(max_leaf_nodes, partial_timelimit, &params);
  eventHandler->setOriginalSolver(solver);

  // This sets branching decision, event handling, etc.
  setCbcParametersForPartialBB(params, cbc_model, eventHandler, num_strong,
      num_before_trusted, std::numeric_limits<double>::max());

#ifdef TRACE
  cbc_model->branchAndBound(3);
#else
  cbc_model->branchAndBound(0);
#endif

  // Free
  // When eventHandler is passed, Cbc currently clones it and does not delete the original
  if (eventHandler) { // in case this behavior gets changed in the future
    delete eventHandler;
  }
} /* generatePartialBBTree */
