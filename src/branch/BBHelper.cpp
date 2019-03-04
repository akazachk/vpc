//============================================================================
// Name        : BBHelper.cpp
// Author      : A. M. Kazachkov
// Version     : 2018-Dec-24
// Description : Helper functions for branch-and-bound
//
// This code is licensed under the terms of the Eclipse Public License (EPL).
//============================================================================

#include <cstdio> // for tmpnam
#include <algorithm>  // std::random_shuffle

// Project files 
#include "BBHelper.hpp"
#include "CutHelper.hpp" // applyCutsCustom
#include "PartialBBDisjunction.hpp"
#include "SolverHelper.hpp"
#include "VPCEventHandler.hpp"
#include "VPCParameters.hpp"
#include "TimeStats.hpp"
#include "utility.hpp"

// COIN-OR
#include <CoinTime.hpp>
#include <OsiCuts.hpp>

#ifdef USE_CBC
// Cbc
#include <CbcSolver.hpp>
#include <CbcTree.hpp>

// Variable selection
#include <CbcBranchDefaultDecision.hpp>
#include <CbcBranchStrongDecision.hpp>
#include <OsiChooseStrongCustom.hpp>
#include <CbcBranchDynamic.hpp>

// Node selection
#include <CbcCompareDefault.hpp>
#include <CbcCompareBFS.hpp>
#include <CbcCompareDepth.hpp>
#include <CbcCompareEstimate.hpp>
#include <CbcCompareObjective.hpp>

// General strategy
#include <CbcStrategy.hpp>
#endif /* USE_CBC */

#ifdef USE_GUROBI
#include "GurobiHelper.hpp"
#endif

void runBBTests(const VPCParameters& params, SummaryBBInfo& info_nocuts,
    SummaryBBInfo& info_mycuts, SummaryBBInfo& info_allcuts,
    const std::string fullfilename, OsiSolverInterface* const solver,
    const double best_bound, const OsiCuts& vpcs, const OsiCuts* const gmics) {
  const int num_bb_runs = std::abs(params.get(intParam::BB_RUNS))
      * ((params.get(intParam::DISJ_TERMS) == 0) || vpcs.sizeCuts() > 0);
  if (num_bb_runs == 0)
    return;

  /***********************************************************************************
   * Perform branch-and-bound
   ***********************************************************************************/
  const int numCutsToAddPerRound =
      (gmics != NULL && gmics->sizeCuts() > 0) ? gmics->sizeCuts() : vpcs.sizeCuts();
  const bool should_test_bb_with_gmics = (params.get(BB_RUNS) > 0) && (gmics != NULL) && (gmics->sizeCuts() > 0);
  const bool should_permute_rows_and_cols = (num_bb_runs >= 2)
      && !use_bb_option(params.get(BB_STRATEGY), BB_Strategy_Options::gurobi)
      && !use_bb_option(params.get(BB_STRATEGY), BB_Strategy_Options::cplex);

#ifdef TRACE
  printf("\n## Performing branch-and-bound tests. ##\n");
#endif

  info_nocuts.vec_bb_info.resize(0);
  info_nocuts.vec_bb_info.resize(num_bb_runs);
  info_mycuts.vec_bb_info.resize(0);
  info_mycuts.vec_bb_info.resize(num_bb_runs);
  info_allcuts.vec_bb_info.resize(0);
  info_allcuts.vec_bb_info.resize(num_bb_runs);
  std::vector<int> row_permutation, col_permutation;
  if (should_permute_rows_and_cols) {
    row_permutation.resize(solver->getNumRows());
    col_permutation.resize(solver->getNumRows());
    for (int row = 0; row < solver->getNumRows(); row++) {
      row_permutation[row] = row;
    }
    for (int col = 0; col < solver->getNumCols(); col++) {
      col_permutation[col] = col;
    }
  }

  OsiCuts GMICsAndVPCs = vpcs;
  if (should_test_bb_with_gmics) {
    if (gmics)
      GMICsAndVPCs.insert(*gmics);
  }

  OsiSolverInterface* runSolver = NULL;
  int initial_random_seed = params.get(intConst::RANDOM_SEED);
  int random_seed = params.get(intConst::RANDOM_SEED);
  for (int run_ind = 0; run_ind < num_bb_runs; run_ind++) {
    // For Cbc, for every run after the first, we will randomize the rows and columns of the input
    if (should_permute_rows_and_cols) {
      std::srand(run_ind + 2);  // sets the seed for pseudo-randomness in C++
      std::vector<int> this_row_permutation(row_permutation);
      std::vector<int> this_col_permutation(col_permutation);
      std::random_shuffle(this_row_permutation.begin(), this_row_permutation.end());
      std::random_shuffle(this_col_permutation.begin(), this_col_permutation.end());

      // Create the permuted problem
      std::vector<double> rowLB(solver->getNumRows()), rowUB(solver->getNumRows());
      std::vector<double> colLB(solver->getNumCols()), colUB(solver->getNumCols());
      std::vector<double> obj(solver->getNumCols());
      const CoinPackedMatrix* mat = solver->getMatrixByRow();

      // Set the rows
      CoinPackedMatrix row_mat;
      row_mat.reverseOrdering();  // make it row-ordered
      row_mat.setMinorDim(solver->getNumCols());
      for (int row = 0; row < solver->getNumRows(); row++) {
        const int curr_row = this_row_permutation[row];
        row_mat.appendRow(mat->getVector(curr_row));
        rowLB[row] = solver->getRowLower()[curr_row];
        rowUB[row] = solver->getRowUpper()[curr_row];
        printf("New row: %d. Orig row: %d. Row upper: %s. Row lower: %s.\n",
            row, curr_row, stringValue(rowLB[row]).c_str(),
            stringValue(rowUB[row]).c_str());
      }
      row_mat.reverseOrdering(); // make it col-ordered for use in the next part

      // Set the cols
      CoinPackedMatrix col_mat;
      col_mat.setMinorDim(solver->getNumRows());
      for (int col = 0; col < solver->getNumCols(); col++) {
        const int curr_col = this_col_permutation[col];
        col_mat.appendCol(row_mat.getVector(curr_col));
        colLB[col] = solver->getColLower()[curr_col];
        colUB[col] = solver->getColUpper()[curr_col];
        obj[col] = solver->getObjCoefficients()[curr_col];
      }

      runSolver = new SolverInterface;
      runSolver->setObjSense(solver->getObjSense());
      double objOffset = 0.;
      solver->getDblParam(OsiDblParam::OsiObjOffset, objOffset);
      runSolver->setDblParam(OsiDblParam::OsiObjOffset, objOffset);

      runSolver->loadProblem(col_mat, colLB.data(), colUB.data(), obj.data(),
          rowLB.data(), rowUB.data());

      // Set integers
      for (int col = 0; col < solver->getNumCols(); col++) {
        const int curr_col = this_col_permutation[col];
        if (solver->isInteger(curr_col)) {
          runSolver->setInteger(curr_col);
        }
      }
    } // test whether num_bb_runs >= 2 and permute row/col order
    else {
      runSolver = solver;
    }

    // Change the random seed per run
    random_seed = initial_random_seed * (run_ind + 1);

    // Do branch and bound
    if (use_bb_option(params.get(BB_STRATEGY), BB_Strategy_Options::gurobi)) {
#ifdef USE_GUROBI
      if (params.get(TEMP) == static_cast<int>(TempOptions::CHECK_CUTS_AGAINST_BB_OPT)) {
        // Get the original solution
        BBInfo tmp_bb_info;
        std::vector<double> solution;
        doBranchAndBoundWithGurobi(params, params.get(BB_STRATEGY),
            fullfilename.c_str(), tmp_bb_info, best_bound, &solution);

        // Check cuts
        for (int cut_ind = 0; cut_ind < vpcs.sizeCuts(); cut_ind++) {
          OsiRowCut currCut = vpcs.rowCut(cut_ind);
          const double rhs = currCut.rhs();
          const int num_el = currCut.row().getNumElements();
          const int* ind = currCut.row().getIndices();
          const double* el = currCut.row().getElements();
          const double activity = dotProduct(num_el, ind, el, solution.data());

          if (lessThanVal(activity, rhs)) {
            warning_msg(warnstring, "Cut %d removes optimal solution. Activity: %.10f. Rhs: %.10f.\n", cut_ind, activity, rhs);
          }
        }
      } // checking cuts for violating the IP opt
//      doBranchAndBoundWithGurobi(params, disable_bb_option(params.get(BB_STRATEGY), BB_Strategy_Options::user_cuts),
//          fullfilename.c_str(), info_nocuts.vec_bb_info[run_ind], best_bound);
      doBranchAndBoundWithUserCutsGurobi(params, params.get(BB_STRATEGY),
          fullfilename.c_str(), NULL, info_nocuts.vec_bb_info[run_ind],
          best_bound);
      doBranchAndBoundWithUserCutsGurobi(params, params.get(BB_STRATEGY),
          fullfilename.c_str(), &vpcs, info_mycuts.vec_bb_info[run_ind],
          best_bound);
      if (should_test_bb_with_gmics) {
        doBranchAndBoundWithUserCutsGurobi(params, params.get(BB_STRATEGY),
            fullfilename.c_str(), &GMICsAndVPCs,
            info_allcuts.vec_bb_info[run_ind], best_bound);
      }
#endif // USE_GUROBI
    } else if (use_bb_option(params.get(BB_STRATEGY), BB_Strategy_Options::cplex)) {
#ifdef USE_CPLEX
      doBranchAndBoundWithUserCutsCplexCallable(
          GlobalVariables::in_f_name.c_str(), &structVPC, vec_bb_info_mycuts[run_ind]);
      doBranchAndBoundWithUserCutsCplexCallable(
          GlobalVariables::in_f_name.c_str(), &SICsAndVPCs,
          vec_bb_info_allcuts[run_ind]);
#endif // USE_CPLEX
    } else if (use_bb_option(params.get(BB_STRATEGY), BB_Strategy_Options::cbc)) {
      doBranchAndBoundNoCuts(params, runSolver, info_nocuts.vec_bb_info[run_ind]);
      doBranchAndBoundYesCuts(params, runSolver, info_mycuts.vec_bb_info[run_ind],
          vpcs, false, numCutsToAddPerRound, 1, "\nBB with VPCs.\n");
      if (should_test_bb_with_gmics) {
        doBranchAndBoundYesCuts(params, runSolver, info_allcuts.vec_bb_info[run_ind],
            GMICsAndVPCs, false, numCutsToAddPerRound, 1, "\nBB with VPC+Gomorys.\n");
      }
    }

    updateBestBBInfo(info_nocuts.best_bb_info, info_nocuts.vec_bb_info[run_ind], (run_ind == 0));
    updateBestBBInfo(info_mycuts.best_bb_info, info_mycuts.vec_bb_info[run_ind], (run_ind == 0));
    updateBestBBInfo(info_allcuts.best_bb_info, info_allcuts.vec_bb_info[run_ind], (run_ind == 0));

    // Free memory if necessary
    if (should_permute_rows_and_cols && runSolver && (runSolver != solver)) {
      delete runSolver;
    }
  } /* end iterating over runs */

  averageBBInfo(info_nocuts.avg_bb_info, info_nocuts.vec_bb_info);
  averageBBInfo(info_mycuts.avg_bb_info, info_mycuts.vec_bb_info);
  averageBBInfo(info_allcuts.avg_bb_info, info_allcuts.vec_bb_info);
  writeBBInforToLog(info_mycuts, info_allcuts, params.logfile, 2);
} /* runBBTests */

/** Methods related to BBInfo */
void updateBestBBInfo(BBInfo& best_info, const BBInfo& curr_info, const bool first) {
  best_info.obj = first ? curr_info.obj : CoinMin(best_info.obj, curr_info.obj);
  best_info.bound = first ? curr_info.bound : CoinMax(best_info.bound, curr_info.bound);
  best_info.iters = first ? curr_info.iters : CoinMin(best_info.iters, curr_info.iters);
  best_info.nodes = first ? curr_info.nodes : CoinMin(best_info.nodes, curr_info.nodes);
  best_info.root_passes = first ? curr_info.root_passes : CoinMin(best_info.root_passes, curr_info.root_passes);
  best_info.first_cut_pass = first ? curr_info.first_cut_pass : CoinMax(best_info.first_cut_pass, curr_info.first_cut_pass);
  best_info.last_cut_pass = first ? curr_info.last_cut_pass : CoinMax(best_info.last_cut_pass, curr_info.last_cut_pass);
  best_info.root_time = first ? curr_info.root_time : CoinMin(best_info.root_time, curr_info.root_time);
  best_info.last_sol_time = first ? curr_info.last_sol_time : CoinMin(best_info.last_sol_time, curr_info.last_sol_time);
  best_info.time = first ? curr_info.time : CoinMin(best_info.time, curr_info.time);
} /* updateMinBBInfo */

void averageBBInfo(BBInfo& avg_info, const std::vector<BBInfo>& info) {
  for (BBInfo curr_info : info) {
    avg_info.obj += curr_info.obj;
    avg_info.bound += curr_info.bound;
    avg_info.iters += curr_info.iters;
    avg_info.nodes += curr_info.nodes;
    avg_info.root_passes += curr_info.root_passes;
    avg_info.first_cut_pass += curr_info.first_cut_pass;
    avg_info.last_cut_pass += curr_info.last_cut_pass;
    avg_info.root_time += curr_info.root_time;
    avg_info.last_sol_time += curr_info.last_sol_time;
    avg_info.time += curr_info.time;
  }
  const int num_bb_runs = info.size();
  avg_info.obj /= num_bb_runs;
  avg_info.bound /= num_bb_runs;
  avg_info.iters /= num_bb_runs;
  avg_info.nodes /= num_bb_runs;
  avg_info.root_passes /= num_bb_runs;
  avg_info.first_cut_pass /= num_bb_runs;
  avg_info.last_cut_pass /= num_bb_runs;
  avg_info.root_time /= num_bb_runs;
  avg_info.last_sol_time /= num_bb_runs;
  avg_info.time /= num_bb_runs;
} /* averageBBInfo */

void printBBInfo(const BBInfo& info, FILE* myfile, const bool print_blanks, const char SEP) {
  if (!print_blanks) {
    fprintf(myfile, "%s%c", stringValue(info.obj, "%.20f").c_str(), SEP);
    fprintf(myfile, "%s%c", stringValue(info.bound, "%.20f").c_str(), SEP);
    fprintf(myfile, "%ld%c", info.iters, SEP);
    fprintf(myfile, "%ld%c", info.nodes, SEP);
    fprintf(myfile, "%ld%c", info.root_passes, SEP);
    fprintf(myfile, "%.20f%c", info.first_cut_pass, SEP);
    fprintf(myfile, "%.20f%c", info.last_cut_pass, SEP);
    fprintf(myfile, "%2.3f%c", info.root_time, SEP);
    fprintf(myfile, "%2.3f%c", info.last_sol_time, SEP);
    fprintf(myfile, "%2.3f%c", info.time, SEP);
  } else {
    for (int i = 0; i < (int) BB_INFO_CONTENTS.size(); i++) {
      fprintf(myfile, "%c", SEP);
    }
  }
} /* printBBInfo */

void printBBInfo(const BBInfo& info_mycuts, const BBInfo& info_allcuts,
    FILE* myfile, const bool print_blanks, const char SEP) {
  if (!print_blanks) {
    fprintf(myfile, "%s%c", stringValue(info_mycuts.obj, "%.20f").c_str(), SEP);
    fprintf(myfile, "%s%c", stringValue(info_allcuts.obj, "%.20f").c_str(), SEP);
    fprintf(myfile, "%s%c", stringValue(info_mycuts.bound, "%.20f").c_str(), SEP);
    fprintf(myfile, "%s%c", stringValue(info_allcuts.bound, "%.20f").c_str(), SEP);
    fprintf(myfile, "%ld%c", info_mycuts.iters, SEP);
    fprintf(myfile, "%ld%c", info_allcuts.iters, SEP);
    fprintf(myfile, "%ld%c", info_mycuts.nodes, SEP);
    fprintf(myfile, "%ld%c", info_allcuts.nodes, SEP);
    fprintf(myfile, "%ld%c", info_mycuts.root_passes, SEP);
    fprintf(myfile, "%ld%c", info_allcuts.root_passes, SEP);
    fprintf(myfile, "%.20f%c", info_mycuts.first_cut_pass, SEP);
    fprintf(myfile, "%.20f%c", info_allcuts.first_cut_pass, SEP);
    fprintf(myfile, "%.20f%c", info_mycuts.last_cut_pass, SEP);
    fprintf(myfile, "%.20f%c", info_allcuts.last_cut_pass, SEP);
    fprintf(myfile, "%2.3f%c", info_mycuts.root_time, SEP);
    fprintf(myfile, "%2.3f%c", info_allcuts.root_time, SEP);
    fprintf(myfile, "%2.3f%c", info_mycuts.last_sol_time, SEP);
    fprintf(myfile, "%2.3f%c", info_allcuts.last_sol_time, SEP);
    fprintf(myfile, "%2.3f%c", info_mycuts.time, SEP);
    fprintf(myfile, "%2.3f%c", info_allcuts.time, SEP);
  } else {
    for (int i = 0; i < (int) BB_INFO_CONTENTS.size() * 2; i++) {
      fprintf(myfile, "%c", SEP);
    }
  }
} /* printBBInfo */

void createStringFromBBInfoVec(const std::vector<BBInfo>& vec_info,
    std::vector<std::string>& vec_str) {
  vec_str.resize(BB_INFO_CONTENTS.size());
  for (BBInfo info : vec_info) {
    vec_str[OBJ_BB_INFO_IND] += (!vec_str[OBJ_BB_INFO_IND].empty() ? ";" : "") + stringValue(info.obj, "%.20f");
    vec_str[BOUND_BB_INFO_IND] += (!vec_str[BOUND_BB_INFO_IND].empty() ? ";" : "") + stringValue(info.bound, "%.20f");
    vec_str[ITERS_BB_INFO_IND] += (!vec_str[ITERS_BB_INFO_IND].empty() ? ";" : "") + std::to_string(info.iters);
    vec_str[NODES_BB_INFO_IND] += (!vec_str[NODES_BB_INFO_IND].empty() ? ";" : "") + std::to_string(info.nodes);
    vec_str[ROOT_PASSES_BB_INFO_IND] += (!vec_str[ROOT_PASSES_BB_INFO_IND].empty() ? ";" : "") + std::to_string(info.root_passes);
    vec_str[FIRST_CUT_PASS_BB_INFO_IND] += (!vec_str[FIRST_CUT_PASS_BB_INFO_IND].empty() ? ";" : "") + std::to_string(info.first_cut_pass);
    vec_str[LAST_CUT_PASS_BB_INFO_IND] += (!vec_str[LAST_CUT_PASS_BB_INFO_IND].empty() ? ";" : "") + std::to_string(info.last_cut_pass);
    vec_str[ROOT_TIME_BB_INFO_IND] += (!vec_str[ROOT_TIME_BB_INFO_IND].empty() ? ";" : "") + std::to_string(info.root_time);
    vec_str[LAST_SOL_TIME_BB_INFO_IND] += (!vec_str[LAST_SOL_TIME_BB_INFO_IND].empty() ? ";" : "") + std::to_string(info.last_sol_time);
    vec_str[TIME_BB_INFO_IND] += (!vec_str[TIME_BB_INFO_IND].empty() ? ";" : "") + std::to_string(info.time);
  }
} /* createStringFromBBInfoVec */

/**
 * Creates temporary file (in /tmp) so that it can be read by a different solver
 * It does not delete the file
 */
void createTmpFileCopy(const VPCParameters& params,
    const OsiSolverInterface* const solver, std::string& f_name) {
  // Generate temporary file name
  char template_name[] = "/tmp/tmpmpsXXXXXX";

  mktemp(template_name);
  f_name = template_name;
  if (f_name.empty()) {
    error_msg(errorstring, "Could not generate temp file.\n");
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }
  solver->writeMps(template_name, "mps", solver->getObjSense());
  f_name += ".mps.gz"; // writeMps calls writeMpsNative, which invokes the CoinMpsIO writer with gzip option = 1
  //doBranchAndBoundWithCplex(f_name.c_str(), bb_opt, bb_iters, bb_nodes, bb_time);
  //remove(f_name.c_str());
} /* createTmpFileCopy (Osi) */

#ifdef USE_CBC
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
void generatePartialBBTree(PartialBBDisjunction* const owner, CbcModel* cbc_model,
    const OsiSolverInterface* const solver, const int max_leaf_nodes,
    const int num_strong, const int num_before_trusted) {
  //  const double partial_timelimit = 100 * max_leaf_nodes * solver->getNumCols()
  //      * GlobalVariables::timeStats.get_time(GlobalConstants::INIT_SOLVE_TIME);
  const double partial_timelimit = owner->params.get(PARTIAL_BB_TIMELIMIT); // will be checked manually by the eventHandler

  // Set up options
  VPCEventHandler* eventHandler = new VPCEventHandler(owner, max_leaf_nodes, partial_timelimit);
  eventHandler->setOriginalSolver(solver);

  // This sets branching decision, event handling, etc.
  setCbcParametersForPartialBB(owner->params, cbc_model, eventHandler, num_strong,
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

void setStrategyForBBTestCbc(const VPCParameters& params,
    CbcModel* const cbc_model,
    int seed = -1) {
  if (seed < 0) seed = params.get(intConst::RANDOM_SEED);
  // Parameters that should always be set
  cbc_model->setMaximumSeconds(params.get(doubleConst::BB_TIMELIMIT)); // time limit
  cbc_model->setRandomSeed(seed); // random seed

  int strategy = params.get(intParam::BB_STRATEGY);
  if (strategy <= 0) {
    // Default strategy
    CbcStrategyDefault strategy;
    cbc_model->setStrategy(strategy);
    /*
    CbcStrategyDefault strategy(-1);
    strategy.setupPreProcessing(-1,0);
    cbc_model->setStrategy(strategy);
    cbc_model->setMaximumCutPassesAtRoot(0);
    cbc_model->setMaximumCutPasses(0);
    cbc_model->setWhenCuts(0);
    */
  } else {
    if (use_bb_option(strategy, BB_Strategy_Options::user_cuts)) {
      // Make sure dual reductions are off
    }

    // Turn off all cuts
    if (use_bb_option(strategy, BB_Strategy_Options::all_cuts_off)) {
      cbc_model->setMaximumCutPassesAtRoot(0);
      cbc_model->setMaximumCutPasses(0);
      cbc_model->setWhenCuts(0);
    }

    // Presolve
    if (use_bb_option(strategy, BB_Strategy_Options::presolve_off)) {
      cbc_model->setTypePresolve(0);
    }
  }

  // Check if we should use strong branching
  // Not sure this works when using StrategyDefault as well...
  if (use_bb_option(std::abs(strategy), BB_Strategy_Options::strong_branching_on)) {
    CbcBranchDefaultDecision branch;
    OsiChooseStrong choose;
    choose.setNumberStrong(cbc_model->solver()->getNumCols());
    choose.setNumberBeforeTrusted(std::numeric_limits<int>::max());
    branch.setChooseMethod(choose);
    cbc_model->setBranchingMethod(branch);
  }
} /* setStrategyForBBTestCbc */

/************************************************************/
/**
 * Perform branch-and-bound without cuts
 */
void doBranchAndBoundNoCuts(const VPCParameters& params,
    const OsiSolverInterface* const solver, BBInfo& info) {
#ifdef TRACE
  printf("\nBB with no cuts.\n");
#endif
  const bool test_using_main = false;

  // Set up solver
  SolverInterface* BBSolver;
  BBSolver = dynamic_cast<SolverInterface*>(solver->clone());
  setLPSolverParameters(BBSolver);

  // Set up model
  CbcModel* cbc_model = new CbcModel;
  cbc_model->swapSolver(BBSolver);
  cbc_model->setModelOwnsSolver(true);
  setIPSolverParameters(cbc_model, params.get(VERBOSITY));

  // Set up options and run B&B
  if (!test_using_main) {
    setStrategyForBBTestCbc(params, cbc_model);
#ifdef TRACE
    cbc_model->branchAndBound(3);
#else
    cbc_model->branchAndBound(0);
#endif
  } else {
    CbcMain0(*cbc_model);
    std::string name, logLevel, presolveOnOff, preprocessOnOff, cutsOnOff, heurOnOff, solveOption;
    name = "BBHelper_doBranchAndBoundNoCuts";
    presolveOnOff = "-presolve=off";
    preprocessOnOff = "-preprocess=off";
    cutsOnOff = "-cuts=off";
    heurOnOff = "-heur=off";
    solveOption = "-solve";
#ifdef TRACE
    logLevel = "-loglevel=3";
#else
    logLevel = "-loglevel=0";
#endif

    int argc = 0;
    const char** cbc_options = new const char*[20];
    cbc_options[argc++] = name.c_str();
    cbc_options[argc++] = logLevel.c_str();
    cbc_options[argc++] = presolveOnOff.c_str();
    cbc_options[argc++] = preprocessOnOff.c_str();
    cbc_options[argc++] = cutsOnOff.c_str();
    cbc_options[argc++] = heurOnOff.c_str();
    cbc_options[argc++] = solveOption.c_str();

    CbcMain1(argc, cbc_options, *cbc_model);
    delete[] cbc_options;
  }

  // Collect statistics
  info.time = CoinCpuTime()
      - cbc_model->getDblParam(CbcModel::CbcStartSeconds);
  info.obj = cbc_model->getObjValue();
  info.iters = cbc_model->getIterationCount();
  info.nodes = cbc_model->getNodeCount();

  // Free
  if (cbc_model) {
    delete cbc_model;
    cbc_model = NULL;
  }
} /* doBranchAndBoundNoCuts */

/************************************************************/
/**
 * Perform branch-and-bound using the given cuts, perhaps doing cut selection
 */
void doBranchAndBoundYesCuts(const VPCParameters& params,
    const OsiSolverInterface* const solver, BBInfo& info, const OsiCuts& structCuts,
    const bool doCutSelection, const int numCutsToAddPerRound,
    const int maxRounds, const std::string logstring) {
//#ifdef TRACE
  printf("%s", logstring.c_str());
//#endif
  const bool test_using_main = false;

  // Check that there are cuts
  const int numCuts = structCuts.sizeCuts();
  if (numCuts == 0) {
    info.nodes = 0;
    info.time = 0.;
    return;
  }

  // Set up solver
  SolverInterface* BBSolver;
  BBSolver = dynamic_cast<SolverInterface*>(solver->clone());
  setLPSolverParameters(BBSolver, params.get(VERBOSITY));

  // Apply cuts
//  if (doCutSelection) {
//    applyCutsInRounds(BBSolver, structCuts, numCutsToAddPerRound, maxRounds);
//  } else {
    applyCutsCustom(BBSolver, structCuts);
//  }

  // Set up model
  CbcModel* cbc_model = new CbcModel;
  cbc_model->swapSolver(BBSolver);
  cbc_model->setModelOwnsSolver(true);
  setIPSolverParameters(cbc_model, params.get(VERBOSITY));

  // Set up options and run B&B
  if (!test_using_main) {
    setStrategyForBBTestCbc(params, cbc_model);
#ifdef TRACE
    cbc_model->branchAndBound(3);
#else
    cbc_model->branchAndBound(0);
#endif
  } else {
    CbcMain0(*cbc_model);
    std::string name, logLevel, presolveOnOff, preprocessOnOff, cutsOnOff, heurOnOff, solveOption;
    name = "BBHelper_doBranchAndBoundYesCuts";
    presolveOnOff = "-presolve=off";
    preprocessOnOff = "-preprocess=off";
    cutsOnOff = "-cuts=off";
    heurOnOff = "-heur=off";
    solveOption = "-solve";
#ifdef TRACE
    logLevel = "-loglevel=3";
#else
    logLevel = "-loglevel=0";
#endif

    int argc = 0;
    const char** cbc_options = new const char*[20];
    cbc_options[argc++] = name.c_str();
    cbc_options[argc++] = logLevel.c_str();
    cbc_options[argc++] = presolveOnOff.c_str();
    cbc_options[argc++] = preprocessOnOff.c_str();
    cbc_options[argc++] = cutsOnOff.c_str();
    cbc_options[argc++] = heurOnOff.c_str();
    cbc_options[argc++] = solveOption.c_str();

    CbcMain1(argc, cbc_options, *cbc_model);
    delete[] cbc_options;
  }

  // Collect statistics
  info.time = CoinCpuTime()
      - cbc_model->getDblParam(CbcModel::CbcStartSeconds);
  info.obj = cbc_model->getObjValue();
  info.iters = cbc_model->getIterationCount();
  info.nodes = cbc_model->getNodeCount();

  // Free
  if (cbc_model) {
    delete cbc_model;
    cbc_model = NULL;
  }
} /* doBranchAndBoundYesCuts */
#endif // USE_CBC

/**
 * amountToPrint:
 *  0 = all (newline-separated),
 *  1 = only header (comma-separated),
 *  2 = only values (comma-separated)
 */
void writeBBInforToLog(const SummaryBBInfo& info_mycuts,
    const SummaryBBInfo& info_allcuts, FILE *myfile, const int amountToPrint,
    const char SEP) {
  if (myfile == NULL) {
    return;
  }
  const int countBBInfoEntries = (int) BB_INFO_CONTENTS.size() * 4 * 2;

  //////////////////// BB INFO
  switch (amountToPrint) {
    case 0: {
//      for (std::string name : BB_INFO_CONTENTS) {
//        fprintf(myfile, "%s%c", ("FIRST BB " + name + " VPC").c_str(), SEP);
//        fprintf(myfile, "%s%c", ("FIRST BB " + name + " VPC+GOMORY").c_str(), SEP);
//      }
//      for (std::string name : BB_INFO_CONTENTS) {
//        fprintf(myfile, "%s%c", ("BEST BB " + name + " VPC").c_str(), SEP);
//        fprintf(myfile, "%s%c", ("BEST BB " + name + " VPC+GOMORY").c_str(), SEP);
//      }
//      for (std::string name : BB_INFO_CONTENTS) {
//        fprintf(myfile, "%s%c", ("AVG BB " + name + " VPC").c_str(), SEP);
//        fprintf(myfile, "%s%c", ("AVG BB " + name + " VPC+GOMORY").c_str(), SEP);
//      }
//      for (std::string name : BB_INFO_CONTENTS) {
//        fprintf(myfile, "%s%c", ("ALL BB " + name + " VPC").c_str(), SEP);
//        fprintf(myfile, "%s%c", ("ALL BB " + name + " VPC+GOMORY").c_str(), SEP);
//      }
      break;
    }
    case 1: {
      for (std::string name : BB_INFO_CONTENTS) {
        fprintf(myfile, "%s%c", ("FIRST BB " + name + " VPC").c_str(), SEP);
        fprintf(myfile, "%s%c", ("FIRST BB " + name + " VPC+GOMORY").c_str(), SEP);
      }
      for (std::string name : BB_INFO_CONTENTS) {
        fprintf(myfile, "%s%c", ("BEST BB " + name + " VPC").c_str(), SEP);
        fprintf(myfile, "%s%c", ("BEST BB " + name + " VPC+GOMORY").c_str(), SEP);
      }
      for (std::string name : BB_INFO_CONTENTS) {
        fprintf(myfile, "%s%c", ("AVG BB " + name + " VPC").c_str(), SEP);
        fprintf(myfile, "%s%c", ("AVG BB " + name + " VPC+GOMORY").c_str(), SEP);
      }
      for (std::string name : BB_INFO_CONTENTS) {
        fprintf(myfile, "%s%c", ("ALL BB " + name + " VPC").c_str(), SEP);
        fprintf(myfile, "%s%c", ("ALL BB " + name + " VPC+GOMORY").c_str(), SEP);
      }
      break;
    }
    case 2: {
      if (info_mycuts.vec_bb_info.size() == 0) {
        for (int i = 0; i < countBBInfoEntries; i++) {
          fprintf(myfile, "%c", SEP);
        }
        fflush(myfile);
        break;
      }
      // First
      printBBInfo(info_mycuts.vec_bb_info[0], info_allcuts.vec_bb_info[0], myfile, false, SEP);

      // Min
      printBBInfo(info_mycuts.best_bb_info, info_allcuts.best_bb_info, myfile, false, SEP);

      // Average
      printBBInfo(info_mycuts.avg_bb_info, info_allcuts.avg_bb_info, myfile, false, SEP);

      // All
      std::vector<std::string> vec_str_mycuts, vec_str_allcuts;
      createStringFromBBInfoVec(info_mycuts.vec_bb_info, vec_str_mycuts);
      createStringFromBBInfoVec(info_allcuts.vec_bb_info, vec_str_allcuts);
      for (int i = 0; i < (int) vec_str_mycuts.size(); i++) {
        fprintf(myfile, "%s%c", vec_str_mycuts[i].c_str(), SEP);
        fprintf(myfile, "%s%c", vec_str_allcuts[i].c_str(), SEP);
      }
      break;
    }
    default: {
      // nothing
    }
  }
  fflush(myfile);
} /* writeBBInforToLog */
