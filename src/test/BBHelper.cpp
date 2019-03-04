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
#include <CbcModel.hpp>

// Variable selection
#include <CbcBranchDefaultDecision.hpp>
#include <OsiChooseVariable.hpp>

// General strategy
#include <CbcStrategy.hpp>
#endif /* USE_CBC */

#ifdef USE_GUROBI
#include "GurobiHelper.hpp"
#endif

void runBBTests(const VPCParameters& params, SummaryBBInfo& info_nocuts,
    SummaryBBInfo& info_mycuts, SummaryBBInfo& info_allcuts,
    const std::string fullfilename, OsiSolverInterface* const solver,
    const double best_bound, const OsiCuts* vpcs, const OsiCuts* const gmics) {
  const int num_vpcs = (vpcs != NULL) ? vpcs->sizeCuts() : 0;
  const int num_bb_runs = std::abs(params.get(intParam::BB_RUNS))
      * ((params.get(intParam::DISJ_TERMS) == 0) || (num_vpcs > 0));
  if (num_bb_runs == 0)
    return;

  /***********************************************************************************
   * Perform branch-and-bound
   ***********************************************************************************/
  const int num_gmics = (gmics != NULL) ? gmics->sizeCuts() : 0;
  const int numCutsToAddPerRound = (num_gmics > 0) ? num_gmics : num_vpcs;

  // B&B mode: ones bit = no_cuts, tens bit = w/vpcs, hundreds bit = w/gmics
  const int mode_param = params.get(intParam::BB_MODE);
  const int mode_ones = mode_param % 10;
  const int mode_tens = (mode_param % 100 - (mode_param % 10)) / 10;
  const int mode_hundreds = (mode_param % 1000 - (mode_param % 100)) / 100;
  const bool branch_with_no_cuts = (mode_ones > 0);
  const bool branch_with_vpcs = (mode_tens > 0) && (num_vpcs > 0);
  const bool branch_with_gmics = (mode_hundreds > 0) && (num_gmics > 0);
  if (branch_with_no_cuts + branch_with_vpcs + branch_with_gmics == 0)
    return;

  const bool should_permute_rows_and_cols = (num_bb_runs >= 2)
      && !use_bb_option(params.get(BB_STRATEGY), BB_Strategy_Options::gurobi)
      && !use_bb_option(params.get(BB_STRATEGY), BB_Strategy_Options::cplex);

#ifdef TRACE
  printf("\n## Performing branch-and-bound tests. ##\n");
#endif

  info_nocuts.num_cuts = 0;
  info_mycuts.num_cuts = num_vpcs;
  info_allcuts.num_cuts = num_vpcs + num_gmics;
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

  OsiCuts GMICsAndVPCs;
  if (branch_with_vpcs)
    GMICsAndVPCs.insert(*vpcs);
  if (branch_with_gmics) {
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
      if (params.get(TEMP) == static_cast<int>(TempOptions::CHECK_CUTS_AGAINST_BB_OPT) && num_vpcs > 0) {
        // Get the original solution
        BBInfo tmp_bb_info;
        std::vector<double> solution;
        doBranchAndBoundWithGurobi(params, params.get(BB_STRATEGY),
            fullfilename.c_str(), tmp_bb_info, best_bound, &solution);

        // Check cuts
        for (int cut_ind = 0; cut_ind < num_vpcs; cut_ind++) {
          OsiRowCut currCut = vpcs->rowCut(cut_ind);
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
      if (branch_with_no_cuts) {
        doBranchAndBoundWithUserCutsGurobi(params, params.get(BB_STRATEGY),
            fullfilename.c_str(), NULL, info_nocuts.vec_bb_info[run_ind],
            best_bound);
      }
      if (branch_with_vpcs) {
        doBranchAndBoundWithUserCutsGurobi(params, params.get(BB_STRATEGY),
            fullfilename.c_str(), vpcs, info_mycuts.vec_bb_info[run_ind],
            best_bound);
      }
      if (branch_with_gmics) {
        doBranchAndBoundWithUserCutsGurobi(params, params.get(BB_STRATEGY),
            fullfilename.c_str(), &GMICsAndVPCs,
            info_allcuts.vec_bb_info[run_ind], best_bound);
      }
#endif // USE_GUROBI
    } else if (use_bb_option(params.get(BB_STRATEGY), BB_Strategy_Options::cplex)) {
//#ifdef USE_CPLEX
//      doBranchAndBoundWithUserCutsCplexCallable(
//          GlobalVariables::in_f_name.c_str(), &structVPC, vec_bb_info_mycuts[run_ind]);
//      doBranchAndBoundWithUserCutsCplexCallable(
//          GlobalVariables::in_f_name.c_str(), &SICsAndVPCs,
//          vec_bb_info_allcuts[run_ind]);
//#endif // USE_CPLEX
    } else if (use_bb_option(params.get(BB_STRATEGY), BB_Strategy_Options::cbc)) {
      if (branch_with_no_cuts) {
        doBranchAndBoundNoCuts(params, runSolver, info_nocuts.vec_bb_info[run_ind]);
      }
      if (branch_with_vpcs) {
        doBranchAndBoundYesCuts(params, runSolver, info_mycuts.vec_bb_info[run_ind],
            *vpcs, false, numCutsToAddPerRound, 1, "\nBB with VPCs.\n");
      }
      if (branch_with_gmics) {
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
  if (!myfile)
    return;
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
  if (!myfile)
    return;
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
