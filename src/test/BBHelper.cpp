/**
 * @file BBHelper.cpp
 * @brief Helper functions for branch-and-bound
 * @author A. M. Kazachkov
 * @date 2018-12-24
 */
#include <cstdio> // for tmpnam
#include <chrono>
#include <random> // random_device, mt19937 (mersenne twister engine), uniform_real_distribution

// Project files 
#include "BBHelper.hpp"
#include "CutHelper.hpp" // applyCutsCustom
#include "PartialBBDisjunction.hpp"
#include "SolverHelper.hpp"
#include "VPCEventHandler.hpp"
#include "VPCParameters.hpp"
using namespace VPCParametersNamespace;
#include "TimeStats.hpp"
#include "utility.hpp"

// COIN-OR
#include <CoinTime.hpp>
#include <OsiCuts.hpp>

#ifdef USE_CBC
#include "CbcHelper.hpp"
#endif // USE_CBC
#ifdef USE_GUROBI
#include "GurobiHelper.hpp"
#endif // USE_GUROBI
#ifdef USE_CPLEX
#include "CplexHelper.hpp"
#endif // USE_CPLEX

/**
 * @details Solver is unchanged unless we are doing row/col permutations
 */
void runBBTests(
    /// [in] Initial set of parameters (will be cloned to change random seeds)
    const VPCParametersNamespace::VPCParameters& base_params,
    /// [out] Store information about branch-and-bound results with no (user-generated) cuts added
    SummaryBBInfo* const info_nocuts,
    /// [out] Store information about branch-and-bound results with only VPCs added among user-generated cuts
    SummaryBBInfo* const info_mycuts,
    /// [out] Store information about branch-and-bound results with all user-generated cuts added (including GMICs, for example)
    SummaryBBInfo* const info_allcuts,
    /// [in] Full path to instance
    const std::string fullfilename,
    /// [in] Original model
    const OsiSolverInterface* const solver,
    /// [in] Best bound (to provide to IP solver, if relevant option is selected)
    const double best_bound,
    /// [in] VPCs to add to solver
    const OsiCuts* vpcs,
    /// [in] (Optional) GMICs to test against
    const OsiCuts* const gmics,
    /// [in,out] IP solution to original problem (will usually be empty unless you are debugging; will be computed here if not previously computed and is requested)
    std::vector<double>* const ip_solution) {
  VPCParametersNamespace::VPCParameters params = base_params;
  const int num_vpcs = (vpcs != NULL) ? vpcs->sizeCuts() : 0;
  // Set number of b&b runs to be zero if no cuts generated (unless no disjunctions were requested in the first place)
  const int num_bb_runs = std::abs(params.get(intParam::BB_RUNS))
      * ((params.get(intParam::DISJ_TERMS) == 0) || (num_vpcs > 0));
  if (num_bb_runs == 0)
    return;

  /***********************************************************************************
   * Perform branch-and-bound
   ***********************************************************************************/
  const int num_gmics = (gmics != NULL) ? gmics->sizeCuts() : 0;
  //const int numCutsToAddPerRound = (num_gmics > 0) ? num_gmics : num_vpcs;

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

  printf("\n## Performing branch-and-bound tests. ##\n");

  if (branch_with_no_cuts) {
    assert(info_nocuts != NULL);
    info_nocuts->vec_bb_info.resize(0);
    info_nocuts->vec_bb_info.resize(num_bb_runs);
  }
  if (branch_with_vpcs) {
    assert(info_mycuts != NULL);
    info_mycuts->vec_bb_info.resize(0);
    info_mycuts->num_cuts = num_vpcs;
    info_mycuts->vec_bb_info.resize(num_bb_runs);
  }
  if (branch_with_gmics) {
    assert(info_allcuts != NULL);
    info_allcuts->num_cuts = num_vpcs + num_gmics;
    info_allcuts->vec_bb_info.resize(num_bb_runs);
  }
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
  const int given_seed = params.get(intParam::RANDOM_SEED);
  const auto initial_random_seed = given_seed > 0 ? given_seed : std::chrono::system_clock::now().time_since_epoch().count();
  int random_seed;
  std::mt19937 rng;
  for (int run_ind = 0; run_ind < num_bb_runs; run_ind++) {
    // Change the random seed per run
    random_seed = initial_random_seed * (run_ind + 1);
    params.set(intParam::RANDOM_SEED, random_seed);

    // For Cbc, for every run after the first, we will randomize the rows and columns of the input
    if (should_permute_rows_and_cols) {
      //std::srand(run_ind + 2);  // sets the seed for pseudo-randomness in C++
      std::seed_seq curr_seed{ random_seed };
      rng.seed(curr_seed);
      std::vector<int> this_row_permutation(row_permutation);
      std::vector<int> this_col_permutation(col_permutation);
      std::shuffle(this_row_permutation.begin(), this_row_permutation.end(), rng);
      std::shuffle(this_col_permutation.begin(), this_col_permutation.end(), rng);

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
      setLPSolverParameters(runSolver, params.get(intParam::VERBOSITY), params.get(doubleParam::TIMELIMIT));
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

    // Do branch and bound
    if (use_bb_option(params.get(BB_STRATEGY), BB_Strategy_Options::gurobi)) {
#ifdef USE_GUROBI
      if (use_temp_option(params.get(TEMP), TempOptions::CHECK_CUTS_AGAINST_BB_OPT) && num_vpcs > 0) {
        if (ip_solution && !ip_solution->empty()) {
          // Get the original solution
          BBInfo tmp_bb_info;
          doBranchAndBoundWithGurobi(params, params.get(BB_STRATEGY),
              fullfilename.c_str(), tmp_bb_info, best_bound, ip_solution);
        }

        if (!ip_solution || ip_solution->empty()) {
          error_msg(errorstring, "Unable to obain IP solution.\n");
          writeErrorToLog(errorstring, params.logfile);
          exit(1);
        }

        // Check cuts
        for (int cut_ind = 0; cut_ind < num_vpcs; cut_ind++) {
          OsiRowCut currCut = vpcs->rowCut(cut_ind);
          const double rhs = currCut.rhs();
          const int num_el = currCut.row().getNumElements();
          const int* ind = currCut.row().getIndices();
          const double* el = currCut.row().getElements();
          const double activity = dotProduct(num_el, ind, el, ip_solution->data());

          if (lessThanVal(activity, rhs)) {
            warning_msg(warnstring, "Cut %d removes optimal solution. Activity: %.10f. Rhs: %.10f.\n", cut_ind, activity, rhs);
          }
        }
      } // checking cuts for violating the IP opt
      if (branch_with_no_cuts) {
        // NB: If user cuts is not explicitly set, this is WITHOUT user cuts
        doBranchAndBoundWithGurobi(params, params.get(BB_STRATEGY),
            fullfilename.c_str(), info_nocuts->vec_bb_info[run_ind],
            best_bound);
      }
      if (branch_with_vpcs) {
        doBranchAndBoundWithUserCutsGurobi(params, params.get(BB_STRATEGY),
            fullfilename.c_str(), vpcs, info_mycuts->vec_bb_info[run_ind],
            best_bound);
      }
      if (branch_with_gmics) {
        doBranchAndBoundWithUserCutsGurobi(params, params.get(BB_STRATEGY),
            fullfilename.c_str(), &GMICsAndVPCs,
            info_allcuts->vec_bb_info[run_ind], best_bound);
      }
#endif // USE_GUROBI
    } else if (use_bb_option(params.get(BB_STRATEGY), BB_Strategy_Options::cplex)) {
#ifdef USE_CPLEX
      if (use_temp_option(params.get(TEMP), TempOptions::CHECK_CUTS_AGAINST_BB_OPT) && num_vpcs > 0) {
        // Get the original solution
        BBInfo tmp_bb_info;
        std::vector<double> solution;
        doBranchAndBoundWithCplexCallable(params, params.get(BB_STRATEGY),
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
        // NB: If user cuts is not explicitly set, this is WITHOUT user cuts
        doBranchAndBoundWithCplexCallable(params, params.get(BB_STRATEGY),
            fullfilename.c_str(), info_nocuts->vec_bb_info[run_ind],
            best_bound);
      }
      if (branch_with_vpcs) {
        doBranchAndBoundWithUserCutsCplexCallable(params, params.get(BB_STRATEGY),
            fullfilename.c_str(), vpcs, info_mycuts->vec_bb_info[run_ind],
            best_bound);
      }
      if (branch_with_gmics) {
        doBranchAndBoundWithUserCutsCplexCallable(params, params.get(BB_STRATEGY),
            fullfilename.c_str(), &GMICsAndVPCs,
            info_allcuts->vec_bb_info[run_ind], best_bound);
      }
#endif // USE_CPLEX
    } else if (use_bb_option(params.get(BB_STRATEGY), BB_Strategy_Options::cbc)) {
#ifdef USE_CBC
      if (runSolver == NULL) {
        runSolver = const_cast<OsiSolverInterface*>(solver);
      }
      if (branch_with_no_cuts) {
        // NB: If user cuts is not explicitly set, this is WITHOUT user cuts
        doBranchAndBoundWithCbc(params, params.get(BB_STRATEGY),
            runSolver, info_nocuts->vec_bb_info[run_ind],
            best_bound);
      }
      if (branch_with_vpcs) {
        doBranchAndBoundWithUserCutsCbc(params, params.get(BB_STRATEGY),
            runSolver, vpcs, info_mycuts->vec_bb_info[run_ind],
            best_bound);
      }
      if (branch_with_gmics) {
        doBranchAndBoundWithUserCutsCbc(params, params.get(BB_STRATEGY),
            runSolver, &GMICsAndVPCs,
            info_allcuts->vec_bb_info[run_ind], best_bound);
      }
      /*if (branch_with_no_cuts) {
        doBranchAndBoundNoCuts(params, (should_permute_rows_and_cols ? runSolver : solver), info_nocuts->vec_bb_info[run_ind]);
      }
      if (branch_with_vpcs) {
        doBranchAndBoundYesCuts(params, (should_permute_rows_and_cols ? runSolver : solver), info_mycuts->vec_bb_info[run_ind],
            *vpcs, false, numCutsToAddPerRound, 1, "\nBB with VPCs.\n");
      }
      if (branch_with_gmics) {
        doBranchAndBoundYesCuts(params, (should_permute_rows_and_cols ? runSolver : solver), info_allcuts->vec_bb_info[run_ind],
            GMICsAndVPCs, false, numCutsToAddPerRound, 1, "\nBB with VPC+Gomorys.\n");
      }*/
#endif // USE_CBC
    }

    if (branch_with_no_cuts) updateBestBBInfo(info_nocuts->best_bb_info, info_nocuts->vec_bb_info[run_ind], (run_ind == 0));
    if (branch_with_vpcs) updateBestBBInfo(info_mycuts->best_bb_info, info_mycuts->vec_bb_info[run_ind], (run_ind == 0));
    if (branch_with_gmics) updateBestBBInfo(info_allcuts->best_bb_info, info_allcuts->vec_bb_info[run_ind], (run_ind == 0));

    // Free memory if necessary
    if (should_permute_rows_and_cols && runSolver && (runSolver != solver)) {
      delete runSolver;
    }
  } /* end iterating over runs */

  if (branch_with_no_cuts) {
    info_nocuts->first_bb_info = info_nocuts->vec_bb_info[0];
    averageBBInfo(info_nocuts->avg_bb_info, info_nocuts->vec_bb_info);
  }
  if (branch_with_vpcs) {
    info_mycuts->first_bb_info = info_mycuts->vec_bb_info[0];
    averageBBInfo(info_mycuts->avg_bb_info, info_mycuts->vec_bb_info);
  }
  if (branch_with_gmics) {
    info_allcuts->first_bb_info = info_allcuts->vec_bb_info[0];
    averageBBInfo(info_allcuts->avg_bb_info, info_allcuts->vec_bb_info);
  }
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
  best_info.root_iters = first ? curr_info.root_iters : CoinMin(best_info.root_iters, curr_info.root_iters);
  best_info.root_time = first ? curr_info.root_time : CoinMin(best_info.root_time, curr_info.root_time);
  best_info.last_sol_time = first ? curr_info.last_sol_time : CoinMin(best_info.last_sol_time, curr_info.last_sol_time);
  best_info.time = first ? curr_info.time : CoinMin(best_info.time, curr_info.time);
} /* updateMinBBInfo */

void averageBBInfo(BBInfo& avg_info, const std::vector<BBInfo>& info) {
  for (const BBInfo& curr_info : info) {
    avg_info.obj += curr_info.obj;
    avg_info.bound += curr_info.bound;
    avg_info.iters += curr_info.iters;
    avg_info.nodes += curr_info.nodes;
    avg_info.root_passes += curr_info.root_passes;
    avg_info.first_cut_pass += curr_info.first_cut_pass;
    avg_info.last_cut_pass += curr_info.last_cut_pass;
    avg_info.root_iters += curr_info.root_iters;
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
  avg_info.root_iters /= num_bb_runs;
  avg_info.root_time /= num_bb_runs;
  avg_info.last_sol_time /= num_bb_runs;
  avg_info.time /= num_bb_runs;
} /* averageBBInfo */

void createStringFromBBInfoVec(const std::vector<BBInfo>& vec_info,
    std::vector<std::string>& vec_str) {
  vec_str.resize(BB_INFO_CONTENTS.size());
  for (const BBInfo& info : vec_info) {
    vec_str[OBJ_BB_INFO_IND] += (!vec_str[OBJ_BB_INFO_IND].empty() ? ";" : "") + stringValue(info.obj, "%.20f");
    vec_str[BOUND_BB_INFO_IND] += (!vec_str[BOUND_BB_INFO_IND].empty() ? ";" : "") + stringValue(info.bound, "%.20f");
    vec_str[ITERS_BB_INFO_IND] += (!vec_str[ITERS_BB_INFO_IND].empty() ? ";" : "") + std::to_string(info.iters);
    vec_str[NODES_BB_INFO_IND] += (!vec_str[NODES_BB_INFO_IND].empty() ? ";" : "") + std::to_string(info.nodes);
    vec_str[ROOT_PASSES_BB_INFO_IND] += (!vec_str[ROOT_PASSES_BB_INFO_IND].empty() ? ";" : "") + std::to_string(info.root_passes);
    vec_str[FIRST_CUT_PASS_BB_INFO_IND] += (!vec_str[FIRST_CUT_PASS_BB_INFO_IND].empty() ? ";" : "") + std::to_string(info.first_cut_pass);
    vec_str[LAST_CUT_PASS_BB_INFO_IND] += (!vec_str[LAST_CUT_PASS_BB_INFO_IND].empty() ? ";" : "") + std::to_string(info.last_cut_pass);
    vec_str[ROOT_ITERS_BB_INFO_IND] += (!vec_str[ROOT_ITERS_BB_INFO_IND].empty() ? ";" : "") + std::to_string(info.root_iters);
    vec_str[ROOT_TIME_BB_INFO_IND] += (!vec_str[ROOT_TIME_BB_INFO_IND].empty() ? ";" : "") + std::to_string(info.root_time);
    vec_str[LAST_SOL_TIME_BB_INFO_IND] += (!vec_str[LAST_SOL_TIME_BB_INFO_IND].empty() ? ";" : "") + std::to_string(info.last_sol_time);
    vec_str[TIME_BB_INFO_IND] += (!vec_str[TIME_BB_INFO_IND].empty() ? ";" : "") + std::to_string(info.time);
  }
} /* createStringFromBBInfoVec */

/**
 * Creates temporary file (in /tmp) so that it can be read by a different solver
 * It does not delete the file
 */
void createTmpFileCopy(
    const VPCParametersNamespace::VPCParameters& params,
    const OsiSolverInterface* const solver,
    std::string& f_name) {
  if (f_name.empty()) {
    try {
      createTmpFilename(f_name, "", params.get(stringParam::TMPFOLDER));
    } catch (const std::exception &e) {
      error_msg(errorstring, "Could not generate temp file: %s.\n", e.what());
      writeErrorToLog(errorstring, params.logfile);
      exit(1);
    }
  }

  solver->writeMps(f_name.c_str(), "mps", solver->getObjSense());
  f_name += ".mps"; // writeMps calls writeMpsNative, which invokes the CoinMpsIO writer with gzip option = 1
} /* createTmpFileCopy (Osi) */

