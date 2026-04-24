/**
 * @file main.cpp
 * @author A. M. Kazachkov
 * @date 2018-12-24
 */

// V-Polyhedral Disjunctive Cuts
// Author:   A. M. Kazachkov
// Date:     2018-Dec-24
//-----------------------------------------------------------------------------

#include <iostream>
#include <cstdio>
#include <string>
#include <vector>
#include <chrono> // for timing
#include <limits> // numeric_limits
#include <memory> // for smart pointers

// For option handling
#include <getopt.h> // getopt, getopt_long, no_argument = 0, required_argument = 1, optional_argument = 2

// COIN-OR
#include <OsiCuts.hpp>
#include <CglGMI.hpp>

// Project files
#include "analysis.hpp" // analyzeStrength, analalyzeBB
#include "BBHelper.hpp" // runBBTests
#include "CglVPC.hpp"
#include "CutHelper.hpp"
#include "Disjunction.hpp" // needed to access disjunction properties
#include "DisjunctionHelper.hpp" // for custom disjunctions
#include "main_input.hpp"
#include "main_round_actions.hpp"
#include "main_round_orchestrator.hpp"
#include "main_reporting.hpp"
#include "main_startup.hpp"
#include "preprocess.hpp" // performCleaning
#include "SolverHelper.hpp" // initializeSolver
#include "SolverInterface.hpp"
#include "VPCParameters.hpp"
using namespace VPCParametersNamespace;
#include "TimeStats.hpp"
#include "utility.hpp"

#ifdef CALC_COND_NUM
#include "condition_number.hpp" // compute_condition_number
#endif

#ifdef USE_GUROBI
#include "GurobiHelper.hpp" // for obtaining ip opt
#endif

#ifdef VPC_DEBUG
#include "debug.hpp"
#endif

enum OverallTimeStats {
  INIT_SOLVE_TIME,
  GOMORY_GEN_TIME,
  GOMORY_APPLY_TIME,
  CUT_TOTAL_TIME,
  VPC_GEN_TIME,
  VPC_APPLY_TIME,
  VPC_INIT_SOLVE_TIME,
  VPC_DISJ_SETUP_TIME,
  VPC_DISJ_GEN_TIME,
  VPC_PRLP_SETUP_TIME,
  VPC_PRLP_SOLVE_TIME,
  VPC_GEN_CUTS_TIME,
  BB_TIME,
  TOTAL_APPLY_TIME,
  TOTAL_TIME,
  NUM_TIME_STATS
}; /* OverallTimeStats */
const std::vector<std::string> OverallTimeStatsName {
  "INIT_SOLVE_TIME",
  "GOMORY_GEN_TIME",
  "GOMORY_APPLY_TIME",
  "CUT_TOTAL_TIME",
  "VPC_GEN_TIME",
  "VPC_APPLY_TIME",
  "VPC_INIT_SOLVE_TIME",
  "VPC_DISJ_SETUP_TIME",
  "VPC_DISJ_GEN_TIME",
  "VPC_PRLP_SETUP_TIME",
  "VPC_PRLP_SOLVE_TIME",
  "VPC_GEN_CUTS_TIME",
  "BB_TIME",
  "TOTAL_APPLY_TIME",
  "TOTAL_TIME"
}; /* OverallTimeStatsName */

// Main file variables
VPCParameters params;

OsiSolverInterface *solver;               ///< stores the cuts we actually want to "count" (i.e., the base LP from which we generate cuts in each round)
OsiSolverInterface *origSolver;           ///< original solver in case we wish to come back to it later
OsiSolverInterface* GMICSolver = NULL;    ///< only GMICs
OsiSolverInterface* VPCSolver = NULL;     ///< if GMICs count, this is only VPCs; otherwise, it is both GMICs and VPCs
OsiSolverInterface* allCutsSolver = NULL; ///< all generated cuts (same as solver if no GMICs, in which case not generated)

OsiCuts gmics, vpcs;

std::string dir = "", filename_stub = "", instname = "", in_file_ext = "";

CglVPC::ExitReason exitReason;
TimeStats timer;
std::time_t start_time_t, end_time_t;
char start_time_string[25];

const int HOST_NAME_MAX = 1024;
char hostname[HOST_NAME_MAX];
// For CPU model and number
std::string cpu_model = "";
int cpu_id = -1;
#ifdef __linux__
#include <sched.h> // for cpu_id
#endif

SummaryBoundInfo boundInfo;
std::vector<SummaryBoundInfo> boundInfoVec;
SummaryBBInfo info_nocuts, info_mycuts, info_allcuts;
SummaryDisjunctionInfo disjInfo;
std::vector<SummaryCutInfo> cutInfoVec, cutInfoGMICVec;
SummaryCutInfo cutInfo, cutInfoGMICs;

// For output
std::string cut_output = "", bb_output = "";

// #ifdef CODE_VERSION
// const std::string CODE_VERSION_STRING = x_macro_to_string(CODE_VERSION);
// #endif
#ifdef VPC_VERSION
const std::string VPC_VERSION_STRING = x_macro_to_string(VPC_VERSION);
#endif
#ifdef VPC_CBC_VERSION
const std::string CBC_VERSION_STRING = x_macro_to_string(VPC_CBC_VERSION);
#endif
#ifdef VPC_CLP_VERSION
const std::string CLP_VERSION_STRING = x_macro_to_string(VPC_CLP_VERSION);
#endif
#ifdef USE_GUROBI
#include <gurobi_c++.h>
const std::string GUROBI_VERSION_STRING = std::to_string(GRB_VERSION_MAJOR) + "." + std::to_string(GRB_VERSION_MINOR) + std::to_string(GRB_VERSION_TECHNICAL);
#endif // USE_GUROBI
#ifdef USE_CPLEX
#include "ilcplex/cpxconst.h"
const std::string CPLEX_VERSION_STRING = std::to_string(CPX_VERSION_VERSION) + "." + std::to_string(CPX_VERSION_RELEASE) + "." + std::to_string(CPX_VERSION_MODIFICATION);
#endif // USE_CPLEX

// Catch abort signal if it ever gets sent
/**
 * Catch a signal. You can output debugging info.
 * If you return from this function, and it was called
 * because abort() was called, your program will exit or crash anyway
 * (with a dialog box on Windows).
 */
#include <csignal>
void signal_handler_with_error_msg(int signal_number) {
  std::string time_name = OverallTimeStatsName[OverallTimeStats::TOTAL_TIME];
  const int time_id = timer.get_id(time_name);
  if (time_id >= 0) {
    const double elapsed_time = timer.get_total_time(time_id);
    error_msg(errorstring, "Abort or seg fault message received. Signal number: %d. Elapsed time: %1.2f.\n", signal_number, elapsed_time);
    writeErrorToLog(errorstring, params.logfile);
  } else {
    error_msg(errorstring, "Abort or seg fault message received. Signal number: %d.\n", signal_number);
    writeErrorToLog(errorstring, params.logfile);
  }
  exit(1);
} /* signal_handler_with_error_msg */

int startUp(int argc, char** argv);
int processArgs(int argc, char** argv);
void doCustomRoundOfCuts(int round_ind, OsiCuts& vpcs, CglVPC& gen, int& num_disj);
int wrapUp(int retCode, int argc, char** argv);

/****************** MAIN FUNCTION **********************/
int main(int argc, char** argv) {
  // Do this early in your program's initialization
  std::signal(SIGABRT, signal_handler_with_error_msg);
  std::signal(SIGSEGV, signal_handler_with_error_msg);
  assert( OverallTimeStatsName.size() == OverallTimeStats::NUM_TIME_STATS );

  //====================================================================================================//
  // Set up timing
  for (int t = 0; t < OverallTimeStats::NUM_TIME_STATS; t++) {
    timer.register_name(OverallTimeStatsName[t]);
  }

  //====================================================================================================//
  // Print welcome message, set up logfile
  timer.start_timer(OverallTimeStats::TOTAL_TIME);
  int status = startUp(argc, argv);
  if (status) { return status; }

  //====================================================================================================//
  // Set up solver and get initial solution
  initializeSolver(solver, params.get(stringParam::FILENAME), params.get(VPCParametersNamespace::intParam::VERBOSITY), params.get(VPCParametersNamespace::doubleParam::TIMELIMIT), params.logfile);
  timer.start_timer(OverallTimeStats::INIT_SOLVE_TIME);
  solver->initialSolve();
  if (!checkSolverOptimality(solver, false)) {
    error_msg(errorstring, "Unable to solve initial LP relaxation.\n");
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }
  timer.end_timer(OverallTimeStats::INIT_SOLVE_TIME);
  boundInfo.lp_obj = solver->getObjValue();

  //====================================================================================================//
  // Save original solver in case we wish to come back to it later
  origSolver = solver->clone();
  if (!origSolver->isProvenOptimal()) {
    origSolver->initialSolve();
    checkSolverOptimality(origSolver, false);
  }

  // Also save copies for calculating other objective values
  // We only need these if cuts other than VPCs are generated
  // The "solver" is the one used for generating cuts for each round
  // solver        ::  stores the cuts we actually want to "count"
  // GMICSolver    ::  only GMICs
  // VPCSolver     ::  if GMICs count, this is only VPCs; otherwise, it is both GMICs and VPCs
  // allCutsSolver ::  all generated cuts
  if (params.get(GOMORY) != 0) {
    GMICSolver = solver->clone();
    VPCSolver = solver->clone();
    allCutsSolver = solver->clone();
  }

  //====================================================================================================//
  // Information from each round of cuts will be saved and optionally printed
  int num_rounds = params.get(ROUNDS); // not const in case we do not exhaust the limit
#ifndef CALC_COND_NUM
  const bool COMPUTE_COND_NUM = false;
#else
  const bool COMPUTE_COND_NUM = false;
#endif
  std::vector<OsiCuts> vpcs_by_round(num_rounds);
  cutInfoVec.resize(num_rounds);
  cutInfoGMICVec.resize(num_rounds);
  boundInfoVec.resize(num_rounds);

  std::vector<double> gmic_gen_time_by_round(num_rounds);
  std::vector<double> gmic_apply_time_by_round(num_rounds);
  //std::vector<double> gmic_cond_num1_by_round(num_rounds);
  std::vector<double> gmic_cond_num2_by_round;
  std::vector<double> vpc_gen_time_by_round(num_rounds);
  std::vector<double> vpc_apply_time_by_round(num_rounds);

  //std::vector<double> vpc_cond_num1_by_round(num_rounds);
  std::vector<double> vpc_cond_num2_by_round;
  if (COMPUTE_COND_NUM) {
    gmic_cond_num2_by_round.resize(num_rounds);
    vpc_cond_num2_by_round.resize(num_rounds);
  }

  //====================================================================================================//
  { // Check whether the initial solution is integer-feasible
    const double intTol = params.get(doubleParam::AWAY);
    fprintf(stdout, "\n## Checking integrality of initial solution using tolerance %g. ##\n", intTol);

    // Now check using solver's getFractionalIndices
    std::vector<int> fracCore = solver->getFractionalIndices(intTol);
    if (static_cast<int>(fracCore.size()) == 0) {
      boundInfo.ip_obj = boundInfo.lp_obj;
      params.set(doubleParam::IP_OBJ, boundInfo.ip_obj);
      warning_msg(warnstring, "solver->getFractionalIndices(%g): Initial solution with value %.6f is integer-feasible.\n", intTol, boundInfo.lp_obj);
    } else {
      // fprintf(stdout, "solver->getFractionalIndices(%g): Initial solution is not integer-feasible. %d fractional variables.\n", intTol, (int) fracCore.size());
    }

    // // Double check manually
    // std::vector<int> fracCore2;
    // for (int i = 0; i < solver->getNumCols(); i++) {
    //   if (solver->isInteger(i)) {
    //     if (!isZero(sol[i] - std::round(sol[i]), intTol)) {
    //       fracCore2.push_back(i);
    //     }
    //   }
    // }
    // if (static_cast<int>(fracCore2.size()) == 0) {
    //   boundInfo.ip_obj = boundInfo.lp_obj;
    //   params.set(doubleParam::IP_OBJ, boundInfo.ip_obj);
    //   warning_msg(warnstring, "Manual check: Initial solution is integer-feasible.\n");
    // } else {
    //   // fprintf(stdout, "Manual check: Initial solution is not integer-feasible. %d fractional variables.\n", (int) fracCore2.size());
    // }

    // If no integer variables are fractional, exit early
    if (static_cast<int>(fracCore.size()) == 0) {
      exitReason = CglVPC::ExitReason::OPTIMAL_SOLUTION_FOUND_EXIT;
      return wrapUp(0, argc, argv);
    }
  } // Check whether the initial solution is integer-feasible

  /** DEBUG TESTING BARRIER METHOD {
    OsiClpSolverInterface* interiorSolver = dynamic_cast<OsiClpSolverInterface*>(solver->clone());
    interiorSolver->getModelPtr()->barrier(false);

    if (interiorSolver->isProvenOptimal()) {
      std::cout << "Original\tBarrier\n";
      for (int i = 0; i < solver->getNumCols(); i++) {
        std::cout << solver->getColSolution()[i] << "\t";
        std::cout << interiorSolver->getColSolution()[i] << "\n";
      }
    } else {
      std::cerr << "Barrier method does not result in optimal solution." << std::endl;
    }
    exit(1);
  } **/

  //====================================================================================================//
  // Get IP solution if requested
  timer.end_timer(OverallTimeStats::TOTAL_TIME);
  std::vector<double> ip_solution;
  if (params.get(TEMP) == static_cast<int>(TempOptions::CHECK_CUTS_AGAINST_BB_OPT)) {
    printf("\n## Try to find optimal IP solution to compare against. ##\n");
    BBInfo tmp_bb_info;
#ifdef USE_GUROBI
    int strategy = params.get(BB_STRATEGY);
    // If no solfile provided, check if we can find it in the same directory as the instance
    if (params.get(SOLFILE).empty()) {
      std::string f_name = dir + "/" + instname + "_gurobi.mst.gz";
      if (fexists(f_name.c_str())) {
        printf("No solution provided, but solution file found at %s.\n", f_name.c_str());
        params.set(SOLFILE, f_name);
      }
    }
    // If solfile provided, enable use_bound
    if (!params.get(SOLFILE).empty() && !use_bb_option(strategy, VPCParametersNamespace::BB_Strategy_Options::use_best_bound)) {
      strategy = enable_bb_option(strategy, VPCParametersNamespace::BB_Strategy_Options::use_best_bound);
    }
    doBranchAndBoundWithGurobi(params, strategy,
        params.get(stringParam::FILENAME).c_str(),
        tmp_bb_info, boundInfo.ip_obj, &ip_solution);
#endif
#ifdef TRACE
    /*{/// DEBUG
    bool first = true;
    for (int i = 0; i < (int) ip_solution.size(); i++) {
      const double val = ip_solution[i];
      if (isZero(val)) continue;
      if (!first) printf(", ");
      else first = false;
      printf("(%d,%g)", i, val);
    }
    printf("\n");
    }*/ /// DEBUG
#endif
    if (ip_solution.size() > 0) {
      if (isInfinity(std::abs(boundInfo.ip_obj))) {
        const double ip_obj = dotProduct(ip_solution.data(), solver->getObjCoefficients(), solver->getNumCols());
        boundInfo.ip_obj = ip_obj;
        params.set(doubleParam::IP_OBJ, boundInfo.ip_obj);
        fprintf(stdout, "Best known IP objective value is %s.\n", stringValue(boundInfo.ip_obj, "%g").c_str());
      }
    }
  } // get IP opt
  timer.start_timer(OverallTimeStats::TOTAL_TIME);

  //====================================================================================================//
  // Possibly preprocess instead of doing cuts
  if (params.get(PREPROCESS) != 0) {
    // Cleaning involves running presolve and branching
    printf("\n## Starting preprocessing/cleaning of instance. ##\n");
    params.set(intParam::BB_MODE, 1); // only do no cuts branching
    performCleaning(params, solver, filename_stub, boundInfo.ip_obj, params.get(PREPROCESS));

    printf("\n## Finished preprocessing/cleaning of instance. ##\n");
    return wrapUp(0, argc, argv);
  }

  //====================================================================================================//
  // Process disjunction options, if different by round
  // Parse string DISJ_OPTIONS into vector of ints, splitting by default delimiter
  std::vector<int> disjOptions;
  if (parseDisjOptions(disjOptions, params)) {
    return wrapUp(1, argc, argv);
  }

  //====================================================================================================//
  // If requested, open cutrounds_logfile file then print header and round "0" information
  FILE* cutrounds_logfile = NULL;
  const int NUM_CUTROUND_INFO = 1 + 3 + 2*6 + 2*COMPUTE_COND_NUM; // round + 3 x bound + 2 x (num cuts, 3xdensity, gen time, apply time) + 2 x cond_num2 when computed
  if (use_temp_option(params.get(intParam::TEMP), TempOptions::PRINT_BOUND_BY_ROUND)) {
    std::string fileWithCuts = "cutrounds";
    const std::string logname = params.get(stringParam::LOGFILE);
    if (!logname.empty()) {
      std::string log_dir, log_instname, log_in_file_ext;
      parseFilename(log_dir, log_instname, log_in_file_ext, logname, params.logfile);
      fileWithCuts = log_dir + "/" + fileWithCuts + "_" + instname;
    }
    fileWithCuts += "_d" + stringValue(params.get(intParam::DISJ_TERMS), "%d");
    fileWithCuts += ".csv";
    cutrounds_logfile = fopen(fileWithCuts.c_str(), "w");
    if (!cutrounds_logfile) {
      error_msg(errorstring, "Unable to open file: %s.\n", fileWithCuts.c_str());
      writeErrorToLog(errorstring, params.logfile);
      exit(1);
    }
    printf("\n## Round-by-round information will be saved to file: %s. ##\n", fileWithCuts.c_str());

    int count = 0;
    fprintf(cutrounds_logfile, "%s,", instname.c_str()); count++;
    fprintf(cutrounds_logfile, "bound_gmic,"); count++;
    fprintf(cutrounds_logfile, "bound_vpc,"); count++;
    fprintf(cutrounds_logfile, "bound_all_cuts,"); count++;
    fprintf(cutrounds_logfile, "num_gmic,"); count++;
    fprintf(cutrounds_logfile, "min_density_gmic,"); count++;
    fprintf(cutrounds_logfile, "max_density_gmic,"); count++;
    fprintf(cutrounds_logfile, "avg_density_gmic,"); count++;
    fprintf(cutrounds_logfile, "gmic_gen_time,"); count++;
    fprintf(cutrounds_logfile, "gmic_apply_time,"); count++;
    //fprintf(cutrounds_logfile, "gmic_cond_num1,"); count++;
    if (COMPUTE_COND_NUM) { fprintf(cutrounds_logfile, "gmic_cond_num2,"); count++; }
    fprintf(cutrounds_logfile, "num_vpc,"); count++;
    fprintf(cutrounds_logfile, "min_density_vpc,"); count++;
    fprintf(cutrounds_logfile, "max_density_vpc,"); count++;
    fprintf(cutrounds_logfile, "avg_density_vpc,"); count++;
    fprintf(cutrounds_logfile, "vpc_gen_time,"); count++;
    fprintf(cutrounds_logfile, "vpc_apply_time,"); count++;
    //fprintf(cutrounds_logfile, "vpc_cond_num1,"); count++;
    if (COMPUTE_COND_NUM) { fprintf(cutrounds_logfile, "vpc_cond_num2,"); count++; }
    fprintf(cutrounds_logfile, "\n");
    assert(count == NUM_CUTROUND_INFO);

    //const double init_cond_num1 = compute_condition_number(solver,1);
#ifdef CALC_COND_NUM
    const double init_cond_num2 = COMPUTE_COND_NUM ? compute_condition_number(solver,2) : 0.0;
#endif

    count = 0;
    const double initSolveTime = timer.get_total_time(OverallTimeStats::INIT_SOLVE_TIME);
    fprintf(cutrounds_logfile, "%d,", 0); count++;
    fprintf(cutrounds_logfile, "%s,", stringValue(boundInfo.lp_obj,"%.20f").c_str()); count++;
    fprintf(cutrounds_logfile, "%s,", stringValue(boundInfo.lp_obj,"%.20f").c_str()); count++;
    fprintf(cutrounds_logfile, "%s,", stringValue(boundInfo.lp_obj,"%.20f").c_str()); count++;
    fprintf(cutrounds_logfile, "%d,", 0); count++;
    fprintf(cutrounds_logfile, "%f,", 0.); count++; // min density
    fprintf(cutrounds_logfile, "%f,", 0.); count++; // max density
    fprintf(cutrounds_logfile, "%f,", 0.); count++; // avg density
    fprintf(cutrounds_logfile, "%f,", 0.); count++;
    fprintf(cutrounds_logfile, "%f,", initSolveTime); count++;
    //fprintf(cutrounds_logfile, "%s,", stringValue(init_cond_num1, "%.20f").c_str()); count++;
#ifdef CALC_COND_NUM
    if (COMPUTE_COND_NUM) { fprintf(cutrounds_logfile, "%s,", stringValue(init_cond_num2, "%.20f").c_str()); count++; }
#endif
    fprintf(cutrounds_logfile, "%d,", 0); count++;
    fprintf(cutrounds_logfile, "%f,", 0.); count++; // min density
    fprintf(cutrounds_logfile, "%f,", 0.); count++; // max density
    fprintf(cutrounds_logfile, "%f,", 0.); count++; // avg density
    fprintf(cutrounds_logfile, "%f,", 0.); count++;
    fprintf(cutrounds_logfile, "%f,", initSolveTime); count++;
#ifdef CALC_COND_NUM
    //fprintf(cutrounds_logfile, "%s,", stringValue(init_cond_num1, "%.20f").c_str()); count++;
    if (COMPUTE_COND_NUM) { fprintf(cutrounds_logfile, "%s,", stringValue(init_cond_num2, "%.20f").c_str()); count++; }
#endif
    fprintf(cutrounds_logfile, "\n");
    assert(count == NUM_CUTROUND_INFO);
  } // print header to cutrounds_logfile file

  //====================================================================================================//
  // Now do rounds of cuts, until a limit is reached (e.g., time, number failures, number cuts, or all rounds are exhausted)
  std::string vpc_version_string_loop;
  std::string cbc_version_string_loop;
  std::string clp_version_string_loop;
  std::string gurobi_version_string_loop;
  std::string cplex_version_string_loop;

  MainAppState loopState {
      params,
      timer,
      OverallTimeStatsName,
      exitReason,
      boundInfo,
      info_nocuts,
      info_mycuts,
      cutInfoGMICs,
      cutInfo,
      disjInfo,
      start_time_t,
      start_time_string,
      sizeof(start_time_string),
      hostname,
      HOST_NAME_MAX,
      cpu_model,
      cpu_id,
      dir,
      filename_stub,
      instname,
      in_file_ext,
      cut_output,
      bb_output,
      solver,
      origSolver,
      GMICSolver,
      VPCSolver,
      allCutsSolver,
      vpc_version_string_loop,
      cbc_version_string_loop,
      clp_version_string_loop,
      gurobi_version_string_loop,
      cplex_version_string_loop
  };

  num_rounds = main_round_orchestrator::runCutRounds(
      num_rounds,
      COMPUTE_COND_NUM,
      cutrounds_logfile,
      NUM_CUTROUND_INFO,
      disjOptions,
      vpcs_by_round,
      cutInfoVec,
      cutInfoGMICVec,
      boundInfoVec,
      gmic_gen_time_by_round,
      gmic_apply_time_by_round,
      gmic_cond_num2_by_round,
      vpc_gen_time_by_round,
      vpc_apply_time_by_round,
      vpc_cond_num2_by_round,
      ip_solution,
      gmics,
      vpcs,
      loopState);

  printf(
      "\n## Finished VPC generation with %d VPCs and %d GMICs. ",
      boundInfo.num_vpc, boundInfo.num_gmic);
  printf(" Initial obj value: %s.",
      stringValue(boundInfo.lp_obj, "%1.6f").c_str());
  printf(" Final VPC obj value: %s.",
      stringValue(boundInfo.vpc_obj, "%1.6f").c_str());
  printf(" Final all cuts obj value: %s.",
      stringValue(boundInfo.all_cuts_obj, "%1.6f").c_str());
  printf(" Disj lb: %s.",
      stringValue(boundInfo.best_disj_obj, "%1.6f").c_str());
  printf(" ##\n");
  fflush(stdout);

  //====================================================================================================//
  // Do branch-and-bound experiments (if requested)
  if (params.get(BB_RUNS) != 0) {
    // Collect cuts from all rounds
    timer.start_timer(BB_TIME);
    runBBTests(params, &info_nocuts, &info_mycuts, &info_allcuts,
        params.get(stringParam::FILENAME), origSolver, boundInfo.ip_obj,
        &vpcs, NULL);
    timer.end_timer(BB_TIME);
  }

  timer.end_timer(OverallTimeStats::TOTAL_TIME);

#ifdef PRINT_LP_WITH_CUTS
  // If num rows in solver is > # initial rows, print cut lp file
  if (solver->getNumRows() > origSolver->getNumRows()) {
    std::string fileWithCuts = filename_stub + "_cuts";
    // solver->writeMps(fileWithCuts.c_str());
    solver->writeLp(fileWithCuts.c_str());
  }
#endif

  // Print information from each round of cuts
  if (cutrounds_logfile != NULL) {
    fclose(cutrounds_logfile);
  } // print bound by round

  //====================================================================================================//
  // Do analyses in preparation for printing
  setCutInfo(cutInfo, num_rounds, cutInfoVec.data());
  analyzeStrength(params, 
      params.get(GOMORY) == 0 ? NULL : GMICSolver,
      params.get(GOMORY) <= 0 ? solver : VPCSolver, 
      params.get(GOMORY) >= 0 ? solver : allCutsSolver, 
      cutInfoGMICs, cutInfo, 
      params.get(GOMORY) == 0 ? NULL : &gmics, 
      &vpcs,
      boundInfo, cut_output);
  analyzeBB(params, info_nocuts, info_mycuts, info_allcuts, bb_output);
  
  //====================================================================================================//
  // Finish up
  return wrapUp(0, argc, argv);
} /* main */

/**
 * Call this early to print welcome message, etc.
 */
int startUp(int argc, char** argv) {
  const std::string vpc_version_string;
  const std::string cbc_version_string;
  const std::string clp_version_string;
  const std::string gurobi_version_string;
  const std::string cplex_version_string;

  MainAppState state {
      params,
      timer,
      OverallTimeStatsName,
      exitReason,
      boundInfo,
      info_nocuts,
      info_mycuts,
      cutInfoGMICs,
      cutInfo,
      disjInfo,
      start_time_t,
      start_time_string,
      sizeof(start_time_string),
      hostname,
      HOST_NAME_MAX,
      cpu_model,
      cpu_id,
      dir,
      filename_stub,
      instname,
      in_file_ext,
      cut_output,
      bb_output,
      solver,
      origSolver,
      GMICSolver,
      VPCSolver,
      allCutsSolver,
      vpc_version_string,
      cbc_version_string,
      clp_version_string,
      gurobi_version_string,
      cplex_version_string
  };

  return main_startup::startUp(
      argc,
      argv,
      state);
} /* startUp */

/**
 * Close the logfile and print to it
 */
int wrapUp(int retCode, int argc, char** argv) {
  std::string vpc_version_string;
  std::string cbc_version_string;
  std::string clp_version_string;
  std::string gurobi_version_string;
  std::string cplex_version_string;
#ifdef VPC_VERSION
  vpc_version_string = VPC_VERSION_STRING;
#endif
#ifdef VPC_CBC_VERSION
  cbc_version_string = CBC_VERSION_STRING;
#endif
#ifdef VPC_CLP_VERSION
  clp_version_string = CLP_VERSION_STRING;
#endif
#ifdef USE_GUROBI
  gurobi_version_string = GUROBI_VERSION_STRING;
#endif
#ifdef USE_CPLEX
  cplex_version_string = CPLEX_VERSION_STRING;
#endif

  MainAppState state {
      params,
      timer,
      OverallTimeStatsName,
      exitReason,
      boundInfo,
      info_nocuts,
      info_mycuts,
      cutInfoGMICs,
      cutInfo,
      disjInfo,
      start_time_t,
      start_time_string,
      sizeof(start_time_string),
      hostname,
      HOST_NAME_MAX,
      cpu_model,
      cpu_id,
      dir,
      filename_stub,
      instname,
      in_file_ext,
      cut_output,
      bb_output,
      solver,
      origSolver,
      GMICSolver,
      VPCSolver,
      allCutsSolver,
      vpc_version_string,
      cbc_version_string,
      clp_version_string,
      gurobi_version_string,
      cplex_version_string
  };

  return main_reporting::wrapUp(
      retCode, argc, argv, state);
} /* wrapUp */

void doCustomRoundOfCuts(int round_ind, OsiCuts& vpcs, CglVPC& gen, int& num_disj) {
  const std::string vpc_version_string;
  const std::string cbc_version_string;
  const std::string clp_version_string;
  const std::string gurobi_version_string;
  const std::string cplex_version_string;

  MainAppState state {
      params,
      timer,
      OverallTimeStatsName,
      exitReason,
      boundInfo,
      info_nocuts,
      info_mycuts,
      cutInfoGMICs,
      cutInfo,
      disjInfo,
      start_time_t,
      start_time_string,
      sizeof(start_time_string),
      hostname,
      HOST_NAME_MAX,
      cpu_model,
      cpu_id,
      dir,
      filename_stub,
      instname,
      in_file_ext,
      cut_output,
      bb_output,
      solver,
      origSolver,
      GMICSolver,
      VPCSolver,
      allCutsSolver,
      vpc_version_string,
      cbc_version_string,
      clp_version_string,
      gurobi_version_string,
      cplex_version_string
  };

  main_round_actions::doCustomRoundOfCuts(
      round_ind,
      vpcs,
      gen,
      num_disj,
      state,
      cutInfoVec);
} /* doCustomRoundOfCuts */

/**
 * See params.hpp for descriptions of the parameters
 */
int processArgs(int argc, char** argv) {
  return main_input::processArgs(argc, argv, params);
} /* processArgs */
