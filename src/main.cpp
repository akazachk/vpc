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
  boundInfo.num_vpc = 0, boundInfo.num_gmic = 0;
  int round_ind = 0;
  for (round_ind = 0; round_ind < num_rounds; ++round_ind) {
    if (num_rounds > 1) {
      printf("\n## Starting round %d/%d. ##\n", round_ind+1, num_rounds);
    }
    int num_disj = 0;

    // Initialize all the boundInfo entries for this round
    boundInfoVec[round_ind].lp_obj = solver->getObjValue();
    boundInfoVec[round_ind].num_vpc = 0;
    boundInfoVec[round_ind].num_gmic = 0;

    timer.start_timer(OverallTimeStats::CUT_TOTAL_TIME);

    //====================================================================================================//
    // Generate Gomory cuts
    // > 0: apply cuts to solver (affects disjunctive cuts generated);
    // < 0: apply cuts to strCutSolver only
    // Option 1: GglGMI
    // (to be added) Option 2: custom generate intersection cuts, calculate Farkas certificate, do strengthening
    // (to be added) Option 3: custom generate intersection cuts, calculate Farkas certificate, do closed-form strengthening
    OsiCuts currGMICs;
    if (params.get(GOMORY) == -1 || params.get(GOMORY) == 1) {
      // Generate GMICs
      const double gmic_gen_start_time = timer.get_total_time(OverallTimeStats::GOMORY_GEN_TIME);
      timer.start_timer(OverallTimeStats::GOMORY_GEN_TIME);

      CglGMI GMIGen;
      // Set parameters so that many GMIs are generated
      // CglGMI's MAX_SUPPORT parameter is equivalent to the MIN_SUPPORT_THRESHOLD value in this code
      // CglGMI's sparsity parameter rejects with bad support when # cut coefficient nonzeros >= MAX_SUPPORT + MAX_SUPPORT_REL * # cols
      // This means that there will be a constant offset on the max allowable support for GMICs compared to VPCs (in the favor of GMICs being denser)
      GMIGen.getParam().setMAX_SUPPORT(params.get(intConst::MIN_SUPPORT_THRESHOLD));
      GMIGen.getParam().setMAX_SUPPORT_REL(params.get(doubleParam::MAX_SUPPORT_REL));
      //GMIGen.getParam().setMAXDYN(params.get(doubleConst::MAX_DYN));
      GMIGen.getParam().setMAXDYN(solver->getInfinity());
      GMIGen.getParam().setMINVIOL(0.0);
      GMIGen.generateCuts(*solver, currGMICs);

      timer.end_timer(OverallTimeStats::GOMORY_GEN_TIME);
      const double gmic_gen_end_time = timer.get_total_time(OverallTimeStats::GOMORY_GEN_TIME);
      gmic_gen_time_by_round[round_ind] = gmic_gen_end_time - gmic_gen_start_time;

      // Apply GMICs to solvers
      boundInfo.num_gmic += currGMICs.sizeCuts();
      boundInfoVec[round_ind].num_gmic += currGMICs.sizeCuts();
      gmics.insert(currGMICs);

      // For timing, only count first application, to solver with only GMICs,
      // to get a sense of how a single LP relaxation is affected
      timer.start_timer(OverallTimeStats::GOMORY_APPLY_TIME);
      const double gmic_apply_start_time = timer.get_total_time(OverallTimeStats::GOMORY_APPLY_TIME);
      applyCutsCustom(GMICSolver, currGMICs, params.logfile);
      const double gmic_apply_end_time = timer.get_total_time(OverallTimeStats::GOMORY_APPLY_TIME);
      timer.end_timer(OverallTimeStats::GOMORY_APPLY_TIME);
      gmic_apply_time_by_round[round_ind] = gmic_apply_end_time - gmic_apply_start_time;

      // Update bound
      boundInfo.gmic_obj = GMICSolver->getObjValue();
      boundInfoVec[round_ind].gmic_obj = boundInfo.gmic_obj;

      // Update condition numbers
#ifdef CALC_COND_NUM
      //gmic_cond_num1_by_round[round_ind] = compute_condition_number(GMICSolver, 1);
      if (COMPUTE_COND_NUM) {
        gmic_cond_num2_by_round[round_ind] = compute_condition_number(GMICSolver, 2);
      }
#endif

      if (params.get(GOMORY) > 0) {
        applyCutsCustom(solver, currGMICs, params.logfile);
      }
      applyCutsCustom(allCutsSolver, currGMICs, params.logfile);

      updateGMICInfo(cutInfoGMICVec[round_ind], &currGMICs, params.get(EPS) / 2.);
    } // Gomory cut generation

    //====================================================================================================//
    // Now for more general cuts
    if ((params.get(VPCParametersNamespace::intParam::MODE) != static_cast<int>(VPCParametersNamespace::VPCMode::DISJ_SET_PBB))
        && ((int) disjOptions.size() > round_ind)) {
      params.set(intParam::DISJ_TERMS, disjOptions[round_ind]);
    }
    const bool SHOULD_GENERATE_VPCS = (params.get(intParam::DISJ_TERMS) != 0 || disjOptions.size() > 0);
    const bool USE_CUSTOM = (params.get(intParam::MODE) == static_cast<int>(VPCParametersNamespace::VPCMode::CUSTOM));

    if (SHOULD_GENERATE_VPCS) {
      const double vpc_start_time = timer.get_total_time(OverallTimeStats::VPC_GEN_TIME);
      timer.start_timer(OverallTimeStats::VPC_GEN_TIME);

      CglVPC gen(params, round_ind);

      // Store the initial solve time in order to set a baseline for the PRLP resolve time
      gen.timer.add_value(
          CglVPC::VPCTimeStatsName[static_cast<int>(CglVPC::VPCTimeStats::INIT_SOLVE_TIME)],
          timer.get_value(OverallTimeStats::INIT_SOLVE_TIME));

      // Proceed with custom disjunctions if specified; otherwise, the disjunction will be set up in the CglVPC class
      if (USE_CUSTOM) {
        doCustomRoundOfCuts(round_ind, vpcs_by_round[round_ind], gen, num_disj);
      } else {
        // Generate disjunctive cuts
        gen.generateCuts(*solver, vpcs_by_round[round_ind]); // solution may change slightly due to enable factorization called in getProblemData...

        // Update statistics about the disjunction objective value and cuts
        exitReason = gen.exitReason;
        if (gen.disj() || gen.disjSet()) {
          DisjunctionSet* disjSet = gen.disjSet();
          if (gen.disjSet() == NULL) {
            disjSet = new DisjunctionSet;
            disjSet->addDisjunction(gen.disj());
          }

          boundInfo.num_vpc += gen.num_cuts; // TODO: does this need to be in this "if", and do we need to subtract initial num cuts?
          boundInfoVec[round_ind].num_vpc += gen.num_cuts;

          for (const Disjunction* const disj : disjSet->disjunctions) {
            num_disj++;

            if (boundInfo.best_disj_obj < disj->best_obj) {
              boundInfo.best_disj_obj = disj->best_obj;
            }
            if (boundInfoVec[round_ind].best_disj_obj < disj->best_obj) {
              boundInfoVec[round_ind].best_disj_obj = disj->best_obj;
            }
            if (boundInfo.worst_disj_obj < disj->worst_obj) {
              boundInfo.worst_disj_obj = disj->worst_obj;
            }
            if (boundInfoVec[round_ind].worst_disj_obj < disj->worst_obj) {
              boundInfoVec[round_ind].worst_disj_obj = disj->worst_obj;
            }

            if (round_ind == 0) {
              // The following should be the same across all disjunctions
              boundInfo.num_root_bounds_changed = disj->common_changed_var.size();
              boundInfo.root_obj = disj->root_obj;
            }
            boundInfoVec[round_ind].num_root_bounds_changed = disj->common_changed_var.size();
            boundInfoVec[round_ind].root_obj = disj->root_obj;
          } // loop over disjunctions used for this round

          if (gen.disjSet() == NULL) {
            delete disjSet;
            disjSet = NULL;
          }
        } // check if any disjunctions were used
        updateDisjInfo(disjInfo, num_disj, gen);
        updateCutInfo(cutInfoVec[round_ind], gen, &vpcs_by_round[round_ind], params.get(EPS) / 2.);
      } // check if mode is _not_ CUSTOM

      // Update timing from underlying generator
      const TimeStats& vpctimer = gen.timer;
      std::vector<CglVPC::VPCTimeStats> vpc_stats = {
        CglVPC::VPCTimeStats::INIT_SOLVE_TIME,
        CglVPC::VPCTimeStats::DISJ_SETUP_TIME,
        CglVPC::VPCTimeStats::DISJ_GEN_TIME,
        CglVPC::VPCTimeStats::PRLP_SETUP_TIME,
        CglVPC::VPCTimeStats::PRLP_SOLVE_TIME,
        CglVPC::VPCTimeStats::GEN_CUTS_TIME,
      };
      std::vector<OverallTimeStats> overall_stats = {
        OverallTimeStats::VPC_INIT_SOLVE_TIME,
        OverallTimeStats::VPC_DISJ_SETUP_TIME,
        OverallTimeStats::VPC_DISJ_GEN_TIME,
        OverallTimeStats::VPC_PRLP_SETUP_TIME,
        OverallTimeStats::VPC_PRLP_SOLVE_TIME,
        OverallTimeStats::VPC_GEN_CUTS_TIME
      };
      for (int i = 0; i < (int) vpc_stats.size(); i++) {
        const CglVPC::VPCTimeStats stat = vpc_stats[i];
        const OverallTimeStats overall_stat = overall_stats[i];
        const clock_t currtimevalue = vpctimer.get_value(CglVPC::VPCTimeStatsName[static_cast<int>(stat)]);
        timer.add_value(overall_stat, currtimevalue);
      }

      timer.end_timer(OverallTimeStats::VPC_GEN_TIME);
      const double vpc_end_time = timer.get_total_time(OverallTimeStats::VPC_GEN_TIME);
      vpc_gen_time_by_round[round_ind] = vpc_end_time - vpc_start_time;
    } // check if SHOULD_GENERATE_VPCS (num disj terms requested != 0)
    else { // else update cutInfo with blanks
      setCutInfo(cutInfoVec[round_ind], 0, NULL);
    }

    //====================================================================================================//
    // Check cuts against IP solution
    if (params.get(TEMP) == static_cast<int>(TempOptions::CHECK_CUTS_AGAINST_BB_OPT) && !ip_solution.empty()) {
      const int num_violated = checkCutsAgainstFeasibleSolution(vpcs_by_round[round_ind], ip_solution);
      printf("\n## Number of cuts violating the IP solution: %d ##\n", num_violated);
    }

    // Insert cuts generated in this round into set of all cuts
    vpcs.insert(vpcs_by_round[round_ind]);

    timer.end_timer(OverallTimeStats::CUT_TOTAL_TIME);
    
    //====================================================================================================//
    // Apply cuts
    printf("\n## Applying disjunctive cuts (# cuts = %d). ##\n", (int) vpcs_by_round[round_ind].sizeCuts());
    timer.start_timer(OverallTimeStats::TOTAL_APPLY_TIME);

    // To track time to apply VPCs, need to be careful of which solver they are being added to
    // If GOMORY = 0, then no GMICs generated; solver is only one that exists
    //   (VPCSolver and allCutsSolver are *not* created; but if they were, it would hold that solver = VPCSolver = allCutsSolver)
    // If GOMORY = 1, then solver has GMICs from this round applied already (will be identical to allCutsSolver)
    // If GOMORY = -1, then solver has only VPCs from prior rounds (will be identical to VPCSolver)
    // Recall that "solver" is the one used for generating cuts for each round
    if (vpcs_by_round[round_ind].sizeCuts() > 0) {
      if (params.get(GOMORY) != 0) { // GMICs generated, so need to track VPC-only objective with VPCSolver
        // Apply to "VPCSolver", which tracks effect of VPCs without other cuts
        const double vpc_apply_start_time = timer.get_total_time(OverallTimeStats::VPC_APPLY_TIME);
        timer.start_timer(OverallTimeStats::VPC_APPLY_TIME);
        applyCutsCustom(VPCSolver, vpcs_by_round[round_ind], params.logfile);
        timer.end_timer(OverallTimeStats::VPC_APPLY_TIME);
        const double vpc_apply_end_time = timer.get_total_time(OverallTimeStats::VPC_APPLY_TIME);
        vpc_apply_time_by_round[round_ind] = vpc_apply_end_time - vpc_apply_start_time;

        // Apply to "solver" which has either only VPCs or also GMICs depending on GOMORY parameter
        // When GOMORY = -1, VPCSolver and solver are identical
        applyCutsCustom(solver, vpcs_by_round[round_ind], params.logfile);

        // Apply to "allCutsSolver", which has all cuts
        // When GOMORY = 1, allCutsSolver and solver are identical
        applyCutsCustom(allCutsSolver, vpcs_by_round[round_ind], params.logfile);
      } // check if GMICs generated (params.get(GOMORY) != 0)
      else {
        // Else, no GMICs generated; only solver needed
        const double vpc_apply_start_time = timer.get_total_time(OverallTimeStats::VPC_APPLY_TIME);
        timer.start_timer(OverallTimeStats::VPC_APPLY_TIME);
        applyCutsCustom(solver, vpcs_by_round[round_ind], params.logfile);
        timer.end_timer(OverallTimeStats::VPC_APPLY_TIME);
        const double vpc_apply_end_time = timer.get_total_time(OverallTimeStats::VPC_APPLY_TIME);
        vpc_apply_time_by_round[round_ind] = vpc_apply_end_time - vpc_apply_start_time;
      }
    } // apply VPCs if any were generated
    timer.end_timer(OverallTimeStats::TOTAL_APPLY_TIME);

    // Update bound info
    boundInfo.vpc_obj      = (params.get(GOMORY) != 0) ? VPCSolver->getObjValue() : solver->getObjValue();
    boundInfo.gmic_vpc_obj = (params.get(GOMORY) != 0) ? allCutsSolver->getObjValue() : solver->getObjValue();
    boundInfo.all_cuts_obj = boundInfo.gmic_vpc_obj;

    // Also update the bound info for this round
    boundInfoVec[round_ind].gmic_vpc_obj = boundInfo.gmic_vpc_obj;
    boundInfoVec[round_ind].vpc_obj      = boundInfo.vpc_obj;
    boundInfoVec[round_ind].all_cuts_obj = boundInfo.all_cuts_obj;

    // Update condition number for VPCs
#ifdef CALC_COND_NUM
    //vpc_cond_num1_by_round[round_ind] = compute_condition_number((params.get(GOMORY) != 0) ? VPCSolver : solver, 1);
    if (COMPUTE_COND_NUM) {
      vpc_cond_num2_by_round[round_ind] = compute_condition_number((params.get(GOMORY) != 0) ? VPCSolver : solver, 2);
    }
#endif

    printf("\n*** INFO:\n");
    printf(
        "Round %d/%d: Completed round of VPC generation (exit reason: %s).\n", 
        round_ind + 1, params.get(ROUNDS),
        CglVPC::ExitReasonName[static_cast<int>(exitReason)].c_str());
    printf(
        "VPCs generated = %d. GMICs generated = %d.\n",
        vpcs_by_round[round_ind].sizeCuts(),
        boundInfoVec[round_ind].num_gmic);
    printf("Obj value: %s.\n",
        stringValue(solver->getObjValue(), "%1.6f").c_str());
    printf("Initial obj value: %s.\n",
        stringValue(boundInfo.lp_obj, "%1.6f").c_str());
    if (round_ind > 0) {
      printf("Round start obj value: %s.\n",
          stringValue(boundInfoVec[round_ind-1].all_cuts_obj, "%1.6f").c_str());
    }
    printf("Disj lb: %s.\n",
        stringValue(boundInfo.best_disj_obj, "%1.6f").c_str());
    printf("Elapsed time: %1.2f seconds. Time limit = %1.2f seconds.\n",
        timer.get_total_time(OverallTimeStatsName[OverallTimeStats::TOTAL_TIME]),
        params.get(TIMELIMIT));
    printf("***\n");
    fflush(stdout);

    // Print information from this round of cuts
    if (use_temp_option(params.get(intParam::TEMP), TempOptions::PRINT_BOUND_BY_ROUND)) {
      int count = 0;
      fprintf(cutrounds_logfile, "%d,", round_ind+1); count++;
      fprintf(cutrounds_logfile, "%s,", stringValue(boundInfoVec[round_ind].gmic_obj,"%.20f").c_str()); count++;
      fprintf(cutrounds_logfile, "%s,", stringValue(boundInfoVec[round_ind].vpc_obj,"%.20f").c_str()); count++;
      fprintf(cutrounds_logfile, "%s,", stringValue(boundInfoVec[round_ind].all_cuts_obj,"%.20f").c_str()); count++;
      fprintf(cutrounds_logfile, "%d,", boundInfoVec[round_ind].num_gmic); count++;
      fprintf(cutrounds_logfile, "%s,", 
        stringValue(
          (boundInfoVec[round_ind].num_gmic > 0)
            ? 100. * cutInfoGMICVec[round_ind].min_support / solver->getNumCols()
            : 0.,
          "%.6f",
          1e100
        ).c_str()
      ); count++;
      fprintf(cutrounds_logfile, "%s,", 
        stringValue(
          (boundInfoVec[round_ind].num_gmic > 0)
            ? 100. * cutInfoGMICVec[round_ind].max_support / solver->getNumCols()
            : 0.,
          "%.6f",
          1e100
        ).c_str()
      ); count++;
      fprintf(cutrounds_logfile, "%s,", 
        stringValue(
          (boundInfoVec[round_ind].num_gmic > 0)
            ? 100. * cutInfoGMICVec[round_ind].avg_support / solver->getNumCols()
            : 0.,
          "%.6f",
          1e100
        ).c_str()
      ); count++;
      fprintf(cutrounds_logfile, "%.20f,", gmic_gen_time_by_round[round_ind]); count++;
      fprintf(cutrounds_logfile, "%.20f,", gmic_apply_time_by_round[round_ind]); count++;
#ifdef CALC_COND_NUM
      //fprintf(cutrounds_logfile, "%s,", stringValue(gmic_cond_num1_by_round[round_ind], "%.20f").c_str()); count++;
      if (COMPUTE_COND_NUM) { fprintf(cutrounds_logfile, "%s,", stringValue(gmic_cond_num2_by_round[round_ind], "%.20f").c_str()); count++; }
#endif
      fprintf(cutrounds_logfile, "%d,", boundInfoVec[round_ind].num_vpc); count++;
      fprintf(cutrounds_logfile, "%s,", 
        stringValue(
          (boundInfoVec[round_ind].num_vpc > 0)
            ? 100. * cutInfoVec[round_ind].min_support / solver->getNumCols()
            : 0.,
          "%.6f",
          1e100
        ).c_str()
      ); count++;
      fprintf(cutrounds_logfile, "%s,", 
        stringValue(
          (boundInfoVec[round_ind].num_vpc > 0)
            ? 100. * cutInfoVec[round_ind].max_support / solver->getNumCols()
            : 0.,
          "%.6f",
          1e100
        ).c_str()
      ); count++;
      fprintf(cutrounds_logfile, "%s,", 
        stringValue(
          (boundInfoVec[round_ind].num_vpc > 0)
            ? 100. * cutInfoVec[round_ind].avg_support / solver->getNumCols()
            : 0.,
          "%.6f",
          1e100
        ).c_str()
      ); count++;
      fprintf(cutrounds_logfile, "%.20f,", vpc_gen_time_by_round[round_ind]); count++;
      fprintf(cutrounds_logfile, "%.20f,", vpc_apply_time_by_round[round_ind]); count++;
#ifdef CALC_COND_NUM
      //fprintf(cutrounds_logfile, "%s,", stringValue(vpc_cond_num1_by_round[round_ind], "%.20f").c_str()); count++;
      if (COMPUTE_COND_NUM) { fprintf(cutrounds_logfile, "%s,", stringValue(vpc_cond_num2_by_round[round_ind], "%.20f").c_str()); count++; }
#endif
      fprintf(cutrounds_logfile, "\n");
      assert(count == NUM_CUTROUND_INFO);
    } // print info from this round

    // Exit early if reached time limit
    if (timer.reachedTimeLimit(OverallTimeStatsName[OverallTimeStats::TOTAL_TIME], params.get(TIMELIMIT))) {
      printf("\n*** INFO: Reached time limit (current time = %f, time limit = %f).\n",
          timer.get_total_time(OverallTimeStatsName[OverallTimeStats::TOTAL_TIME]), params.get(TIMELIMIT));
      break;
    }

    // Exit early from rounds of cuts if no cuts generated or solver is not optimal
    if ((SHOULD_GENERATE_VPCS && vpcs_by_round[round_ind].sizeCuts() == 0)
        || (!SHOULD_GENERATE_VPCS && boundInfoVec[round_ind].num_gmic == 0)
        || !solver->isProvenOptimal()
        || isInfinity(std::abs(solver->getObjValue()))
        || exitReason == CglVPC::ExitReason::OPTIMAL_SOLUTION_FOUND_EXIT) {
      break;
    }
  } // loop over rounds of cuts
  if (round_ind < num_rounds) {
    num_rounds = round_ind+1;
  }

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
  return main_startup::startUp(
      argc,
      argv,
      params,
      exitReason,
      boundInfo,
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
      OverallTimeStatsName);
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

  return main_reporting::wrapUp(
      retCode, argc, argv,
      params, timer, OverallTimeStatsName, exitReason,
      start_time_t, start_time_string,
      hostname, cpu_model, cpu_id,
      instname, cut_output, bb_output,
      boundInfo, info_nocuts, info_mycuts,
      cutInfoGMICs, cutInfo, disjInfo,
      solver, origSolver, GMICSolver, VPCSolver, allCutsSolver,
      vpc_version_string, cbc_version_string, clp_version_string,
      gurobi_version_string, cplex_version_string);
} /* wrapUp */

void doCustomRoundOfCuts(int round_ind, OsiCuts& vpcs, CglVPC& gen, int& num_disj) {
  printf("\n## Setting up disjunction(s) ##\n");

  gen.timer.start_timer(CglVPC::VPCTimeStatsName[static_cast<int>(CglVPC::VPCTimeStats::DISJ_GEN_TIME)]);
  DisjunctionSet* disjSet = new DisjunctionSet();
  CglVPC::ExitReason setDisjExitReason = setDisjunctions(disjSet, solver, params, VPCParametersNamespace::VPCMode::SPLITS);
  gen.timer.end_timer(CglVPC::VPCTimeStatsName[static_cast<int>(CglVPC::VPCTimeStats::DISJ_GEN_TIME)]);
  
  const int numDisj = disjSet->size();

  // If integer-optimal solution was found, all disjunctions but one will have been deleted
  if (setDisjExitReason == CglVPC::ExitReason::OPTIMAL_SOLUTION_FOUND_EXIT) {
    warning_msg(warnstr,
        "An integer (optimal) solution was found prior while getting disjunction. "
        "We will generate between n and 2n cuts, restricting the value of each variable.\n");
    const double* solution = disjSet->integer_sol.data();
    if (!solution) {
      error_msg(errorstring,
          "Though status is that optimal integer solution found, unable to find this solution.\n");
      writeErrorToLog(errorstring, params.logfile);
      exit(1);
    }
    for (int col = 0; col < solver->getNumCols(); col++) {
      const double val = solution[col];

      // Check which of the bounds needs to be fixed
      for (int b = 0; b < 2; b++) {
        if ((b == 0 && greaterThanVal(val, solver->getColLower()[col]))
            || (b == 1 && lessThanVal(val, solver->getColUpper()[col]))) {
          const double mult = (b == 0) ? 1. : -1.;
          const double el = mult * 1.;

          OsiRowCut currCut;
          currCut.setLb(mult * val);
          currCut.setRow(1, &col, &el, false);
          gen.addCut(currCut, vpcs,
              CglVPC::CutType::OPTIMALITY_CUT,
              CglVPC::ObjectiveType::ONE_SIDED,
              0);
        }
      }
    } // iterate over columns and add optimality cut if needed
    exitReason = CglVPC::ExitReason::OPTIMAL_SOLUTION_FOUND_EXIT;
    boundInfo.num_vpc += gen.num_cuts;
    updateCutInfo(cutInfoVec[round_ind], gen, &vpcs, params.get(EPS) / 2.);
  } // check if integer-optimal solution
  else if (setDisjExitReason == CglVPC::ExitReason::SUCCESS_EXIT && numDisj > 0) {
    const int cutLimit = std::ceil(
        gen.getCutLimit(params.get(CUTLIMIT),
            solver->getFractionalIndices().size()) / numDisj); // distribute cut limit over the disjunctions
    gen.setupRepeatedUse(true);
    for (Disjunction* disj : disjSet->disjunctions) {
      if (!disj)
        continue;
#ifdef TRACE
      printf("\n## Generating cuts from disj %s ##\n", disj->name.c_str());
#endif
      num_disj++;
      gen.params.set(CUTLIMIT, cutLimit);
      gen.setDisjunction(disj, false);
      gen.generateCuts(*solver, vpcs); // solution may change slightly due to enable factorization called in getProblemData...
      exitReason = gen.exitReason;
      boundInfo.num_vpc += gen.num_cuts;
      if (boundInfo.best_disj_obj < gen.disj()->best_obj)
        boundInfo.best_disj_obj = gen.disj()->best_obj;
      if (boundInfo.worst_disj_obj < gen.disj()->worst_obj)
        boundInfo.worst_disj_obj = gen.disj()->worst_obj;
      updateDisjInfo(disjInfo, num_disj, gen);
      updateCutInfo(cutInfoVec[round_ind], gen, &vpcs, params.get(EPS) / 2.);
    }
    gen.setupRepeatedUse(false);
  }  // if successful generation of disjunction, generate cuts
  else if (setDisjExitReason == CglVPC::ExitReason::NO_DISJUNCTION_EXIT) {
    // Do nothing
    exitReason = CglVPC::ExitReason::NO_DISJUNCTION_EXIT;
  }
  else {
    error_msg(errorstring, "Unknown exit reason (%s) from setDisjunctions.\n", CglVPC::ExitReasonName[static_cast<int>(setDisjExitReason)].c_str());
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }

  if (disjSet) {
    delete disjSet;
  }
} /* doCustomRoundOfCuts */

/**
 * See params.hpp for descriptions of the parameters
 */
int processArgs(int argc, char** argv) {
  return main_input::processArgs(argc, argv, params);
} /* processArgs */
