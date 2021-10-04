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
#include "preprocess.hpp" // performCleaning
#include "SolverHelper.hpp" // initializeSolver
#include "SolverInterface.hpp"
#include "VPCParameters.hpp"
using namespace VPCParametersNamespace;
#include "TimeStats.hpp"
#include "utility.hpp"

enum OverallTimeStats {
  INIT_SOLVE_TIME,
  VPC_GEN_TIME,
  VPC_APPLY_TIME,
  BB_TIME,
  TOTAL_TIME,
  NUM_TIME_STATS
}; /* OverallTimeStats */
const std::vector<std::string> OverallTimeStatsName {
  "INIT_SOLVE_TIME",
  "VPC_GEN_TIME",
  "VPC_APPLY_TIME",
  "BB_TIME",
  "TOTAL_TIME"
};

// Main file variables
VPCParameters params;
OsiSolverInterface *solver, *origSolver;
OsiSolverInterface* GMICSolver = NULL;
OsiSolverInterface* VPCSolver = NULL;
OsiCuts gmics, vpcs;
std::string dir = "", filename_stub = "", instname = "", in_file_ext = "";
CglVPC::ExitReason exitReason;
TimeStats timer;
std::time_t start_time_t, end_time_t;
char start_time_string[25];

SummaryBoundInfo boundInfo;
SummaryBBInfo info_nocuts, info_mycuts, info_allcuts;
SummaryDisjunctionInfo disjInfo;
std::vector<SummaryCutInfo> cutInfoVec;
SummaryCutInfo cutInfo, cutInfoGMICs;

// For output
std::string cut_output = "", bb_output = "";

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
  error_msg(errorstring, "Abort or seg fault message received. Signal number: %d.\n", signal_number);
  writeErrorToLog(errorstring, params.logfile);
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
  initializeSolver(solver, params.get(stringParam::FILENAME));
  timer.start_timer(OverallTimeStats::INIT_SOLVE_TIME);
  solver->initialSolve();
  if (!checkSolverOptimality(solver, false)) {
    error_msg(errorstring, "Unable to solve initial LP relaxation.\n");
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }
  timer.end_timer(OverallTimeStats::INIT_SOLVE_TIME);
  boundInfo.lp_obj = solver->getObjValue();

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
  // Save original solver in case we wish to come back to it later
  origSolver = solver->clone();
  if (!origSolver->isProvenOptimal()) {
    origSolver->initialSolve();
    checkSolverOptimality(origSolver, false);
  }

  // Also save copies for calculating other objective values
  // We only need these if cuts other than VPCs are generated
  // solver      ::  stores the cuts we actually want to "count"
  // GMICSolver  ::  only GMICs
  // VPCSolver   ::  if GMICs count, this is only VPCs; otherwise, it is both GMICs and VPCs
  if (params.get(GOMORY) != 0) {
    GMICSolver = solver->clone();
    VPCSolver = solver->clone();
  }

  // Possibly preprocess instead of doing cuts
  if (params.get(PREPROCESS) != 0) {
    // Cleaning involves running presolve and branching
    printf("\n## Starting preprocessing/cleaning of instance. ##\n");
    params.set(intParam::BB_MODE, 1); // only do no cuts branching
    performCleaning(params, solver, filename_stub, boundInfo.ip_obj, params.get(PREPROCESS));

    printf("\n## Finished preprocessing/cleaning of instance. ##\n");
    return wrapUp(0, argc, argv);
  }

  // Now do rounds of cuts, until a limit is reached (e.g., time, number failures, number cuts, or all rounds are exhausted)
  boundInfo.num_vpc = 0, boundInfo.num_gmic = 0;
  int num_rounds = params.get(ROUNDS);
  std::vector<OsiCuts> vpcs_by_round(num_rounds);
  cutInfoVec.resize(num_rounds);
  int round_ind = 0;
  int num_disj = 0;
  std::vector<int> disjOptions;
  if (!params.get(stringParam::DISJ_OPTIONS).empty()) {
    // TODO
    error_msg(errorstring, "Setting DISJ_OPTIONS is currently not implemented. Received from command line: %s.\n", params.get(stringParam::DISJ_OPTIONS).c_str());
    writeErrorToLog(errorstring, params.logfile);
    return wrapUp(1, argc, argv);
  }
  for (round_ind = 0; round_ind < num_rounds; ++round_ind) {
    if (num_rounds > 1) {
      printf("\n## Starting round %d/%d. ##\n", round_ind+1, num_rounds);
    }
    timer.start_timer(OverallTimeStats::VPC_GEN_TIME);

    if (params.get(GOMORY) == -1 || params.get(GOMORY) == 1) {
      OsiCuts currGMICs;
      CglGMI GMIGen;
      // Set parameters so that many GMIs are generated
      GMIGen.getParam().setMAX_SUPPORT(solver->getNumCols());
      GMIGen.getParam().setMAX_SUPPORT_REL(0.5);
      //GMIGen.getParam().setMAXDYN(params.get(doubleConst::MAX_DYN));
      GMIGen.getParam().setMAXDYN(solver->getInfinity());
      GMIGen.generateCuts(*solver, currGMICs);
      gmics.insert(currGMICs);
      boundInfo.num_gmic += currGMICs.sizeCuts();
      applyCutsCustom(GMICSolver, currGMICs, params.logfile);
      boundInfo.gmic_obj = GMICSolver->getObjValue();
      if (params.get(GOMORY) > 0) {
        applyCutsCustom(solver, currGMICs, params.logfile);
      }
      if (params.get(GOMORY) < 0) {
        applyCutsCustom(VPCSolver, currGMICs, params.logfile);
      }
    }

    if (disjOptions.size() > round_ind) {
      params.set(intParam::DISJ_TERMS, disjOptions[round_ind]);
    }
    CglVPC gen(params, round_ind);

    // Store the initial solve time in order to set a baseline for the PRLP resolve time
    gen.timer.add_value(
        CglVPC::VPCTimeStatsName[static_cast<int>(CglVPC::VPCTimeStats::INIT_SOLVE_TIME)],
        timer.get_value(OverallTimeStats::INIT_SOLVE_TIME));

    // Proceed with custom disjunctions if specified; otherwise, the disjunction will be set up in the CglVPC class
    if (params.get(MODE) != static_cast<int>(CglVPC::VPCMode::CUSTOM)) {
      gen.generateCuts(*solver, vpcs_by_round[round_ind]); // solution may change slightly due to enable factorization called in getProblemData...
      exitReason = gen.exitReason;
      if (gen.disj()) {
        num_disj++;
        boundInfo.num_vpc += gen.num_cuts; // TODO: does this need to be in this if, and do we need to subtract initial num cuts?
        if (boundInfo.best_disj_obj < gen.disj()->best_obj)
          boundInfo.best_disj_obj = gen.disj()->best_obj;
        if (boundInfo.worst_disj_obj < gen.disj()->worst_obj)
          boundInfo.worst_disj_obj = gen.disj()->worst_obj;
      }
      updateDisjInfo(disjInfo, num_disj, gen);
      updateCutInfo(cutInfoVec[round_ind], gen);
    } // check if mode is _not_ CUSTOM
    else {
      doCustomRoundOfCuts(round_ind, vpcs_by_round[round_ind], gen, num_disj);
    } // check if mode is CUSTOM
    timer.end_timer(OverallTimeStats::VPC_GEN_TIME);

    timer.start_timer(OverallTimeStats::VPC_APPLY_TIME);
    applyCutsCustom(solver, vpcs_by_round[round_ind], params.logfile);
    if (params.get(GOMORY) > 0) { // GMICs added to solver, so VPCSolver tracks values without GMICs
      boundInfo.gmic_vpc_obj = solver->getObjValue();
      applyCutsCustom(VPCSolver, vpcs_by_round[round_ind], params.logfile);
      boundInfo.vpc_obj = VPCSolver->getObjValue();
      boundInfo.all_cuts_obj = boundInfo.gmic_vpc_obj;
    }
    else if (params.get(GOMORY) < 0) {
      boundInfo.vpc_obj = solver->getObjValue();
      applyCutsCustom(VPCSolver, vpcs_by_round[round_ind], params.logfile);
      boundInfo.gmic_vpc_obj = VPCSolver->getObjValue();
      boundInfo.all_cuts_obj = boundInfo.gmic_vpc_obj;
    } else {
      boundInfo.vpc_obj = solver->getObjValue();
      boundInfo.all_cuts_obj = boundInfo.vpc_obj;
    }
    timer.end_timer(OverallTimeStats::VPC_APPLY_TIME);

    vpcs.insert(vpcs_by_round[round_ind]);

    printf(
        "\n## Round %d/%d: Completed round of VPC generation (exit reason: %s). # cuts generated = %d.\n",
        round_ind + 1, params.get(ROUNDS),
        CglVPC::ExitReasonName[static_cast<int>(exitReason)].c_str(),
        vpcs_by_round[round_ind].sizeCuts());
    fflush(stdout);
    printf("Initial obj value: %1.6f. New obj value: %s. Disj lb: %s. ##\n",
        boundInfo.lp_obj, stringValue(solver->getObjValue(), "%1.6f").c_str(),
        stringValue(boundInfo.best_disj_obj, "%1.6f").c_str());

    // Exit early from rounds of cuts if no cuts generated or solver is not optimal
    if (gen.num_cuts == 0 || !solver->isProvenOptimal()
        || isInfinity(std::abs(solver->getObjValue()))
        || exitReason == CglVPC::ExitReason::OPTIMAL_SOLUTION_FOUND_EXIT)
      break;
  } // loop over rounds of cuts
  if (round_ind < num_rounds)
    num_rounds = round_ind+1;

  printf(
      "\n## Finished VPC generation with %d cuts. Initial obj value: %s. Final obj value: %s. Disj lb: %s. ##\n",
      boundInfo.num_vpc,
      stringValue(boundInfo.lp_obj, "%1.6f").c_str(),
      stringValue(boundInfo.vpc_obj, "%1.6f").c_str(),
      stringValue(boundInfo.best_disj_obj, "%1.6f").c_str());

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
  if (boundInfo.num_vpc > 0) {
    std::string fileWithCuts = filename_stub + "_cuts";
    solver->writeMps(fileWithCuts.c_str());
  }
#endif

  // Do analyses in preparation for printing
  setCutInfo(cutInfo, num_rounds, cutInfoVec.data());
  analyzeStrength(params, 
      params.get(GOMORY) == 0 ? NULL : GMICSolver,
      params.get(GOMORY) <= 0 ? solver : VPCSolver, 
      params.get(GOMORY) >= 0 ? solver : VPCSolver, 
      cutInfoGMICs, cutInfo, 
      params.get(GOMORY) == 0 ? NULL : &gmics, 
      &vpcs,
      boundInfo, cut_output);
  analyzeBB(params, info_nocuts, info_mycuts, info_allcuts, bb_output);
  return wrapUp(0, argc, argv);
} /* main */

/**
 * Call this early to print welcome message, etc.
 */
int startUp(int argc, char** argv) {
  int status = 0;

  // Input handling
  printf("## V-Polyhedral Disjunctive Cuts ##\n");
  printf("# Aleksandr M. Kazachkov\n");
  printf("# Based on joint work with Egon Balas\n");
  for (int i = 0; i < argc; i++) {
    std::cout << argv[i] << " ";
  }
  std::cout << std::endl;
  time(&start_time_t);
  struct tm* start_timeinfo = localtime(&start_time_t);
  snprintf(start_time_string, sizeof(start_time_string) / sizeof(char), "%s", asctime(start_timeinfo));
  printf("Start time: %s\n", start_time_string);

  /*
  printf("\n## Version Information ##\n");
#ifdef VPC_VERSION
  printf("# VPC Version %s\n", VPC_VERSION_STRING.substr(0,8).c_str());
#endif
#ifdef VPC_CBC_VERSION
  printf("# Cbc Version %s\n", CBC_VERSION_STRING.substr(0,8).c_str());
#endif
#ifdef VPC_CLP_VERSION
  printf("# Clp Version %s\n", CLP_VERSION_STRING.substr(0,8).c_str());
#endif
#ifdef USE_GUROBI
  printf("# Gurobi Version %s\n", GUROBI_VERSION_STRING.c_str());
#endif
#ifdef USE_CPLEX
  printf("# CPLEX Version %s\n", CPLEX_VERSION_STRING.c_str());
#endif
  printf("\n");
  */

  status = processArgs(argc, argv);
  if (status) { return status; }

  // Get instance file
  if (params.get(stringParam::FILENAME).empty()) {
    error_msg(errorstring,
        "No instance file provided. Use -f /path/to/instance.mps to specify the instance.\n");
    exit(1);
  }
  printf("Instance file: %s\n", params.get(stringParam::FILENAME).c_str());
  
  if (parseFilename(dir, instname, in_file_ext, params.get(stringParam::FILENAME), params.logfile) != 0) {
    error_msg(errorstring,
        "Unable to parse filename: %s. Found: dir=\"%s\", instname=\"%s\",ext=\"%s\".\n",
        params.get(stringParam::FILENAME).c_str(), dir.c_str(),
        instname.c_str(), in_file_ext.c_str());
    exit(1);
  }
  filename_stub = dir + "/" + instname;

  // Prepare logfile
  const std::string logname = params.get(stringParam::LOGFILE);
  bool logexists = false; 
  if (!logname.empty()) {
    logexists = fexists(logname.c_str());
    params.logfile = fopen(logname.c_str(), "a");
  }

  // Read opt value (if not yet inputted)
  if (!isInfinity(params.get(doubleParam::IP_OBJ))) {
    boundInfo.ip_obj = params.get(doubleParam::IP_OBJ);
  }
  if (isInfinity(boundInfo.ip_obj) && !params.get(stringParam::OPTFILE).empty()) {
    std::string optfile = params.get(stringParam::OPTFILE);
    std::string csvext = ".csv";
    if (optfile.size() > csvext.size() && optfile.compare(optfile.size()-csvext.size(), csvext.size(), csvext) == 0) {
  #ifdef TRACE
      fprintf(stdout, "Reading objective information from \"%s\".\n", optfile.c_str());
      //std::cout << "Reading objective information from \"" + params.get(stringParam::OPTFILE) + "\"" << std::endl;
  #endif
      boundInfo.ip_obj = getObjValueFromFile(optfile, params.get(stringParam::FILENAME), params.logfile);
      params.set(doubleParam::IP_OBJ, boundInfo.ip_obj);
      fprintf(stdout, "Best known objective value is %s.\n", stringValue(boundInfo.ip_obj, "%f").c_str());
      //std::cout << "Best known objective value is " << boundInfo.ip_obj << std::endl;
      if (isInfinity(boundInfo.ip_obj)) {
        warning_msg(warnstring, "Did not find objective value.\n");
      }
    }
  }
  if (params.logfile != NULL) {
    if (!logexists) {
      if (params.get(PREPROCESS) == 0) {
        printHeader(params, OverallTimeStatsName);
      } else {
        printPreprocessingHeader(params);
      }
    }
    fprintf(params.logfile, "%s,", instname.c_str());
    if (params.get(PREPROCESS) == 0) {
      printParams(params, params.logfile, 2); // only values
    }
    fflush(params.logfile);
  }

  return status;
} /* startUp */

/**
 * Close the logfile and print to it
 */
int wrapUp(int retCode, int argc, char** argv) {
  const int exitReasonInt = static_cast<int>(exitReason);

  time(&end_time_t);
  struct tm* end_timeinfo = localtime(&end_time_t);
  char end_time_string[25];
  snprintf(end_time_string, sizeof(end_time_string) / sizeof(char), "%s", asctime(end_timeinfo));

  FILE* logfile = params.logfile;
  if (logfile != NULL) {
    if (params.get(PREPROCESS) == 0) {
      // Bound and gap info
      printBoundAndGapInfo(boundInfo, params.logfile);
      // B&B info
      printSummaryBBInfo({info_nocuts, info_mycuts}, params.logfile, params.get(BB_RUNS) == 0);
      // Orig prob
      printOrigProbInfo(origSolver, params.logfile);
      // Post-cut prob
      printPostCutProbInfo(solver, cutInfoGMICs, cutInfo, params.logfile);
      // Disj info
      printDisjInfo(disjInfo, params.logfile);
      // Cut, obj, fail info
      printCutInfo(cutInfoGMICs, cutInfo, params.logfile);
      // Full B&B info
      printFullBBInfo({info_nocuts, info_mycuts}, params.logfile, params.get(BB_RUNS) == 0);
      // Print time info
      timer.print(params.logfile, 2); // only values
    } else {
    }

#ifdef VPC_VERSION
    fprintf(logfile, "%s,", VPC_VERSION_STRING.substr(0,8).c_str());
#else
    fprintf(logfile, ",");
#endif
#ifdef VPC_CBC_VERSION
    fprintf(logfile, "%s,", CBC_VERSION_STRING.substr(0,8).c_str());
#else
    fprintf(logfile, ",");
#endif
#ifdef VPC_CLP_VERSION
    fprintf(logfile, "%s,", CLP_VERSION_STRING.substr(0,8).c_str());
#else
    fprintf(logfile, ",");
#endif
#ifdef USE_GUROBI
    fprintf(logfile, "%s,", GUROBI_VERSION_STRING.c_str());
#else
    fprintf(logfile, ",");
#endif
#ifdef USE_CPLEX
    fprintf(logfile, "%s,", CPLEX_VERSION_STRING.c_str());
#else
    fprintf(logfile, ",");
#endif

    fprintf(logfile, "%s,", CglVPC::ExitReasonName[exitReasonInt].c_str());
    fprintf(logfile, "%s,", end_time_string);
    fprintf(logfile, "%.f,", difftime(end_time_t, start_time_t));
    fprintf(logfile, "%s,", instname.c_str());
    fprintf(logfile, "DONE\n");
    fclose(logfile); // closes params.logfile
  }

#ifdef TRACE
  // Print parameters
  printf("\n## Parameter values ##\n");
  printParams(params, stdout, 1);
  printf("\n");
  printParams(params, stdout, 2);
  printf("\n");

  int NAME_WIDTH = 25;
  // Print time info
  printf("\n## Time information ##\n");
  for (int i = 0; i < NUM_TIME_STATS; i++) {
    std::string name = OverallTimeStatsName[i];
    if (params.get(intParam::DISJ_TERMS) == 0 && name.compare(0, 3, "VPC") == 0)
      continue;
    if (params.get(intParam::BB_RUNS) == 0 && name.compare(0, 2, "BB") == 0)
      continue;
    printf("%-*.*s%s\n", NAME_WIDTH, NAME_WIDTH, (name).c_str(),
        stringValue(timer.get_time(name), "%.3f").c_str());
  }
#endif

  // Print results from adding cuts
  printf("%s", cut_output.c_str());

  // Print branch-and-bound results
  printf("%s", bb_output.c_str());

  printf("\n## Exiting VPC generation with reason %s. ##\n", CglVPC::ExitReasonName[exitReasonInt].c_str());
#ifdef VPC_VERSION
  printf("VPC Version: %s\n", VPC_VERSION_STRING.substr(0,8).c_str());
#endif
#ifdef VPC_CBC_VERSION
  printf("Cbc Version: %s\n", CBC_VERSION_STRING.substr(0,8).c_str());
#endif
#ifdef VPC_CLP_VERSION
  printf("Clp Version: %s\n", CLP_VERSION_STRING.substr(0,8).c_str());
#endif
#ifdef USE_GUROBI
  printf("Gurobi Version: %s\n", GUROBI_VERSION_STRING.c_str());
#endif
#ifdef USE_CPLEX
  printf("CPLEX Version: %s\n", CPLEX_VERSION_STRING.c_str());
#endif
  printf("Instance: %s\n", instname.c_str());
  if (!params.get(stringParam::LOGFILE).empty()) {
    printf("Log: %s\n", params.get(stringParam::LOGFILE).c_str());
  } else {
    printf("Log: stdout\n");
  }
  printf("Start time: %s\n", start_time_string);
  printf("End time: %s\n", end_time_string);
  printf("Elapsed time: %.f seconds\n", difftime(end_time_t, start_time_t));
  { // Print command used to repeat this run
    if (argc > 0) printf("Command: ");
    for (int i = 0; i < argc; i++) {
      printf("%s ", argv[i]);
    }
    if (argc > 0) printf("\n");
  }

  if (solver) {
    delete solver;
  }
  if (origSolver) {
    delete origSolver;
  }
  if (GMICSolver) {
    delete GMICSolver;
  }
  if (VPCSolver) {
    delete VPCSolver;
  }

  return retCode;
} /* wrapUp */

void doCustomRoundOfCuts(int round_ind, OsiCuts& vpcs, CglVPC& gen, int& num_disj) {
  std::vector<Disjunction*> disjVec;
  printf("\n## Setting up disjunction(s) ##\n");
  gen.timer.start_timer(CglVPC::VPCTimeStatsName[static_cast<int>(CglVPC::VPCTimeStats::DISJ_GEN_TIME)]);
  CglVPC::ExitReason setDisjExitReason = setDisjunctions(disjVec, solver, params, CglVPC::VPCMode::SPLITS);
  gen.timer.end_timer(CglVPC::VPCTimeStatsName[static_cast<int>(CglVPC::VPCTimeStats::DISJ_GEN_TIME)]);
  const int numDisj = disjVec.size();

  // If integer-optimal solution was found, all disjunctions but one will have been deleted
  if (setDisjExitReason == CglVPC::ExitReason::OPTIMAL_SOLUTION_FOUND_EXIT) {
    warning_msg(warnstr,
        "An integer (optimal) solution was found prior while getting disjunction. "
        "We will generate between n and 2n cuts, restricting the value of each variable.\n");
    const double* solution = disjVec[0]->integer_sol.data();
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
              CglVPC::ObjectiveType::ONE_SIDED);
        }
      }
    } // iterate over columns and add optimality cut if needed
    exitReason = CglVPC::ExitReason::OPTIMAL_SOLUTION_FOUND_EXIT;
    boundInfo.num_vpc += gen.num_cuts;
    updateCutInfo(cutInfoVec[round_ind], gen);
  } // check if integer-optimal solution
  else if (setDisjExitReason == CglVPC::ExitReason::SUCCESS_EXIT && numDisj > 0) {
    const int cutLimit = std::ceil(
        gen.getCutLimit(params.get(CUTLIMIT),
            solver->getFractionalIndices().size()) / numDisj); // distribute cut limit over the disjunctions
    gen.setupRepeatedUse(true);
    for (Disjunction* disj : disjVec) {
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
      updateCutInfo(cutInfoVec[round_ind], gen);
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

  // Delete disjunctions
  for (Disjunction* disj : disjVec) {
    if (disj) {
      delete disj;
      disj = NULL;
    }
  }
} /* doCustomRoundOfCuts */

/**
 * See params.hpp for descriptions of the parameters
 */
int processArgs(int argc, char** argv) {
  // Handle inputs
  // struct option declared in getopt.h
  // name: name of the long option
  // has_arg: 0,1,2 for none, required, or optional
  // *flag: how results are returned; if NULL, getopt_long() returns val (e.g., can be the equivalent short option character), and o/w getopt_long() returns 0, and flag points to a var which is set to val if the option is found, but left unchanged if the option is not found
  // val: value to return, or to load into the variable pointed to by flag
  const char* const short_opts = "b:B:c:d:D:f:g:hi:l:m:o:r:R:s:S:t:T:v:";
  const struct option long_opts[] =
  {
      {"bb_runs",               required_argument, 0, 'b'},
      {"bb_mode",               required_argument, 0, 'b'*'1'},
      {"bb_strategy",           required_argument, 0, 'B'},
      {"bb_timelimit",          required_argument, 0, 'B'*'0'},
      {"cbc",                   required_argument, 0, 'B'*'1'*'1'},
      {"cplex",                 required_argument, 0, 'B'*'1'*'2'},
      {"gurobi",                required_argument, 0, 'B'*'1'*'3'},
      {"user_cuts",             required_argument, 0, 'B'*'2'},
      {"use_best_bound",        required_argument, 0, 'B'*'3'},
      {"bb_all_cuts",           required_argument, 0, 'B'*'4'},
      {"bb_gmics",              required_argument, 0, 'B'*'5'},
      {"bb_heuristics",         required_argument, 0, 'B'*'6'},
      {"bb_presolve",           required_argument, 0, 'B'*'7'},
      {"bb_strong_branching",   required_argument, 0, 'B'*'8'},
      {"cutlimit",              required_argument, 0, 'c'},
      {"disj_terms",            required_argument, 0, 'd'},
      {"disj_options",          required_argument, 0, 'D'},
      {"file",                  required_argument, 0, 'f'},
      {"gomory",                required_argument, 0, 'g'},
      {"help",                  no_argument,       0, 'h'},
      {"ip_obj",                required_argument, 0, 'i'},
      {"logfile",               required_argument, 0, 'l'},
      {"mode",                  required_argument, 0, 'm'},
      {"optfile",               required_argument, 0, 'o'},
      {"partial_bb_strategy",   required_argument, 0, 's'},
      {"partial_bb_num_strong", required_argument, 0, 'S'},
      {"partial_bb_timelimit",  required_argument, 0, 'T'},
      {"preprocess",            required_argument, 0, 'p'*'1'},
      {"rounds",                required_argument, 0, 'r'},
      {"prlp_timelimit",        required_argument, 0, 'R'},
      {"solfile",               required_argument, 0, 's'*'1'},
      {"temp",                  required_argument, 0, 't'*'1'},
      {"tikz",                  required_argument, 0, 't'*'2'},
      {"timelimit",             required_argument, 0, 't'},
      {"use_all_ones",          required_argument, 0, 'u'*'1'},
      {"use_disj_lb",           required_argument, 0, 'u'*'2'},
      {"use_iter_bilinear",     required_argument, 0, 'u'*'3'},
      {"use_tight_points",      required_argument, 0, 'u'*'4'},
      {"use_tight_rays",        required_argument, 0, 'u'*'5'},
      {"use_unit_vectors",      required_argument, 0, 'u'*'6'},
      {"verbosity",             required_argument, 0, 'v'},
      {nullptr, no_argument, nullptr, 0}
  };

  int inp;
  int status = 0;
  while ((inp = getopt_long(argc, argv, short_opts, long_opts, nullptr)) != -1) {
    switch (inp) {
      case 'b': {
                  int val;
                  intParam param = intParam::BB_RUNS;
                  if (!parseInt(optarg, val)) {
                    error_msg(errorstring, "Error reading %s. Given value: %s.\n", params.name(param).c_str(), optarg);
                    exit(1);
                  }
                  params.set(param, val);
                  break;
                }
      case 'b'*'1': {
                  int val;
                  intParam param = intParam::BB_MODE;
                  if (!parseInt(optarg, val)) {
                    error_msg(errorstring, "Error reading %s. Given value: %s.\n", params.name(param).c_str(), optarg);
                    exit(1);
                  }
                  params.set(param, val);
                  break;
                }
      case 'B': {
                  int val;
                  intParam param = intParam::BB_STRATEGY;
                  if (!parseInt(optarg, val)) {
                    error_msg(errorstring, "Error reading %s. Given value: %s.\n", params.name(param).c_str(), optarg);
                    exit(1);
                  }
                  params.set(param, val);
                  if (use_bb_option(val, BB_Strategy_Options::gurobi)) {
#ifndef USE_GUROBI
                    error_msg(errorstring, "Requesting to use Gurobi, but macro USE_GUROBI not set.\n");
                    exit(1);
#endif
                  }
                  if (use_bb_option(val, BB_Strategy_Options::cplex)) {
#ifndef USE_CPLEX
                    error_msg(errorstring, "Requesting to use CPLEX, but macro USE_CPLEX not set.\n");
                    exit(1);
#endif
                  }
                  break;
                }
      case 'B'*'0': {
                  double val;
                  doubleParam param = doubleParam::BB_TIMELIMIT;
                  if (!parseDouble(optarg, val)) {
                    error_msg(errorstring, "Error reading %s. Given value: %s.\n", params.name(param).c_str(), optarg);
                    exit(1);
                  }
                  params.set(param, val);
                  break;
                }
      case 'B'*'1'*'1': {
                  intParam param = intParam::BB_STRATEGY;
                  int val;
                  if (!parseInt(optarg, val) || (val != 0 && val != 1)) {
                    error_msg(errorstring, "Error reading binary parameter cbc. Given value: %s.\n", optarg);
                    exit(1);
                  }
                  switch (val) {
                    case 0:
                      val = disable_bb_option(params.get(param), BB_Strategy_Options::cbc);
                      break;
                    case 1:
                      val = enable_bb_option(params.get(param), BB_Strategy_Options::cbc);
                      val = disable_bb_option(val, BB_Strategy_Options::cplex);
                      val = disable_bb_option(val, BB_Strategy_Options::gurobi);
                      break;
                  }
                  params.set(param, val);
                  break;
                }
      case 'B'*'1'*'2': {
                  intParam param = intParam::BB_STRATEGY;
                  int val;
                  if (!parseInt(optarg, val) || (val != 0 && val != 1)) {
                    error_msg(errorstring, "Error reading binary parameter cplex. Given value: %s.\n", optarg);
                    exit(1);
                  }
                  switch (val) {
                    case 0:
                      val = disable_bb_option(params.get(param), BB_Strategy_Options::cplex);
                      break;
                    case 1:
                      val = enable_bb_option(params.get(param), BB_Strategy_Options::cplex);
                      val = disable_bb_option(val, BB_Strategy_Options::cbc);
                      val = disable_bb_option(val, BB_Strategy_Options::gurobi);
                      break;
                  }
                  params.set(param, val);
                  break;
                }
      case 'B'*'1'*'3': {
                  intParam param = intParam::BB_STRATEGY;
                  int val;
                  if (!parseInt(optarg, val) || (val != 0 && val != 1)) {
                    error_msg(errorstring, "Error reading binary parameter gurobi. Given value: %s.\n", optarg);
                    exit(1);
                  }
                  switch (val) {
                    case 0:
                      val = disable_bb_option(params.get(param), BB_Strategy_Options::gurobi);
                      break;
                    case 1:
                      val = enable_bb_option(params.get(param), BB_Strategy_Options::gurobi);
                      val = disable_bb_option(val, BB_Strategy_Options::cbc);
                      val = disable_bb_option(val, BB_Strategy_Options::cplex);
                      break;
                  }
                  params.set(param, val);
                  break;
                }
      case 'B'*'2': {
                  intParam param = intParam::BB_STRATEGY;
                  int val;
                  if (!parseInt(optarg, val) || (val != 0 && val != 1)) {
                    error_msg(errorstring, "Error reading binary parameter user_cuts. Given value: %s.\n", optarg);
                    exit(1);
                  }
                  switch (val) {
                    case 0:
                      val = disable_bb_option(params.get(param), BB_Strategy_Options::user_cuts);
                      break;
                    case 1:
                      val = enable_bb_option(params.get(param), BB_Strategy_Options::user_cuts);
                      break;
                  }
                  params.set(param, val);
                  break;
                }
      case 'B'*'3': {
                  intParam param = intParam::BB_STRATEGY;
                  int val;
                  if (!parseInt(optarg, val) || (val != 0 && val != 1)) {
                    error_msg(errorstring, "Error reading binary parameter use_best_bound. Given value: %s.\n", optarg);
                    exit(1);
                  }
                  switch (val) {
                    case 0:
                      val = disable_bb_option(params.get(param), BB_Strategy_Options::use_best_bound);
                      break;
                    case 1:
                      val = enable_bb_option(params.get(param), BB_Strategy_Options::use_best_bound);
                      break;
                  }
                  val = enable_bb_option(val, BB_Strategy_Options::use_best_bound);
                  params.set(param, val);
                  break;
                }
      case 'B'*'4': {
                  int val;
                  intParam param = intParam::BB_STRATEGY;
                  if (!parseInt(optarg, val) || (val != 0 && val != 1)) {
                    error_msg(errorstring, "Error reading binary parameter bb_all_cuts. Given value: %s.\n", optarg);
                    exit(1);
                  }
                  switch (val) {
                    case 0:
                      val = enable_bb_option(params.get(param), BB_Strategy_Options::all_cuts_off);
                      val = disable_bb_option(val, BB_Strategy_Options::all_cuts_on);
                      break;
                    case 1:
                      val = enable_bb_option(params.get(param), BB_Strategy_Options::all_cuts_on);
                      val = disable_bb_option(val, BB_Strategy_Options::all_cuts_off);
                      break;
                  }
                  params.set(param, val);
                  break;
                }
      case 'B'*'5': {
                  int val;
                  intParam param = intParam::BB_STRATEGY;
                  if (!parseInt(optarg, val) || (val != 0 && val != 1)) {
                    error_msg(errorstring, "Error reading binary parameter bb_gmics. Given value: %s.\n", optarg);
                    exit(1);
                  }
                  switch (val) {
                    case 0:
                      val = enable_bb_option(params.get(param), BB_Strategy_Options::gmics_off);
                      val = disable_bb_option(val, BB_Strategy_Options::gmics_on);
                      break;
                    case 1:
                      val = enable_bb_option(params.get(param), BB_Strategy_Options::gmics_on);
                      val = disable_bb_option(val, BB_Strategy_Options::gmics_off);
                      break;
                  }
                  params.set(param, val);
                  break;
                }
      case 'B'*'6': {
                  int val;
                  intParam param = intParam::BB_STRATEGY;
                  if (!parseInt(optarg, val) || (val != 0 && val != 1)) {
                    error_msg(errorstring, "Error reading binary parameter bb_heuristics. Given value: %s.\n", optarg);
                    exit(1);
                  }
                  switch (val) {
                    case 0:
                      val = enable_bb_option(params.get(param), BB_Strategy_Options::presolve_off);
                      val = disable_bb_option(val, BB_Strategy_Options::presolve_on);
                      break;
                    case 1:
                      val = enable_bb_option(params.get(param), BB_Strategy_Options::presolve_on);
                      val = disable_bb_option(val, BB_Strategy_Options::presolve_off);
                      break;
                  }
                  params.set(param, val);
                  break;
                }
      case 'B'*'7': {
                  int val;
                  intParam param = intParam::BB_STRATEGY;
                  if (!parseInt(optarg, val) || (val != 0 && val != 1)) {
                    error_msg(errorstring, "Error reading binary parameter bb_presolve. Given value: %s.\n", optarg);
                    exit(1);
                  }
                  switch (val) {
                    case 0:
                      val = enable_bb_option(params.get(param), BB_Strategy_Options::heuristics_off);
                      val = disable_bb_option(val, BB_Strategy_Options::heuristics_on);
                      break;
                    case 1:
                      val = enable_bb_option(params.get(param), BB_Strategy_Options::heuristics_on);
                      val = disable_bb_option(val, BB_Strategy_Options::heuristics_off);
                      break;
                  }
                  params.set(param, val);
                  break;
                }
      case 'B'*'8': {
                  int val;
                  intParam param = intParam::BB_STRATEGY;
                  if (!parseInt(optarg, val) || (val != 0 && val != 1)) {
                    error_msg(errorstring, "Error reading binary parameter bb_strong_branching. Given value: %s.\n", optarg);
                    exit(1);
                  }
                  switch (val) {
                    case 0:
                      val = disable_bb_option(params.get(param), BB_Strategy_Options::strong_branching_on);
                      break;
                    case 1:
                      val = enable_bb_option(params.get(param), BB_Strategy_Options::strong_branching_on);
                      break;
                  }
                  params.set(param, val);
                  break;
                }
      case 'c': {
                  int val;
                  intParam param = intParam::CUTLIMIT;
                  if (!parseInt(optarg, val)) {
                    error_msg(errorstring, "Error reading %s. Given value: %s.\n", params.name(param).c_str(), optarg);
                    exit(1);
                  }
                  params.set(param, val);
                  break;
                }
      case 'd': {
                  int val;
                  intParam param = intParam::DISJ_TERMS;
                  if (!parseInt(optarg, val)) {
                    error_msg(errorstring, "Error reading %s. Given value: %s.\n", params.name(param).c_str(), optarg);
                    exit(1);
                  }
                  params.set(param, val);
                  break;
                }
      case 'D': {
                  params.set(stringParam::DISJ_OPTIONS, optarg);
                  break;
                }
      case 'f': {
                  params.set(stringParam::FILENAME, optarg);
                  break;
                }
      case 'g': {
                  int val;
                  intParam param = intParam::GOMORY;
                  if (!parseInt(optarg, val)) {
                    error_msg(errorstring, "Error reading %s. Given value: %s.\n", params.name(param).c_str(), optarg);
                    exit(1);
                  }
                  params.set(param, val);
                  break;
                }
      case 'i': {
                 double val;
                 doubleParam param = doubleParam::IP_OBJ;
                 if (!parseDouble(optarg, val)) {
                   error_msg(errorstring, "Error reading %s. Given value: %s.\n", params.name(param).c_str(), optarg);
                   exit(1);
                 }
                 params.set(param, val);
                 break;
               }
      case 'l': {
                  params.set(stringParam::LOGFILE, optarg);
                  break;
                }
      case 'm': {
                 int val;
                 intParam param = intParam::MODE;
                 if (!parseInt(optarg, val)) {
                   error_msg(errorstring, "Error reading %s. Given value: %s.\n", params.name(param).c_str(), optarg);
                   exit(1);
                 }
                 params.set(param, val);
                 break;
               }
      case 'o': {
                  params.set(stringParam::OPTFILE, optarg);
                  break;
                }
      case 's'*'1': {
                      params.set(stringParam::SOLFILE, optarg);
                      break;
                    }
      case 's': {
                  int val;
                  intParam param = intParam::PARTIAL_BB_STRATEGY;
                  if (!parseInt(optarg, val)) {
                    error_msg(errorstring, "Error reading %s. Given value: %s.\n", params.name(param).c_str(), optarg);
                    exit(1);
                  }
                  params.set(param, val);
                  break;
                }
      case 'S': {
                  int val;
                  intParam param = intParam::PARTIAL_BB_NUM_STRONG;
                  if (!parseInt(optarg, val)) {
                    error_msg(errorstring, "Error reading %s. Given value: %s.\n", params.name(param).c_str(), optarg);
                    exit(1);
                  }
                  params.set(param, val);
                  break;
                }
      case 'T': {
                  double val;
                  doubleParam param = doubleParam::PARTIAL_BB_TIMELIMIT;
                  if (!parseDouble(optarg, val)) {
                    error_msg(errorstring, "Error reading %s. Given value: %s.\n", params.name(param).c_str(), optarg);
                    exit(1);
                  }
                  params.set(param, val);
                  break;
                }
      case 'p'*'1': {
                  int val;
                  intParam param = intParam::PREPROCESS;
                  if (!parseInt(optarg, val)) {
                    error_msg(errorstring, "Error reading %s. Given value: %s.\n", params.name(param).c_str(), optarg);
                    exit(1);
                  }
                  params.set(intParam::PREPROCESS, val);
                  params.set(intParam::BB_RUNS, 1);
                  params.set(intParam::BB_MODE, 001);
                  const int base_bb_strategy =
                    get_bb_option_value({
                        BB_Strategy_Options::user_cuts,
                        BB_Strategy_Options::presolve_off
                    });
#ifdef USE_GUROBI
                  params.set(intParam::BB_STRATEGY,
                      base_bb_strategy + get_bb_option_value({BB_Strategy_Options::gurobi})
                  ); // previously 10776
#elif USE_CPLEX
                  params.set(intParam::BB_STRATEGY,
                      base_bb_strategy + get_bb_option_value({BB_Strategy_Options::cplex})
                  ); // previously 10772
#elif USE_CBC
                  params.set(intParam::BB_STRATEGY,
                      base_bb_strategy + get_bb_option_value({BB_Strategy_Options::cbc})
                  ); // previously 10770
#else
                  params.set(intParam::BB_STRATEGY,
                      base_bb_strategy + get_bb_option_value({BB_Strategy_Options::off})
                  ); // previously 10768
#endif
                  params.set(doubleParam::TIMELIMIT, 7200);
                  break;
                }
      case 'R': {
                  double val;
                  doubleParam param = doubleParam::PRLP_TIMELIMIT;
                  if (!parseDouble(optarg, val)) {
                    error_msg(errorstring, "Error reading %s. Given value: %s.\n", params.name(param).c_str(), optarg);
                    exit(1);
                  }
                  params.set(param, val);
                  break;
                }
      case 'r': {
                   int val;
                   intParam param = intParam::ROUNDS;
                   if (!parseInt(optarg, val)) {
                     error_msg(errorstring, "Error reading %s. Given value: %s.\n", params.name(param).c_str(), optarg);
                     exit(1);
                   }
                   params.set(param, val);
                   break;
                 }
      case 't'*'1': {
                      int val;
                      intParam param = intParam::TEMP;
                      if (!parseInt(optarg, val)) {
                        error_msg(errorstring, "Error reading %s. Given value: %s.\n", params.name(param).c_str(), optarg);
                        exit(1);
                      }

                      // Validate temp parameter when tikz_string is being used
                      if (use_temp_option(std::abs(val), TempOptions::GEN_TIKZ_STRING_WITH_VPCS)
                          || use_temp_option(std::abs(val), TempOptions::GEN_TIKZ_STRING_WITH_GMICS)
                          || use_temp_option(std::abs(val), TempOptions::GEN_TIKZ_STRING_AND_RETURN)
                          || use_temp_option(std::abs(val), TempOptions::GEN_TIKZ_STRING_AND_EXIT)) {
                        val = enable_temp_option(val, TempOptions::GEN_TIKZ_STRING); // TODO ensure this works with negative values
                      }

                      params.set(param, val);
                      break;
                    }
      case 't'*'2': {
                      intParam param = intParam::TEMP;
                      int val = params.get(param); // old param value, so we do not overwrite it
                      if (strcmp(optarg, "GEN_TIKZ_STRING") == 0) {
                        val = enable_temp_option(val, TempOptions::GEN_TIKZ_STRING);
                      }
                      else if (strcmp(optarg, "GEN_TIKZ_STRING_WITH_VPCS") == 0) {
                        val = enable_temp_option(val, TempOptions::GEN_TIKZ_STRING_WITH_VPCS);
                        val = enable_temp_option(val, TempOptions::GEN_TIKZ_STRING);
                      }
                      else if (strcmp(optarg, "GEN_TIKZ_STRING_WITH_GMICS") == 0) {
                        val = enable_temp_option(val, TempOptions::GEN_TIKZ_STRING_WITH_GMICS);
                        val = enable_temp_option(val, TempOptions::GEN_TIKZ_STRING);
                      }
                      else if (strcmp(optarg, "GEN_TIKZ_STRING_AND_RETURN") == 0) {
                        val = enable_temp_option(val, TempOptions::GEN_TIKZ_STRING_AND_RETURN);
                        val = enable_temp_option(val, TempOptions::GEN_TIKZ_STRING);
                      }
                      else if (strcmp(optarg, "GEN_TIKZ_STRING_AND_EXIT") == 0) {
                        val = enable_temp_option(val, TempOptions::GEN_TIKZ_STRING_AND_EXIT);
                        val = enable_temp_option(val, TempOptions::GEN_TIKZ_STRING);
                      }
                      else {
                        error_msg(errorstring,
                            "Error reading %s. Given value: %s. Possible values can be found in VPCParameters.hpp under the enum class TempOptions.\n",
                            "tikz parameter", optarg);
                        exit(1);
                      }

                      params.set(param, val);
                      break;
                    }
      case 't': {
                  double val;
                  doubleParam param = doubleParam::TIMELIMIT;
                  if (!parseDouble(optarg, val)) {
                    error_msg(errorstring, "Error reading %s. Given value: %s.\n", params.name(param).c_str(), optarg);
                    exit(1);
                  }
                  params.set(param, val);
                  break;
                }
      case 'u'*'1': {
                      int val;
                      intParam param = intParam::USE_ALL_ONES;
                      if (!parseInt(optarg, val)) {
                        error_msg(errorstring, "Error reading %s. Given value: %s.\n", params.name(param).c_str(), optarg);
                        exit(1);
                      }
                      params.set(param, val);
                      break;
                    }
      case 'u'*'2': {
                      int val;
                      intParam param = intParam::USE_DISJ_LB;
                      if (!parseInt(optarg, val)) {
                        error_msg(errorstring, "Error reading %s. Given value: %s.\n", params.name(param).c_str(), optarg);
                        exit(1);
                      }
                      params.set(param, val);
                      break;
                    }
      case 'u'*'3': {
                      int val;
                      intParam param = intParam::USE_ITER_BILINEAR;
                      if (!parseInt(optarg, val)) {
                        error_msg(errorstring, "Error reading %s. Given value: %s.\n", params.name(param).c_str(), optarg);
                        exit(1);
                      }
                      params.set(param, val);
                      break;
                    }
      case 'u'*'4': {
                      int val;
                      intParam param = intParam::USE_TIGHT_POINTS;
                      if (!parseInt(optarg, val)) {
                        error_msg(errorstring, "Error reading %s. Given value: %s.\n", params.name(param).c_str(), optarg);
                        exit(1);
                      }
                      params.set(param, val);
                      break;
                    }
      case 'u'*'5': {
                      int val;
                      intParam param = intParam::USE_TIGHT_RAYS;
                      if (!parseInt(optarg, val)) {
                        error_msg(errorstring, "Error reading %s. Given value: %s.\n", params.name(param).c_str(), optarg);
                        exit(1);
                      }
                      params.set(param, val);
                      break;
                    }
      case 'u'*'6': {
                      int val;
                      intParam param = intParam::USE_UNIT_VECTORS;
                      if (!parseInt(optarg, val)) {
                        error_msg(errorstring, "Error reading %s. Given value: %s.\n", params.name(param).c_str(), optarg);
                        exit(1);
                      }
                      params.set(param, val);
                      break;
                    }
      case 'v': {
                   int val;
                   intParam param = intParam::VERBOSITY;
                   if (!parseInt(optarg, val)) {
                     error_msg(errorstring, "Error reading %s. Given value: %s.\n", params.name(param).c_str(), optarg);
                     exit(1);
                   }
                   params.set(param, val);
                   break;
                 }
      case 'h':
      case '?':
      default: {
                 // print help
                std::string helpstring;
                helpstring += "\n## DESCRIPTION ##\n";
                helpstring += "Code for generating V-polyhedral disjunctive cuts.\n";
                helpstring += "\n## OPTIONS ##\n";
                helpstring += "-h, --help\n\tPrint this help message.\n";
                helpstring += "--preprocess\n\tSet default preprocessing options (temp=1, bb_mode=001, bb_runs=7, timelimit=7200).\n";
                helpstring += "--temp\n\tSet temporary options (e.g., value of 1 = do preprocessing on instance).\n";
                helpstring += "\n# Input/output #\n";
                helpstring += "-f file, --file=file\n\tFilename.\n";
                helpstring += "-i val, --ip_obj=val\n\tValue of integer optimum for this instance (takes precedence over -o/--optfile).\n";
                helpstring += "-l logfile, --logfile=logfile\n\tWhere to print log messages.\n";
                helpstring += "-o optfile, --optfile=optfile\n\tWhere to find integer optimum value information (a csv file formatted as \"instance_name,value\" on each row, or a .sol/.mst file to be read by CPLEX/Gurobi).\n";
                helpstring += "--solfile=solfile\n\tWhere to find integer optimum solution information (e.g., a mst/sol file produced by Gurobi/CPLEX/etc).\n";
                helpstring += "-v level, --verbosity=level\n\tVerbosity level (0: print little, 1: let solver output be visible).\n";
                helpstring += "\n# General VPC options #\n";
                helpstring += "-c num cuts, --cutlimit=num cuts\n\tMaximum number of cuts to generate (0+ = as given, -k = k * # fractional variables at root).\n";
                helpstring += "-d num terms, --disj_terms=num terms\n\tMaximum number of disjunctive terms or disjunctions to generate (depending on mode).\n";
                helpstring += "-g -1/0/1, --gomory=-1/0/1\n\t0: do not use Gomory cuts before generating VPCs, +/-1: generate Gomory cuts before generating VPCs (-1: only gen, +1: also apply to LP).\n";
                helpstring += "-m mode, --mode=mode\n\tMode for generating disjunction(s). 0: partial b&b tree, 1: splits, 2: crosses (not implemented), 3: custom.\n";
                helpstring += "-r num rounds, --rounds=num rounds\n\tNumber of rounds of cuts to apply.\n";
                helpstring += "-R num seconds, --prlp_timelimit=num seconds\n\tNumber of seconds allotted for solving the PRLP.\n";
                helpstring += "-t num seconds, --timelimit=num seconds\n\tTotal number of seconds allotted for cut generation.\n";
                helpstring += "\n# Partial branch-and-bound options #\n";
                helpstring += "-s strategy, --partial_bb_strategy=strategy\n\tPartial branch-and-bound strategy; this is a complicated parameter, and the user should check params.hpp for the description.\n";
                helpstring += "-S num strong, --partial_bb_num_strong=num strong\n\tNumber of candidates for strong branching to consider during the creation of the partial branch-and-bound tree.\n";
                helpstring += "-T num seconds, --partial_bb_timelimit=num seconds\n\tTotal number of seconds allotted for generating the partial branch-and-bound tree.\n";
                helpstring += "--tikz=[GEN_TIKZ_STRING | GEN_TIKZ_STRING_WITH_VPCS | ... ]\n\tWhether to generate code that prints tree used in partial b&b tree.\n";
                helpstring += "\n# Objective options #\n";
                helpstring += "--use_all_ones=0/1\n\tUse all ones objective.\n";
                helpstring += "--use_disj_lb=0/1\n\tUse disjunctive lower bound objective.\n";
                helpstring += "--use_iter_bilinear=num iters to do\n\tNumber of iterations to do in iterative bilinear procedure (1 = cut off the optimal post-SIC point).\n";
                helpstring += "--use_tight_points=0/1\n\tUse objectives for being tight on points in collection.\n";
                helpstring += "--use_tight_rays=0/1\n\tUse objectives for being tight on rays in collection.\n";
                helpstring += "--use_unit_vectors=0/1\n\tUse unit vectors in nonbasic space.\n";
                helpstring += "\n# Branch-and-bound options #\n";
                helpstring += "-b 0+ --bb_runs=0+\n\tNumber of branch-and-bound repeats.\n";
                helpstring += "-B strategy --bb_strategy=strategy\n\tBranch-and-bound strategy (see VPCParameters.hpp; default = 536, corresponding to gurobi: 8 + user_cuts: 16 + presolve_off: 512; another common setting is 10776, which in addition enables heuristics_off: 2048, use_best_bound: 8192). Can be set individually.\n";
                helpstring += "--cbc=0/1\n";
                helpstring += "--cplex=0/1\n";
                helpstring += "--gurobi=0/1\n";
                helpstring += "--use_best_bound=0/1\n";
                helpstring += "--user_cuts=0/1\n";
                helpstring += "--bb_all_cuts=0/1\n";
                helpstring += "--bb_gmics=0/1\n";
                helpstring += "--bb_heuristics=0/1\n";
                helpstring += "--bb_presolve=0/1\n";
                helpstring += "--bb_strong_branching=0/1\n";
                helpstring += "--bb_mode={0,1,10,11,100,...,111}\n\tWhich branch-and-bound experiments to run (ones = no cuts, tens = vpcs, hundreds = gmics).\n";
                helpstring += "--bb_timelimit=num seconds\n\tNumber of seconds allotted to the solver to run branch-and-bound tests.\n";
                helpstring += "## END OF HELP ##\n";
                std::cout << helpstring << std::endl;
                status = 1;
               }
    } // switch statement for input
  } // process args
  //for (int i = 0; i < argc; i++) {
  //  std::cout << argv[i] << " ";
  //}
  //std::cout << std::endl;
  return status;
} /* processArgs */
