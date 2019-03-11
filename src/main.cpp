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

// For option handling
#include <getopt.h> // getopt, getopt_long
#define no_argument 0
#define required_argument 1
#define optional_argument 2

// COIN-OR
#include <OsiCuts.hpp>

// Project files
#include "analysis.hpp" // analyzeStrength, analalyzeBB
#include "BBHelper.hpp" // runBBTests
#include "CglVPC.hpp"
#include "CutHelper.hpp"
#include "Disjunction.hpp" // needed to access disjunction properties
#include "DisjunctionHelper.hpp" // for custom disjunctions
#include "SolverHelper.hpp" // defines SolverInterface
#include "VPCParameters.hpp"
#include "TimeStats.hpp"
#include "utility.hpp"

enum OverallTimeStats {
  TOTAL_TIME,
  INIT_SOLVE_TIME,
  VPC_GEN_TIME,
  VPC_APPLY_TIME,
  BB_TIME,
  NUM_TIME_STATS
}; /* OverallTimeStats */
const std::vector<std::string> OverallTimeStatsName {
  "TOTAL_TIME",
  "INIT_SOLVE_TIME",
  "VPC_GEN_TIME",
  "VPC_APPLY_TIME",
  "BB_TIME"
};

// Main file variables
VPCParameters params;
OsiSolverInterface* solver, *origSolver;
OsiCuts gmics, vpcs;
std::string instname = "", in_file_ext = "";
ExitReason exitReason;
TimeStats timer;
std::time_t start_time_t, end_time_t;

SummaryBoundInfo boundInfo;
SummaryBBInfo info_nocuts, info_mycuts, info_allcuts;
SummaryDisjunctionInfo disjInfo;
std::vector<SummaryCutInfo> cutInfoVec;
SummaryCutInfo cutInfo, cutInfoGMICs;

// For output
std::string cut_output = "", bb_output = "";

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

void startUp(int argc, char** argv);
void processArgs(int argc, char** argv);
void initializeSolver(OsiSolverInterface* &solver);
void doCustomRoundOfCuts(int round_ind, OsiCuts& vpcs, CglVPC& gen, int& num_disj);
int wrapUp(int retCode);

/****************** MAIN FUNCTION **********************/
int main(int argc, char** argv) {
  // Do this early in your program's initialization
  std::signal(SIGABRT, signal_handler_with_error_msg);
  std::signal(SIGSEGV, signal_handler_with_error_msg);

  // Set up timing
  for (int t = 0; t < OverallTimeStats::NUM_TIME_STATS; t++) {
    timer.register_name(OverallTimeStatsName[t]);
  }

  // Print welcome message, set up logfile
  timer.start_timer(OverallTimeStats::TOTAL_TIME);
  startUp(argc, argv);

  // Set up solver and get initial solution
  initializeSolver(solver);
  timer.start_timer(OverallTimeStats::INIT_SOLVE_TIME);
  solver->initialSolve();
  checkSolverOptimality(solver, false);
  timer.end_timer(OverallTimeStats::INIT_SOLVE_TIME);
  boundInfo.lp_obj = solver->getObjValue();

  // Save original solver in case we wish to come back to it later
  origSolver = solver->clone();
  if (!origSolver->isProvenOptimal()) {
    origSolver->initialSolve();
    checkSolverOptimality(origSolver, false);
  }

  // Possibly preprocess instead of doing cuts
  if (params.get(TEMP) == static_cast<int>(TempOptions::PREPROCESS)) {
    // Cleaning involves running Gurobi presolve
//    performCleaning(solver, params, CLEANING_MODE_OPTION);

    printf("\n## Finished cleaning. ##\n");
    return wrapUp(0);
  }

  // Now do rounds of cuts, until a limit is reached (e.g., time, number failures, number cuts, or all rounds are exhausted)
  boundInfo.num_vpc = 0, boundInfo.num_gmic = 0;
  int num_rounds = params.get(ROUNDS);
  std::vector<OsiCuts> vpcs_by_round(num_rounds);
  cutInfoVec.resize(num_rounds);
  int round_ind = 0;
  int num_disj = 0;
  for (round_ind = 0; round_ind < num_rounds; ++round_ind) {
    if (num_rounds > 1) {
      printf("\n## Starting round %d/%d. ##\n", round_ind+1, num_rounds);
    }
    timer.start_timer(OverallTimeStats::VPC_GEN_TIME);
    CglVPC gen(params);

    // Store the initial solve time in order to set a baseline for the PRLP resolve time
    gen.timer.add_value(
        CglVPC::VPCTimeStatsName[static_cast<int>(CglVPC::VPCTimeStats::INIT_SOLVE_TIME)],
        timer.get_value(OverallTimeStats::INIT_SOLVE_TIME));

    // Proceed with custom disjunctions if specified; otherwise, the disjunction will be set up in the CglVPC class
    if (params.get(MODE) != static_cast<int>(CglVPC::VPCMode::CUSTOM)) {
      gen.generateCuts(*solver, vpcs_by_round[round_ind]); // solution may change slightly due to enable factorization called in getProblemData...
      exitReason = gen.exitReason;
      updateDisjInfo(disjInfo, num_disj, gen);
      updateCutInfo(cutInfoVec[round_ind], gen);
      if (gen.disj()) {
        num_disj++;
        boundInfo.num_vpc += gen.num_cuts;
        if (boundInfo.best_disj_obj < gen.disj()->best_obj)
          boundInfo.best_disj_obj = gen.disj()->best_obj;
        if (boundInfo.worst_disj_obj < gen.disj()->worst_obj)
          boundInfo.worst_disj_obj = gen.disj()->worst_obj;
        boundInfo.ip_obj = gen.ip_obj;
      }
    } // check if mode is _not_ CUSTOM
    else {
      doCustomRoundOfCuts(round_ind, vpcs_by_round[round_ind], gen, num_disj);
    } // check if mode is CUSTOM
    timer.end_timer(OverallTimeStats::VPC_GEN_TIME);

    timer.start_timer(OverallTimeStats::VPC_APPLY_TIME);
    applyCutsCustom(solver, vpcs_by_round[round_ind]);
    boundInfo.vpc_obj = solver->getObjValue();
    timer.end_timer(OverallTimeStats::VPC_APPLY_TIME);

    vpcs.insert(vpcs_by_round[round_ind]);

    printf(
        "\n## Round %d/%d: Completed round of VPC generation (exit reason: %s). # cuts generated = %d.\n",
        round_ind + 1, params.get(ROUNDS),
        ExitReasonName[static_cast<int>(exitReason)].c_str(),
        vpcs_by_round[round_ind].sizeCuts());
    fflush(stdout);
    printf("Initial obj value: %1.6f. New obj value: %s. Disj lb: %s. ##\n",
        boundInfo.lp_obj, stringValue(solver->getObjValue(), "%1.6f").c_str(),
        stringValue(boundInfo.best_disj_obj, "%1.6f").c_str());

    // Exit early from rounds of cuts if no cuts generated or solver is not optimal
    if (gen.num_cuts == 0 || !solver->isProvenOptimal()
        || isInfinity(std::abs(solver->getObjValue())))
      break;
  } // loop over rounds of cuts
  if (round_ind < num_rounds)
    num_rounds = round_ind+1;

  // Do branch-and-bound experiments (if requested)
  if (params.get(BB_RUNS) != 0) {
    // Collect cuts from all rounds
    timer.start_timer(BB_TIME);
    runBBTests(params, info_nocuts, info_mycuts, info_allcuts,
        params.get(stringParam::FILENAME), solver, boundInfo.ip_obj, &vpcs, NULL);
    timer.end_timer(BB_TIME);
  }

  printf(
      "\n## Finished VPC generation with %d cuts. Initial obj value: %s. Final obj value: %s. Disj lb: %s. ##\n",
      boundInfo.num_vpc,
      stringValue(boundInfo.lp_obj, "%1.6f").c_str(),
      stringValue(boundInfo.vpc_obj, "%1.6f").c_str(),
      stringValue(boundInfo.best_disj_obj, "%1.6f").c_str());
  timer.end_timer(OverallTimeStats::TOTAL_TIME);

  // Do analyses in preparation for printing
  setCutInfo(cutInfo, num_rounds, cutInfoVec.data());
  analyzeStrength(params, cutInfoGMICs, cutInfo, solver, &gmics, &vpcs,
      boundInfo, cut_output);
  analyzeBB(params, info_nocuts, info_mycuts, info_allcuts, bb_output);
  return wrapUp(0);
} /* main */

/**
 * Call this early to print welcome message, etc.
 */
void startUp(int argc, char** argv) {
  // Input handling
  printf("## V-Polyhedral Disjunctive Cuts ##\n");
  printf("Aleksandr M. Kazachkov\n");
  printf("Based on joint work with Egon Balas\n");
  for (int i = 0; i < argc; i++) {
    std::cout << argv[i] << " ";
  }
  std::cout << std::endl;

  time(&start_time_t);
  struct tm* timeinfo = localtime(&start_time_t);
  char time_string[25];
  snprintf(time_string, sizeof(time_string) / sizeof(char), "%s", asctime(timeinfo));
  printf("Start time: %s\n", time_string);

  processArgs(argc, argv);

  // Get instance file
  printf("Instance file: %s.\n", params.get(stringParam::FILENAME).c_str());

  std::string fullfilename = params.get(stringParam::FILENAME);
  // Get file name stub
  size_t found_dot = fullfilename.find_last_of(".");
  std::string filename = fullfilename.substr(0, found_dot);

  // Put string after last '.' into string in_file_ext
  if (found_dot >= fullfilename.length()) {
    error_msg(errorstring, "Cannot find the file extension (no '.' in input file name: %s).\n", fullfilename.c_str());
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }

  // Check if archived file
  in_file_ext = fullfilename.substr(found_dot + 1);
  if (in_file_ext.compare("gz") == 0 || in_file_ext.compare("bz2") == 0) {
    unsigned found_dot_tmp = filename.find_last_of('.');

    // Put string after last '.' into string in_file_ext
    if (found_dot_tmp >= filename.length()) {
      error_msg(errorstring,
          "Other than gz or bz2, cannot find the file extension (no '.' in input file name: %s).\n",
          fullfilename.c_str());
      writeErrorToLog(errorstring, params.logfile);
      exit(1);
    }

    in_file_ext = filename.substr(found_dot_tmp + 1);
    filename = filename.substr(0, found_dot_tmp);
  }

//  size_t slashindex = fullfilename.find_last_of("/\\");
//  const std::string dir = (slashindex != std::string::npos) ? fullfilename.substr(0, slashindex+1) : "./";
//  const std::string filestub = (slashindex != std::string::npos) ? fullfilename.substr(slashindex+1) : fullfilename;
  size_t slashindex = filename.find_last_of("/\\");
  instname = (slashindex != std::string::npos) ? filename.substr(slashindex+1) : filename;

  // Prepare logfile
  const std::string logname = params.get(stringParam::LOGFILE);
  if (!logname.empty()) {
    const bool logexists = fexists(logname.c_str());
    params.logfile = fopen(logname.c_str(), "a");
    if (!logexists) {
      printHeader(params, OverallTimeStatsName);
    }
    // Print instance name and parameters
    fprintf(params.logfile, "%s,", instname.c_str());
    printParams(params, params.logfile, 2); // only values
    fflush(params.logfile);
  }
} /* startUp */

/**
 * Close the logfile and print to it
 */
int wrapUp(int retCode /*= 0*/) {
  const int exitReasonInt = static_cast<int>(exitReason);

  time(&end_time_t);
  struct tm* start_timeinfo = localtime(&start_time_t);
  struct tm* end_timeinfo = localtime(&end_time_t);
  char start_time_string[25], end_time_string[25];
  snprintf(start_time_string, sizeof(start_time_string) / sizeof(char), "%s", asctime(start_timeinfo));
  snprintf(end_time_string, sizeof(end_time_string) / sizeof(char), "%s", asctime(end_timeinfo));

  FILE* logfile = params.logfile;
  if ((params.get(TEMP) != static_cast<int>(TempOptions::PREPROCESS))
      && (logfile != NULL)) {
    // Bound and gap info
    printBoundAndGapInfo(boundInfo, params.logfile);
    // B&B info
    printSummaryBBInfo({info_nocuts, info_mycuts}, params.logfile);
    // Orig prob
    printOrigProbInfo(origSolver, params.logfile);
    // Post-cut prob
    printPostCutProbInfo(solver, cutInfoGMICs, cutInfo, params.logfile);
    // Disj info
    printDisjInfo(disjInfo, params.logfile);
    // Cut, obj, fail info
    printCutInfo(cutInfoGMICs, cutInfo, params.logfile);
    // Full B&B info
    printFullBBInfo({info_nocuts, info_mycuts}, params.logfile);
    // Print time info
    timer.print(params.logfile, 2); // only values
    // Print exit reason and finish
    fprintf(logfile, "%s,", ExitReasonName[exitReasonInt].c_str());
    fprintf(logfile, "%s,", end_time_string);
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

  printf("\n## Exiting VPC generation with reason %s. ##\n", ExitReasonName[exitReasonInt].c_str());
  printf("Instance: %s\n", instname.c_str());
  printf("Start time: %s\n", start_time_string);
  printf("End time: %s\n", end_time_string);
  printf("Elapsed time: %.f seconds\n", difftime(end_time_t, start_time_t));

  if (solver) {
    delete solver;
  }
  if (origSolver) {
    delete origSolver;
  }
  return retCode;
} /* wrapUp */

void initializeSolver(OsiSolverInterface* &solver) {
  // Generate cuts
  solver = new SolverInterface;
  setLPSolverParameters(solver, params.get(VERBOSITY));

  if (in_file_ext.compare("lp") == 0) {
#ifdef TRACE
    printf("\n## Reading LP file. ##\n");
#endif
    solver->readLp(params.get(stringParam::FILENAME).c_str());
  } else {
    if (in_file_ext.compare("mps") == 0) {
#ifdef TRACE
      printf("\n## Reading MPS file. ##\n");
#endif
      solver->readMps(params.get(stringParam::FILENAME).c_str());
    } else {
      error_msg(errorstring, "Unrecognized extension: %s.\n",
          in_file_ext.c_str());
      exit(1);
    }
  } // read file
} /* initializeSolver */

void doCustomRoundOfCuts(int round_ind, OsiCuts& vpcs, CglVPC& gen, int& num_disj) {
  std::vector<Disjunction*> disjVec;
  printf("\n## Setting up disjunction(s) ##\n");
  gen.timer.start_timer(CglVPC::VPCTimeStatsName[static_cast<int>(CglVPC::VPCTimeStats::DISJ_GEN_TIME)]);
  ExitReason setDisjExitReason = setDisjunctions(disjVec, solver, params, CglVPC::VPCMode::SPLITS);
  gen.timer.end_timer(CglVPC::VPCTimeStatsName[static_cast<int>(CglVPC::VPCTimeStats::DISJ_GEN_TIME)]);
  const int numDisj = disjVec.size();

  // If integer-optimal solution was found, all disjunctions but one will have been deleted
  if (setDisjExitReason == ExitReason::PARTIAL_BB_OPTIMAL_SOLUTION_FOUND_EXIT) {
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
    exitReason = ExitReason::PARTIAL_BB_OPTIMAL_SOLUTION_FOUND_EXIT;
    boundInfo.num_vpc += gen.num_cuts;
    updateCutInfo(cutInfoVec[round_ind], gen);
  } // check if integer-optimal solution
  else if (setDisjExitReason == ExitReason::SUCCESS_EXIT && numDisj > 0) {
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
      boundInfo.ip_obj = gen.ip_obj;
      updateDisjInfo(disjInfo, num_disj, gen);
      updateCutInfo(cutInfoVec[round_ind], gen);
    }
    gen.setupRepeatedUse(false);
  }  // if successful generation of disjunction, generate cuts
  else if (setDisjExitReason == ExitReason::NO_DISJUNCTION_EXIT) {
    // Do nothing
    exitReason = ExitReason::NO_DISJUNCTION_EXIT;
  }
  else {
    error_msg(errorstring, "Unknown exit reason (%s) from setDisjunctions.\n", ExitReasonName[static_cast<int>(setDisjExitReason)].c_str());
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
void processArgs(int argc, char** argv) {
  // Handle inputs
  // struct option declared in getopt.h
  // name: name of the long option
  // has_arg: 0,1,2 for none, required, or optional
  // *flag: how results are returned; if NULL, getopt_long() returns val (e.g., can be the equivalent short option character), and o/w getopt_long() returns 0, and flag points to a var which is set to val if the option is found, but left unchanged if the option is not found
  // val: value to return, or to load into the variable pointed to by flag
  const char* const short_opts = "b:B:c:d:f:hl:m:o:r:R:s:S:t:T:v:";
  const struct option long_opts[] =
  {
      {"bb_runs", required_argument, 0, 'b'},
      {"bb_mode", required_argument, 0, 'b'*'2'},
      {"bb_strategy", required_argument, 0, 'B'},
      {"cutlimit", required_argument, 0, 'c'},
      {"disj_terms", required_argument, 0, 'd'},
      {"file", required_argument, 0, 'f'},
      {"help", no_argument, 0, 'h'},
      {"logfile", required_argument, 0, 'l'},
      {"mode", required_argument, 0, 'm'},
      {"optfile", required_argument, 0, 'o'},
      {"partial_bb_strategy", required_argument, 0, 's'},
      {"partial_bb_num_strong", required_argument, 0, 'S'},
      {"partial_bb_timelimit", required_argument, 0, 'T'},
      {"rounds", required_argument, 0, 'r'},
      {"prlp_timelimit", required_argument, 0, 'R'},
      {"timelimit", required_argument, 0, 't'},
      {"use_all_ones", required_argument, 0, 'u'*'1'},
      {"use_disj_lb", required_argument, 0, 'u'*'2'},
      {"use_iter_bilinear", required_argument, 0, 'u'*'3'},
      {"use_tight_points", required_argument, 0, 'u'*'4'},
      {"use_tight_rays", required_argument, 0, 'u'*'5'},
      {"use_unit_vectors", required_argument, 0, 'u'*'6'},
      {"verbosity", required_argument, 0, 'v'},
      {nullptr, no_argument, nullptr, 0}
  };

  int inp;
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
      case 'b'*'2': {
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
      case 'f': {
                  params.set(stringParam::FILENAME, optarg);
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
                helpstring += "-t, --temp\n\tSet temporary options (e.g., value of 1 = do preprocessing on instance).\n";
                helpstring += "\n# Input/output #\n";
                helpstring += "-f file, --file=file\n\tFilename.\n";
                helpstring += "-l logfile, --logfile=logfile\n\tWhere to print log messages.\n";
                helpstring += "-o optfile, --optfile=optfile\n\tWhere to find integer optimum value information (a csv file formatted as \"instance_name,value\" on each row).\n";
                helpstring += "-v level, --verbosity=level\n\tVerbosity level (0: print little, 1: let solver output be visible).\n";
                helpstring += "\n# General VPC options #\n";
                helpstring += "-c num cuts, --cutlimit=num cuts\n\tMaximum number of cuts to generate (0 = no limit, -k = k * # fractional variables at root).\n";
                helpstring += "-d num terms, --disj_terms=num terms\n\tMaximum number of disjunctive terms or disjunctions to generate (depending on mode).\n";
                helpstring += "-m mode, --mode=mode\n\tMode for generating disjunction(s). 0: partial b&b tree, 1: splits, 2: crosses (not implemented), 3: custom.\n";
                helpstring += "-r num rounds, --rounds=num rounds\n\tNumber of rounds of cuts to apply.\n";
                helpstring += "-R num seconds, --prlp_timelimit=num seconds\n\tNumber of seconds allotted for solving the PRLP.\n";
                helpstring += "-t num seconds, --timelimit=num seconds\n\tTotal number of seconds allotted for cut generation.\n";
                helpstring += "\n# Partial branch-and-bound options #\n";
                helpstring += "-s strategy, --partial_bb_strategy=strategy\n\tPartial branch-and-bound strategy; this is a complicated parameter, and the user should check params.hpp for the description.\n";
                helpstring += "-S num strong, --partial_bb_num_strong=num strong\n\tNumber of candidates for strong branching to consider during the creation of the partial branch-and-bound tree.\n";
                helpstring += "-T num seconds, --partial_bb_timelimit=num seconds\n\tTotal number of seconds allotted for generating the partial branch-and-bound tree.\n";
                helpstring += "\n# Objective options #\n";
                helpstring += "--use_all_ones=0/1\n\tUse all ones objective.\n";
                helpstring += "--use_disj_lb=0/1\n\tUse disjunctive lower bound objective.\n";
                helpstring += "--use_iter_bilinear=num iters to do\n\tNumber of iterations to do in iterative bilinear procedure (1 = cut off the optimal post-SIC point).\n";
                helpstring += "--use_tight_points=0/1\n\tUse objectives for being tight on points in collection.\n";
                helpstring += "--use_tight_rays=0/1\n\tUse objectives for being tight on rays in collection.\n";
                helpstring += "--use_unit_vectors=0/1\n\tUse unit vectors in nonbasic space.\n";
                helpstring += "\n# Branch-and-bound options #\n";
                helpstring += "-b 0+ --bb_runs=0+\n\tNumber of branch-and-bound repeats.\n";
                helpstring += "-B strategy --bb_strategy=strategy\n\tBranch-and-bound strategy (see BBHelper.hpp).\n";
                helpstring += "--bb_mode={0,1,10,11,100,...,111}\n\tWhich branch-and-bound experiments to run (ones = no cuts, tens = vpcs, hundreds = gmics).\n";
                helpstring += "## END OF HELP ##\n";
                std::cout << helpstring << std::endl;
                exit(1);
               }
    } // switch statement for input
  } // process args
  //for (int i = 0; i < argc; i++) {
  //  std::cout << argv[i] << " ";
  //}
  //std::cout << std::endl;
} /* processArgs */
