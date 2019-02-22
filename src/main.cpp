// V-Polyhedral Disjunctive Cuts
// Author:   A. M. Kazachkov
// Date:     2018-Dec-24
//-----------------------------------------------------------------------------

#include <iostream>
#include <cstdio>
#include <string>
#include <vector>
#include <chrono> // for timing

// For option handling
#include <getopt.h> // getopt, getopt_long
#define no_argument 0
#define required_argument 1
#define optional_argument 2

// COIN-OR
#include <OsiCuts.hpp>
//#include <CoinPackedMatrix.hpp>
//#include <CglGomory.hpp>
//#include <CglGMI.hpp>
//#include <CglLandP.hpp>
//#include <OsiCuts.hpp>

// Project files
#include "CglVPC.hpp" // defines SolverInterface
#include "SolverHelper.hpp"
#include "VPCParameters.hpp"
#include "TimeStats.hpp"
#include "utility.hpp"

// Main file variables
VPCParameters params;
OsiSolverInterface* solver;
//ExitReason exitReason = ExitReason::UNKNOWN;
std::string instname;

enum OverallTimeStats {
  INITIAL_SOLVE_TIME,
  GEN_VPC_TIME,
  APPLY_VPC_TIME,
  BB_NOVPC_TIME,
  BB_VPC_TIME,
  NUM_TIME_STATS
}; /* OverallTimeStats */
const std::vector<std::string> OverallTimeStatsName {
  "INITIAL_SOLVE_TIME",
  "GEN_VPC_TIME",
  "APPLY_VPC_TIME",
  "BB_NOVPC_TIME",
  "BB_VPC_TIME"
};
TimeStats timer;

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

void processArgs(int argc, char** argv);
void initializeSolver(OsiSolverInterface* &solver);
int wrapUp(int retCode);

int main(int argc, char** argv) {
  // Do this early in your program's initialization
  std::signal(SIGABRT, signal_handler_with_error_msg);
  std::signal(SIGSEGV, signal_handler_with_error_msg);

  // Input handling
  printf("## V-Polyhedral Disjunctive Cuts ##\n");
  printf("Aleksandr M. Kazachkov\n");
  printf("Based on joint work with Egon Balas\n");
  processArgs(argc, argv);

  printf("Instance file: %s.\n", params.get(stringParam::FILENAME).c_str());
  for (int i = 0; i < argc; i++) {
    std::cout << argv[i] << " ";
  }
  std::cout << std::endl;

#ifdef TRACE
  // Print parameters
  printf("\n## Parameter values ##\n");
  printParams(params, stdout);
#endif

  // Prepare logfile
  const std::string logname = params.get(stringParam::LOGFILE);
  if (!logname.empty()) {
    const bool logexists = fexists(logname.c_str());
    params.logfile = fopen(logname.c_str(), "a");
    if (!logexists) {
      fprintf(params.logfile, "Instance,");
      fprintf(params.logfile, "\n");
      fprintf(params.logfile, "%s,", instname.c_str());
      fflush(params.logfile);
    }
  }

  // Set up timing
  for (int t = 0; t < OverallTimeStats::NUM_TIME_STATS; t++) {
    timer.register_name(OverallTimeStatsName[t]);
  }

  initializeSolver(solver);

  timer.start_timer(OverallTimeStats::INITIAL_SOLVE_TIME);
  solver->initialSolve();
  timer.end_timer(OverallTimeStats::INITIAL_SOLVE_TIME);

  timer.start_timer(OverallTimeStats::GEN_VPC_TIME);
  CglVPC gen(params);
  OsiCuts vpcs;
  gen.generateCuts(*solver, vpcs); // solution may change slightly due to enable factorization called in getProblemData...
  timer.end_timer(OverallTimeStats::GEN_VPC_TIME);

  printf("\n## Finished VPC generation using partial tree with %d disjunctive terms, generating %d cuts. ##\n", gen.num_disj_terms, vpcs.sizeCuts());
  return wrapUp(0);
} /* main */

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
  const char* const short_opts = "d:f:hl:o:R:s:S:t:T:";
  const struct option long_opts[] =
  {
    {"disj_terms", required_argument, 0, 'd'},
    {"file", required_argument, 0, 'f'},
    {"help", no_argument, 0, 'h'},
    {"logfile", required_argument, 0, 'l'},
    {"optfile", required_argument, 0, 'o'},
    {"partial_bb_strategy", required_argument, 0, 's'},
    {"partial_bb_num_strong", required_argument, 0, 'S'},
    {"partial_bb_timelimit", required_argument, 0, 'T'},
    {"prlp_timelimit", required_argument, 0, 'R'},
    {"timelimit", required_argument, 0, 't'},
    {"use_all_ones", required_argument, 0, 'u'*'1'},
    {"use_disj_lb", required_argument, 0, 'u'*'2'},
    {"use_iter_bilinear", required_argument, 0, 'u'*'3'},
    {"use_tight_points", required_argument, 0, 'u'*'4'},
    {"use_tight_rays", required_argument, 0, 'u'*'5'},
    {"use_unit_vectors", required_argument, 0, 'u'*'6'},
    {nullptr, no_argument, nullptr, 0}
  };

  int inp;
  while ((inp = getopt_long(argc, argv, short_opts, long_opts, nullptr)) != -1) {
    switch (inp) {
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
      case 'h':
      case '?':
      default: {
                 // print help
								std::string helpstring;
								helpstring += "\n## DESCRIPTION ##\n";
								helpstring += "Code for generating V-polyhedral disjunctive cuts.\n";
								helpstring += "\n## OPTIONS ##\n";
								helpstring += "-h, --help\n\tPrint this help message.\n";
								helpstring += "\n# Input/output #\n";
								helpstring += "-f file, --file=file\n\tFilename.\n";
								helpstring += "-l logfile, --logfile=logfile\n\tWhere to print log messages.\n";
								helpstring += "-o optfile, --optfile=optfile\n\tWhere to find integer optimum value information (a csv file formatted as \"instance_name,value\" on each row).\n";
								helpstring += "\n# General VPC options #\n";
								helpstring += "-d num terms, --disj_terms=num terms\n\tMaximum number of disjunctive terms to generate.\n";
								helpstring += "-R num seconds, --prlp_timelimit=num seconds\n\tNumber of seconds allotted for solving the PRLP.\n";
								helpstring += "-t num seconds, --timelimit=num seconds\n\tTotal number of seconds allotted for cut generation.\n";
								helpstring += "\n# Partial branch-and-bound options #\n";
								helpstring += "-s strategy, --partial_bb_strategy=strategy\n\tPartial branch-and-bound strategy; this is a complicated parameter, and the user should check params.hpp for the description.\n";
								helpstring += "-S num strong, --partial_bb_num_strong=num strong\n\tNumber of candidates for strong branching to consider during the creation of the partial branch-and-bound tree.\n";
								helpstring += "-T num seconds, --partial_bb_timelimit=num seconds\n\tTotal number of seconds allotted for generating the partial branch-and-bound tree.\n";
								helpstring += "\n# Objective options #\n";
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

void initializeSolver(OsiSolverInterface* &solver) {
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
  std::string in_file_ext(fullfilename.substr(found_dot + 1));
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

  // Generate cuts
  solver = new SolverInterface;
  setLPSolverParameters(solver);

  if (in_file_ext.compare("lp") == 0) {
#ifdef TRACE
    printf("\n## Reading LP file. ##\n");
#endif
    solver->readLp(fullfilename.c_str());
  } else {
    if (in_file_ext.compare("mps") == 0) {
#ifdef TRACE
      printf("\n## Reading MPS file. ##\n");
#endif
      solver->readMps(fullfilename.c_str());
    } else {
      error_msg(errorstring, "Unrecognized extension: %s.\n",
          in_file_ext.c_str());
      exit(1);
    }
  } // read file
} /* initializeSolver */

/**
 * Close the logfile and print to it
 */
int wrapUp(int retCode /*= 0*/) {
  FILE* logfile = params.logfile;
  if (logfile != NULL) {
    printParams(params, logfile);
//    const int exitReasonInt = static_cast<int>(exitReason);
//    fprintf(logfile, "%s,", ExitReasonName[exitReasonInt].c_str());
    std::time_t end_time_t;
    time(&end_time_t);
    struct tm* timeinfo = localtime(&end_time_t);
    char end_time_string[25];
    snprintf(end_time_string, sizeof(end_time_string) / sizeof(char), "%s",
        asctime(timeinfo));
    fprintf(logfile, "%s", end_time_string);
    fprintf(logfile, "\n");

    fclose(logfile); // closes params.logfile
  }
  if (solver) {
    delete solver;
  }
  return retCode;
} /* wrapUp */
