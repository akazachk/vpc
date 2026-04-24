#include "main_input.hpp"
#include <cstdlib>
#include <cstring>
#include <getopt.h>
#include <iostream>
#include "utility.hpp"

using namespace VPCParametersNamespace;

namespace main_input {
int processArgs(int argc, char** argv, VPCParametersNamespace::VPCParameters& params) {
  // Handle inputs
  // struct option declared in getopt.h
  // name: name of the long option
  // has_arg: 0,1,2 for none, required, or optional
  // *flag: how results are returned; if NULL, getopt_long() returns val (e.g., can be the equivalent short option character), and o/w getopt_long() returns 0, and flag points to a var which is set to val if the option is found, but left unchanged if the option is not found
  // val: value to return, or to load into the variable pointed to by flag
  const char* const short_opts = "b:B:c:d:D:f:g:hi:l:m:o:r:R:s:S:t:T:v:";
  const struct option long_opts[] =
  {
      {"away",                  required_argument, 0, 'a'*'2'},
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
      {"max_support_rel",       required_argument, 0, 'm'*'1'},
      {"max_num_hot_start_viol",required_argument, 0, 'm'*'2'},
      {"mode",                  required_argument, 0, 'm'},
      {"optfile",               required_argument, 0, 'o'},
      {"partial_bb_strategy",   required_argument, 0, 's'},
      {"partial_bb_num_strong", required_argument, 0, 'S'},
      {"partial_bb_keep_pruned",required_argument, 0, 'S'*'1'},
      {"partial_bb_timelimit",  required_argument, 0, 'T'},
      {"preprocess",            required_argument, 0, 'p'*'1'},
      {"random_seed",           required_argument, 0, 'r'*'1'},
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
      case 'a'*'2': {
                      double val;
                      doubleParam param = doubleParam::AWAY;
                      if (!parseDouble(optarg, val)) {
                        error_msg(errorstring, "Error reading %s. Given value: %s.\n", params.name(param).c_str(), optarg);
                        exit(1);
                      }
                      params.set(param, val);
                      break;
                    }
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
      case 'm'*'1': {
                      double val;
                      doubleParam param = doubleParam::MAX_SUPPORT_REL;
                      if (!parseDouble(optarg, val)) {
                        error_msg(errorstring, "Error reading %s. Given value: %s.\n", params.name(param).c_str(), optarg);
                        exit(1);
                      }
                      params.set(param, val);
                      break;
                    }
      case 'm'*'2': {
                 int val;
                 intParam param = intParam::MAX_NUM_HOT_START_VIOL;
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
      case 'S'*'1': {
                  int val;
                  intParam param = intParam::PARTIAL_BB_KEEP_PRUNED_NODES;
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
      case 'r'*'1': {
                      int val;
                      intParam param = intParam::RANDOM_SEED;
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
                helpstring += "--disj_options={num_terms1, num_terms2,...}\n\tNumber of terms to use in each round.\n";
                helpstring += "-g -1/0/1, --gomory=-1/0/1\n\t0: do not use Gomory cuts before generating VPCs, +/-1: generate Gomory cuts before generating VPCs (-1: only gen, +1: also apply to LP).\n";
                helpstring += "-m mode, --mode=mode\n\tMode for generating disjunction(s). 0: partial b&b tree, 1: splits, 2: crosses (not implemented), 3: custom, 4: do set of disj_terms.\n";
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
                helpstring += "\n# Other options #\n";
                helpstring += "--max_num_hot_start_viol=num viols\n\tMaximum number of violations to allow before hot start is disabled and use solver resolve instead.\n";
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
} // namespace main_input
