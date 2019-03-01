// Name:     VPCParameters.hpp
// Author:   A. M. Kazachkov
// Date:     2019-Feb-20
//-----------------------------------------------------------------------------
#pragma once

/********************************************************************************************************************
 * This file contains the parameters/constants that are used in the CglVPC class
 * (also in the PRLP and Disjunction classes)
 *
 * To add a new parameter/constant:
 * 1. Add it to the relevant enum (intParam, doubleParam, stringParam, intConst, doubleConst) or create a new enum
 * 2. Add the name of the parameter/constant (in the same position) in *ParamName
 * 3. In the struct VPCParameters, add the default value of the parameter/constant in the relevant place
 * 4. Optionally, add a way to set the parameter in the option handling part of the code
 ********************************************************************************************************************/

#include <map>
#include <string>
#include <vector>
#include <cstdio> // fprintf
#include <fstream>
#include <sstream>

#include "utility.hpp" // for parseInt/Double, stringValue, and lowerCaseString

// Define SolverInterface that we can later change in a way that will be used in all VPC files
#ifdef USE_CLP
  #include <OsiClpSolverInterface.hpp>
  using SolverInterface = OsiClpSolverInterface;
#else
  #include <OsiSolverInterface.hpp>
  using SolverInterface = OsiSolverInterface;
#endif

/********** PARAMETERS **********/
enum intParam {
  CUTLIMIT, // max number of cuts generated; 0 = no limit, -k = k * # fractional variables at root
  DISJ_TERMS, // number of disjunctive terms or number of disjunctions, depending on MODE
  MODE, // 0: partial b&b tree, 1: splits, 2: crosses (not implemented), 3: custom
  // PARTIAL_BB_STRATEGY:
  // Total used to decide the choose:
  // variable decision => hundreds digit: 0: default, 1: default+second criterion, 2: max min change+second (max max change), 3: second-best default, 4: second-best max-min change, -x: -1 * (1+x);
  // branch decision => tens digit: 0: default, 1: dynamic, 2: strong, 3: none;
  // node comparison decision => ones digit: 0: default: 1: bfs, 2: depth, 3: estimate, 4: objective
  PARTIAL_BB_STRATEGY,
  PARTIAL_BB_NUM_STRONG, // -1: num cols, -2: sqrt(num cols), >= 0: that many
  PRLP_FLIP_BETA, // controls rhs in nb space, -1: do not cut away LP opt, 0: cut away LP opt, 1: both
  ROUNDS, // number of VPC rounds to do
  STRENGTHEN, // 0: no, 1: yes, when possible, 2: same as 1 plus add GMICs to strengthen each disjunctive term
  TEMP, // useful for various temporary parameter changes
  // Objective options
  USE_ALL_ONES, // 0: do not use, 1: use
  USE_DISJ_LB, // 0: do not use, 1: use
  USE_ITER_BILINEAR, // 0: do not use, 1+: number of iterations to do
  USE_TIGHT_POINTS, // 0: do not use, 1+: num points to try
  USE_TIGHT_RAYS, // 0: do not use, 1+: num rays to try, -1: try the first sqrt(# rays) rays
  USE_UNIT_VECTORS, // 0: do not use, 1+: num to try, <0: abs(val) * sqrt(n)
  // Other options
  VERBOSITY,
  NUM_INT_PARAMS
}; /* intParam */
const std::vector<std::string> intParamName {
  "CUTLIMIT",
  "DISJ_TERMS",
  "MODE",
  "PARTIAL_BB_STRATEGY",
  "PARTIAL_BB_NUM_STRONG",
  "PRLP_FLIP_BETA",
  "ROUNDS",
  "STRENGTHEN",
  "TEMP",
  "USE_ALL_ONES",
  "USE_DISJ_LB",
  "USE_ITER_BILINEAR",
  "USE_TIGHT_POINTS",
  "USE_TIGHT_RAYS",
  "USE_UNIT_VECTORS",
  "VERBOSITY"
}; /* intParamName */

enum doubleParam {
  EPS,
  MIN_ORTHOGONALITY, // minimum orthogonality between cuts added to the collection
  PARTIAL_BB_TIMELIMIT,
  PRLP_TIMELIMIT, // -1 = unlimited time for initial solve, and then half that time subsequently; -2 = unlimited time always
  TIMELIMIT,
  NUM_DOUBLE_PARAMS
}; /* doubleParam */
const std::vector<std::string> doubleParamName {
  "EPS",
  "MIN_ORTHOGONALITY",
  "PARTIAL_BB_TIMELIMIT",
  "PRLP_TIMELIMIT",
  "TIMELIMIT"
}; /* doubleParamName */

enum stringParam {
  FILENAME,
  LOGFILE,
  OPTFILE,
  NUM_STRING_PARAMS
}; /* stringParam */
const std::vector<std::string> stringParamName {
  "FILE",
  "LOGFILE",
  "OPTFILE"
}; /* stringParamName */

/********** CONSTANTS **********/
enum class intConst {
  CHECK_DUPLICATES, // do not add duplicate cuts
  LUB, // value for var upper bound considered "large"
  MAX_SUPPORT_ABS,
  // VPC objective-related constants
  MODE_OBJ_PER_POINT,
  NUM_OBJ_PER_POINT, // # cuts to try to generate for the strong branching lb point (and others)
  NB_SPACE, // whether to generate cuts in the nonbasic space (currently must be set to true)
  PRLP_PRESOLVE, // 0: no presolve, 1: only initial solve, 2: both initial solve and resolve
  RANDOM_SEED,
  NUM_INT_CONST
}; /* intConst */
const std::vector<std::string> intConstName {
  "CHECK_DUPLICATES",
  "LUB",
  "MAX_SUPPORT_ABS",
  // VPC objective-related constants
  // For each point, MODE_OBJ_PER_POINT says which objectives to try
  // 0 = all of the rows that are not tight, and subtract as they get tight
  // 1 = one point/ray at time
  // 2 = keep trying even if the point/ray has become tight in the process?
  // 0 is more expensive at each step and not clear that it would be better
  // 0x: only rays
  // 1x: points+rays
  // 2x: points+rays+variables
  // 0xx: small to large angle with obj (ascending)
  // 1xx: large to small angle with obj (descending)
  // 2xx: small to large slack (ascending)
  // 3xx: large to small slack (descending)
  "MODE_OBJ_PER_POINT",
  "NUM_OBJ_PER_POINT",
  "NB_SPACE",
  "PRLP_PRESOLVE",
  "RANDOM_SEED"
}; /* intConstName */
enum class doubleConst {
  AWAY,
  DIFFEPS, // to check whether something is different enough to throw an error
  INF, // infinity (INFINITY is taken as a macro from math header)
  RAYEPS, // value for which a ray coefficient will be treated as zero
  // Time limits
  BB_TIMELIMIT, // time limit for doing branch-and-bound
  MIN_PRLP_TIMELIMIT, // minimum amount of time allotted for solving/resolving PRLP
  // Safety related constants:
  EPS_COEFF, // any cut coefficient smaller than this will be replaced by zero
  EPS_COEFF_LUB, // for variables with large upper bound, any cut coefficient smaller than this will be replaced by zero
  MIN_VIOL_ABS,
  MIN_VIOL_REL,
  MAX_DYN, // |alpha_max| / |alpha_min| upper bound; Maximum ratio between largest and smallest non zero coefficients in a cut
  MAX_DYN_LUB, // Same as MAX_DYN but when some of the variables involved in the cut have a large upper bound; should be >= MAX_DYN logically
  MAX_SUPPORT_REL,
  NUM_DOUBLE_CONST
}; /* doubleConst */
const std::vector<std::string> doubleConstName {
  "AWAY",
  "DIFFEPS",
  "INF",
  "RAYEPS",
  "BB_TIMELIMIT",
  "MIN_PRLP_TIMELIMIT",
  "EPS_COEFF",
  "EPS_COEFF_LUB",
  "MIN_VIOL_ABS",
  "MIN_VIOL_REL",
  "MAX_DYN",
  "MAX_DYN_LUB",
  "MAX_SUPPORT_REL"
}; /* doubleConstName */

/********** VPC PARAMETERS STRUCT **********/
struct VPCParameters {
  FILE* logfile = NULL; // NB: right now this is a shallow copy if this struct gets copied

  // Mutable parameters (of int, double, and string types)
  std::map<intParam, int> intParamValues {
    {intParam::CUTLIMIT, -1}, // -1 = limit is set as number of fractional integer variables at root
    {intParam::DISJ_TERMS, 0}, // no disjunction (=> no cuts)
    {intParam::MODE, 0}, // disjunction from a partial b&b tree
    {intParam::PARTIAL_BB_STRATEGY, 4}, // 004 => default var & branch decisions, and choose next node by min objective
    {intParam::PARTIAL_BB_NUM_STRONG, 5}, // consider 5 strong branching candidates
    {intParam::PRLP_FLIP_BETA, 0}, // cut away the LP optimum
    {intParam::ROUNDS, 1}, // number VPC rounds to do
    {intParam::STRENGTHEN, 1}, // strengthen GMICs but not VPCs
    {intParam::TEMP, 0}, // do not enable any temporary options
    {intParam::USE_ALL_ONES, 1},
    {intParam::USE_DISJ_LB, 1},
    {intParam::USE_ITER_BILINEAR, 0},
    {intParam::USE_TIGHT_POINTS, 0},
    {intParam::USE_TIGHT_RAYS, 0},
    {intParam::USE_UNIT_VECTORS, 0},
#ifdef TRACE
    {intParam::VERBOSITY, 1},
#else
    {intParam::VERBOSITY, 0},
#endif
  }; /* intParamValues */
  std::map<doubleParam, double> doubleParamValues {
    {doubleParam::EPS, 1e-7},
    {doubleParam::MIN_ORTHOGONALITY, 0.000},
    {doubleParam::PARTIAL_BB_TIMELIMIT, 3600},
    {doubleParam::PRLP_TIMELIMIT, -1.},
    {doubleParam::TIMELIMIT, 60}
  }; /* doubleParamValues */
  std::map<stringParam, std::string> stringParamValues {
    {stringParam::FILENAME, ""},
    {stringParam::LOGFILE, ""},
    {stringParam::OPTFILE, ""}
  }; /* stringParamValues */

  // Constants (amk: though I am not setting to const in case user wishes to change)
  std::map<intConst, double> intConstValues {
    {intConst::CHECK_DUPLICATES, 1},
    {intConst::LUB, 1e3},
    {intConst::MAX_SUPPORT_ABS, std::numeric_limits<int>::max()},
    {intConst::MODE_OBJ_PER_POINT, 121}, // one ray at a time, points+rays+variables, large to small angle (descending)
    {intConst::NUM_OBJ_PER_POINT, -2}, // sqrt(n)
    {intConst::NB_SPACE, 1}, // currently only works with true
    {intConst::PRLP_PRESOLVE, 2}, // use presolve when solving PRLP (either initial or resolve)
    {intConst::RANDOM_SEED, 628}
  }; /* intConstValues */
  std::map<doubleConst, double> doubleConstValues {
    {doubleConst::AWAY, 1e-3},
    {doubleConst::DIFFEPS, 1e-3}, // to check whether something is different enough to throw an error
    {doubleConst::INF, std::numeric_limits<double>::max()},
    {doubleConst::RAYEPS, 1e-7},
//    {doubleConst::BB_STRATEGY, 0}, // see BBHelper.hpp; 10776 = 010101000011000 => gurobi: 1, user_cuts: 1, presolve_off: 1, heuristics_off: 1, use_best_bound: 1
    {doubleConst::BB_TIMELIMIT, 3600.},
    {doubleConst::MIN_PRLP_TIMELIMIT, 5.},
    {doubleConst::EPS_COEFF, 1e-5},
    {doubleConst::EPS_COEFF_LUB, 1e-13},
    {doubleConst::MIN_VIOL_ABS, 1e-7},
    {doubleConst::MIN_VIOL_REL, 1e-7},
    {doubleConst::MAX_DYN, 1e8},
    {doubleConst::MAX_DYN_LUB, 1e13},
    {doubleConst::MAX_SUPPORT_REL, 0.9}
  }; /* doubleConstValues */

  // Methods (name/get/set)
  std::string name(intParam param) const { return intParamName[static_cast<int>(param)]; }
  std::string name(doubleParam param) const { return doubleParamName[static_cast<int>(param)]; }
  std::string name(stringParam param) const { return stringParamName[static_cast<int>(param)]; }
  std::string name(intConst param) const { return intConstName[static_cast<int>(param)]; }
  std::string name(doubleConst param) const { return doubleConstName[static_cast<int>(param)]; }

  int         get(intParam param) const { return intParamValues.find(param)->second; }
  double      get(doubleParam param) const { return doubleParamValues.find(param)->second; }
  std::string get(stringParam param) const { return stringParamValues.find(param)->second; }
  int         get(intConst param) const { return intConstValues.find(param)->second; }
  double      get(doubleConst param) const { return doubleConstValues.find(param)->second; }

  void        set(intParam param, int value) { intParamValues[param] = value; }
  void        set(doubleParam param, double value) { doubleParamValues[param] = value; }
  void        set(stringParam param, std::string value) { stringParamValues[param] = value; }

  /** Set parameters by name */
  bool        set(std::string tmpname, const std::string tmp) {
    std::string name = upperCaseStringNoUnderscore(tmpname);
    for (unsigned i = 0; i < intParamName.size(); i++) {
      std::string str1 = upperCaseStringNoUnderscore(intParamName[i]);
      if (str1.compare(name) == 0) {
        int value;
        if (tmp.compare("+inf") == 0) {
          value = std::numeric_limits<int>::max();
        } else if (tmp.compare("-inf") == 0) {
          value = std::numeric_limits<int>::min();
        } else {
          value = parseInt(tmp.c_str(), value);
        }
        set(static_cast<intParam>(i), value);
        return true;
      }
    }
    for (unsigned i = 0; i < doubleParamName.size(); i++) {
      std::string str1 = upperCaseStringNoUnderscore(doubleParamName[i]);
      if (str1.compare(name) == 0) {
        double value;
        if (tmp.compare("+inf") == 0) {
          value = std::numeric_limits<double>::max();
        } else if (tmp.compare("-inf") == 0) {
          value = std::numeric_limits<double>::lowest();
        } else {
          value = parseDouble(tmp.c_str(), value);
        }
        set(static_cast<doubleParam>(i), value);
        return true;
      }
    }
    for (unsigned i = 0; i < stringParamName.size(); i++) {
      std::string str1 = upperCaseStringNoUnderscore(stringParamName[i]);
      if (str1.compare(name) == 0) {
        set(static_cast<stringParam>(i), tmp);
        return true;
      }
    }
    return false;
  } /* set parameters by name */
}; /* struct VPCParameters */

/**
 * Print parameters and constants
 */
inline void readParams(VPCParameters& params, std::string infilename) {
  std::ifstream infile(infilename.c_str());
  if (infile.is_open()) {
    std::string line;
    while (std::getline(infile, line)) {
      std::istringstream iss(line);
      if (line.empty()) {
        continue;
      }
      std::string param_name;
      if (!(std::getline(iss, param_name, ','))) {
        warning_msg(errorstring, "Could not read parameter name. String is %s.\n", line.c_str());
        continue;
      }

      try {
        std::string token;
        if (!(std::getline(iss, token, ','))) {
          throw;
        }
        if (!params.set(param_name, token)) {
          warning_msg(warnstring,
              "Unable to find parameter %s. Value not set.\n",
              param_name.c_str());
          continue;
        }
      } catch (std::exception& e) {
        warning_msg(errorstring, "Could not read parameter value. String is %s.\n", line.c_str());
        continue;
      }
    }
    infile.close();
  } else {
    // If we were not able to open the file, throw an error
    error_msg(errorstring, "Not able to open params file with name %s.\n", infilename.c_str());
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }
} /* readParams */

/**
 * Print parameters and constants
 * amountToPrint: 0 = all (also adds a newline after each param/constant), 1 = only names, 2 = only values
 */
inline void printParams(VPCParameters& params, FILE* logfile = stdout, const int amountToPrint = 0) {
  if (!logfile)
    return;
  for (int i = 0; i < intParam::NUM_INT_PARAMS; i++) {
    intParam param = static_cast<intParam>(i);
    switch (amountToPrint) {
      case 1: {
        fprintf(logfile, "%s,", lowerCaseString(params.name(param)).c_str());
        break;
      }
      case 2: {
        fprintf(logfile, "%s,", stringValue(params.get(param)).c_str());
        break;
      }
      default: {
        fprintf(logfile, "%s,%s\n", lowerCaseString(params.name(param)).c_str(),
            stringValue(params.get(param)).c_str());
      }
    }
  }
  for (int i = 0; i < doubleParam::NUM_DOUBLE_PARAMS; i++) {
    doubleParam param = static_cast<doubleParam>(i);
    switch (amountToPrint) {
      case 1: {
        fprintf(logfile, "%s,", lowerCaseString(params.name(param)).c_str());
        break;
      }
      case 2: {
        fprintf(logfile, "%s,", stringValue(params.get(param)).c_str());
        break;
      }
      default: {
        fprintf(logfile, "%s,%s\n", lowerCaseString(params.name(param)).c_str(),
            stringValue(params.get(param)).c_str());
      }
    }
  }
  for (int i = 0; i < stringParam::NUM_STRING_PARAMS; i++) {
    stringParam param = static_cast<stringParam>(i);
    switch (amountToPrint) {
      case 1: {
        fprintf(logfile, "%s,", lowerCaseString(params.name(param)).c_str());
        break;
      }
      case 2: {
        fprintf(logfile, "%s,", params.get(param).c_str());
        break;
      }
      default: {
        fprintf(logfile, "%s,%s\n", lowerCaseString(params.name(param)).c_str(),
            params.get(param).c_str());
      }
    }
  }
  for (int i = 0; i < static_cast<int>(intConst::NUM_INT_CONST); i++) {
    intConst param = static_cast<intConst>(i);
    switch (amountToPrint) {
      case 1: {
        fprintf(logfile, "%s,", lowerCaseString(params.name(param)).c_str());
        break;
      }
      case 2: {
        fprintf(logfile, "%s,", stringValue(params.get(param)).c_str());
        break;
      }
      default: {
        fprintf(logfile, "%s,%s\n", lowerCaseString(params.name(param)).c_str(),
            stringValue(params.get(param)).c_str());
      }
    }
  }
  for (int i = 0; i < static_cast<int>(doubleConst::NUM_DOUBLE_CONST); i++) {
    doubleConst param = static_cast<doubleConst>(i);
    switch (amountToPrint) {
      case 1: {
        fprintf(logfile, "%s,", lowerCaseString(params.name(param)).c_str());
        break;
      }
      case 2: {
        fprintf(logfile, "%s,", stringValue(params.get(param)).c_str());
        break;
      }
      default: {
        fprintf(logfile, "%s,%s\n", lowerCaseString(params.name(param)).c_str(),
            stringValue(params.get(param)).c_str());
      }
    }
  }
  fflush(logfile);
} /* printParams */
