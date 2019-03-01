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
 * 1. Add it to the relevant enum (intParam, doubleParam, stringParam, intConst, doubleConst)
 * 2. In the struct VPCParameters, add the parameter using the name/default/min/max constructor of the parameter
 * 3. Optionally, add a way to set the parameter in the option handling part of the code (in the main file)
 ********************************************************************************************************************/

#include <map>
#include <string>
#include <vector>
#include <cstdio> // fprintf
#include <fstream>
#include <sstream>
#include <algorithm> // min_element, max_element
//#include <type_traits>

#include "utility.hpp" // parseInt/Double, stringValue, and lowerCaseString

// Define SolverInterface that we can later change in a way that will be used in all VPC files
#ifdef USE_CLP
  #include <OsiClpSolverInterface.hpp>
  using SolverInterface = OsiClpSolverInterface;
#else
  #include <OsiSolverInterface.hpp>
  using SolverInterface = OsiSolverInterface;
#endif

// Definitions
template <class T>
class Parameter {
public:
  Parameter(std::string name, const T& val) : param_name(name), val(val) {}
  virtual ~Parameter() {}
  virtual std::string to_string(const char* fmt = NULL) const = 0;

  virtual const T get() const { return this->val; }
  virtual std::string name() const { return this->param_name; }
  virtual bool set(const T& val) { this->val = val; return check(); }
  virtual bool check() const { return true; }
protected:
  std::string param_name;
  T val;
}; /* Parameter */

//template <class T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
template <class T>
class NumericParameter : public Parameter<T> {
public:
    using Parameter<T>::Parameter;
    NumericParameter(const std::string& name, const T& val, const T& min_val, const T& max_val)
    : Parameter<T>(name, val), min_val(min_val), max_val(max_val) { check(); }
    NumericParameter(const std::string& name, const T& val, const std::vector<T>& allowed_vals)
    : Parameter<T>(name, val), allowed_vals(max_val) {
      this->min_val = *std::min_element(allowed_vals.begin(), allowed_vals.end());
      this->max_val = *std::max_element(allowed_vals.begin(), allowed_vals.end());
      check();
    }
    virtual ~NumericParameter() {}

    virtual std::string to_string(const char* fmt = NULL) const {
      if (fmt) {
        char tmp[25];
        snprintf(tmp, sizeof(tmp) / sizeof(char), fmt, this->val);
        std::string str(tmp);
        return str;
      } else {
        return std::to_string(this->val);
      }
    }

    virtual const T get_min() const { return this->min_val; }
    virtual const T get_max() const { return this->max_val; }
    virtual bool val_is_allowed(const T test_val) const {
      if (!this->allowed_vals.empty()) {
        return !(std::find(allowed_vals.begin(), allowed_vals.end(), test_val) == allowed_vals.end());
      } else {
        return !(lessThanVal(this->val, this->min_val) || greaterThanVal(this->val, this->max_val));
      }
    }
    virtual const std::vector<T>& get_allowed_vals() const { return this->allowed_vals; }
    virtual bool check() const {
      if (val_is_allowed(this->val)) {
        return true;
      } else {
      std::cerr << "*** ERROR: Error setting parameter " << this->name()
          << ": Val = " << stringValue(this->val)
          << ". Min = " << stringValue(this->min_val)
          << ". Max = " << stringValue(this->max_val)
          << "." << std::endl;
        exit(1);
      }
    }

protected:
    T min_val, max_val;
    std::vector<T> allowed_vals;
}; /* NumericParameter */

class IntParameter : public NumericParameter<int> {
public:
  using NumericParameter<int>::NumericParameter;
  virtual std::string to_string(const char* fmt = NULL) const {
    if (fmt) {
      return NumericParameter<int>::to_string(fmt);
    } else {
      return NumericParameter<int>::to_string("%d");
    }
  }
}; /* IntParameter */

class DoubleParameter : public NumericParameter<double> {
public:
  using NumericParameter<double>::NumericParameter;
  virtual std::string to_string(const char* fmt = NULL) const {
    if (fmt) {
      return NumericParameter<double>::to_string(fmt);
    } else {
      return NumericParameter<double>::to_string("%.3e");
    }
  }
}; /* DoubleParameter */

class StringParameter : public Parameter<std::string> {
public:
  using Parameter<std::string>::Parameter;
  virtual ~StringParameter() {}
  virtual std::string to_string(const char* fmt = NULL) const {
    if (fmt) {
      char tmp[25];
      snprintf(tmp, sizeof(tmp) / sizeof(char), fmt, val.c_str());
      std::string str(tmp);
      return str;
    } else {
      return val;
    }
  }
}; /* StringParameter */

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
enum doubleParam {
  EPS,
  MIN_ORTHOGONALITY, // minimum orthogonality between cuts added to the collection
  PARTIAL_BB_TIMELIMIT,
  PRLP_TIMELIMIT, // -1 = unlimited time for initial solve, and then half that time subsequently; -2 = unlimited time always
  TIMELIMIT,
  NUM_DOUBLE_PARAMS
}; /* doubleParam */
enum stringParam {
  FILENAME,
  LOGFILE,
  OPTFILE,
  NUM_STRING_PARAMS
}; /* stringParam */

/********** CONSTANTS **********/
enum class intConst {
  CHECK_DUPLICATES, // do not add duplicate cuts
  LUB, // value for var upper bound considered "large"
  MAX_SUPPORT_ABS,
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
  MODE_OBJ_PER_POINT,
  NUM_OBJ_PER_POINT, // # cuts to try to generate for the strong branching lb point (and others)
  NB_SPACE, // whether to generate cuts in the nonbasic space (currently must be set to true)
  PRLP_PRESOLVE, // 0: no presolve, 1: only initial solve, 2: both initial solve and resolve
  RANDOM_SEED,
  NUM_INT_CONST
}; /* intConst */
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

/********** VPC PARAMETERS STRUCT **********/
struct VPCParameters {
  FILE* logfile = NULL; // NB: right now this is a shallow copy if this struct gets copied

  std::map<intParam, IntParameter> intParamValues {
    {intParam::CUTLIMIT, IntParameter("CUTLIMIT", -1, std::numeric_limits<int>::min(), std::numeric_limits<int>::max())},
    {intParam::DISJ_TERMS, IntParameter("DISJ_TERMS", 0, std::numeric_limits<int>::min(), std::numeric_limits<int>::max())},
    {intParam::MODE, IntParameter("MODE", 0, 0, std::numeric_limits<int>::max())},
    {intParam::PARTIAL_BB_STRATEGY, IntParameter("PARTIAL_BB_STRATEGY", 4, std::numeric_limits<int>::min(), std::numeric_limits<int>::max())},
    {intParam::PARTIAL_BB_NUM_STRONG, IntParameter("PARTIAL_BB_NUM_STRONG", 5, std::numeric_limits<int>::min(), std::numeric_limits<int>::max())},
    {intParam::PRLP_FLIP_BETA, IntParameter("PRLP_FLIP_BETA", 0, -1, 1)},
    {intParam::ROUNDS, IntParameter("ROUNDS", 1, 0, std::numeric_limits<int>::max())},
    {intParam::STRENGTHEN, IntParameter("STRENGTHEN", 1, 0, 2)},
    {intParam::TEMP, IntParameter("TEMP", 0, std::numeric_limits<int>::min(), std::numeric_limits<int>::max())},
    {intParam::USE_ALL_ONES, IntParameter("USE_ALL_ONES", 1, 0, 1)},
    {intParam::USE_DISJ_LB, IntParameter("USE_DISJ_LB", 1, 0, 1)},
    {intParam::USE_ITER_BILINEAR, IntParameter("USE_ITER_BILINEAR", 0, 0, std::numeric_limits<int>::max())},
    {intParam::USE_TIGHT_POINTS, IntParameter("USE_TIGHT_POINTS", 0, std::numeric_limits<int>::min(), std::numeric_limits<int>::max())},
    {intParam::USE_TIGHT_RAYS, IntParameter("USE_TIGHT_RAYS", 0, std::numeric_limits<int>::min(), std::numeric_limits<int>::max())},
    {intParam::USE_UNIT_VECTORS, IntParameter("USE_UNIT_VECTORS", 0, std::numeric_limits<int>::min(), std::numeric_limits<int>::max())},
#ifdef TRACE
    {intParam::VERBOSITY, IntParameter("VERBOSITY", 1, 0, 2)},
#else
    {intParam::VERBOSITY, IntParameter("VERBOSITY", 0, 0, 2)},
#endif
  }; /* intParamValues */
  std::map<doubleParam, DoubleParameter> doubleParamValues {
    {doubleParam::EPS, DoubleParameter("EPS", 1e-7, 0., 1.)},
    {doubleParam::MIN_ORTHOGONALITY, DoubleParameter("MIN_ORTHOGONALITY", 0., 0., 1.)},
    {doubleParam::PARTIAL_BB_TIMELIMIT, DoubleParameter("PARTIAL_BB_TIMELIMIT", 3600, 0., std::numeric_limits<double>::max())},
    {doubleParam::PRLP_TIMELIMIT, DoubleParameter("PRLP_TIMELIMIT", -1, -1., std::numeric_limits<double>::max())},
    {doubleParam::TIMELIMIT, DoubleParameter("TIMELIMIT", 60, 0., std::numeric_limits<double>::max())},
  }; /* doubleParamValues */
  std::map<stringParam, StringParameter> stringParamValues {
    {stringParam::FILENAME, StringParameter("FILENAME", "")},
    {stringParam::LOGFILE, StringParameter("LOGFILE", "")},
    {stringParam::OPTFILE, StringParameter("OPTFILE", "")}
  }; /* stringParamValues */

  // Constants
  std::map<intConst, IntParameter> intConstValues {
    {intConst::CHECK_DUPLICATES, IntParameter("CHECK_DUPLICATES", 1, 1, 1)},
    {intConst::LUB, IntParameter("LUB", 1e3, 1e3, 1e3)},
    {intConst::MAX_SUPPORT_ABS, IntParameter("MAX_SUPPORT_ABS", std::numeric_limits<int>::max(), std::numeric_limits<int>::max(), std::numeric_limits<int>::max())},
    {intConst::MODE_OBJ_PER_POINT, IntParameter("MODE_OBJ_PER_POINT", 121, 121, 121)}, // one ray at a time, points+rays+variables, large to small angle (descending)
    {intConst::NUM_OBJ_PER_POINT, IntParameter("NUM_OBJ_PER_POINT", -2, -2, -2)}, // sqrt(n)
    {intConst::NB_SPACE, IntParameter("NB_SPACE", 1, 1, 1)}, // currently only works with true
    {intConst::PRLP_PRESOLVE, IntParameter("PRLP_PRESOLVE", 2, 2, 2)}, // use presolve when solving PRLP (either initial or resolve)
    {intConst::RANDOM_SEED, IntParameter("RANDOM_SEED", 628, 628, 628)}
  }; /* intConstValues */
  std::map<doubleConst, DoubleParameter> doubleConstValues {
    {doubleConst::AWAY, DoubleParameter("AWAY", 1e-3, 1e-3, 1e-3)},
    {doubleConst::DIFFEPS, DoubleParameter("DIFFEPS", 1e-3, 1e-3, 1e-3)}, // to check whether something is different enough to throw an error
    {doubleConst::INF, DoubleParameter("INF", std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max())},
    {doubleConst::RAYEPS, DoubleParameter("RAYEPS", 1e-7, 1e-7, 1e-7)},
//    {doubleConst::BB_STRATEGY, 0}, // see BBHelper.hpp; 10776 = 010101000011000 => gurobi: 1, user_cuts: 1, presolve_off: 1, heuristics_off: 1, use_best_bound: 1
    {doubleConst::BB_TIMELIMIT, DoubleParameter("BB_TIMELIMIT", 3600., 3600., 3600.)},
    {doubleConst::MIN_PRLP_TIMELIMIT, DoubleParameter("MIN_PRLP_TIMELIMIT", 5., 5., 5.)},
    {doubleConst::EPS_COEFF, DoubleParameter("EPS_COEFF", 1e-5, 1e-5, 1e-5)},
    {doubleConst::EPS_COEFF_LUB, DoubleParameter("EPS_COEFF_LUB", 1e-13, 1e-13, 1e-13)},
    {doubleConst::MIN_VIOL_ABS, DoubleParameter("MIN_VIOL_ABS", 1e-7, 1e-7, 1e-7)},
    {doubleConst::MIN_VIOL_REL, DoubleParameter("MIN_VIOL_REL", 1e-7, 1e-7, 1e-7)},
    {doubleConst::MAX_DYN, DoubleParameter("MAX_DYN", 1e8, 1e8, 1e8)},
    {doubleConst::MAX_DYN_LUB, DoubleParameter("MAX_DYN_LUB", 1e13, 1e13, 1e13)},
    {doubleConst::MAX_SUPPORT_REL, DoubleParameter("MAX_SUPPORT_REL", 0.9, 0.9, 0.9)}
  }; /* doubleConstValues */

  // Methods (name/get/set)
  std::string name(intParam param) const { return intParamValues.find(param)->second.name(); }
  std::string name(doubleParam param) const { return doubleParamValues.find(param)->second.name(); }
  std::string name(stringParam param) const { return stringParamValues.find(param)->second.name(); }
  std::string name(intConst param) const { return intConstValues.find(param)->second.name(); }
  std::string name(doubleConst param) const { return doubleConstValues.find(param)->second.name(); }

  int         get(intParam param) const { return intParamValues.find(param)->second.get(); }
  double      get(doubleParam param) const { return doubleParamValues.find(param)->second.get(); }
  std::string get(stringParam param) const { return stringParamValues.find(param)->second.get(); }
  int         get(intConst param) const { return intConstValues.find(param)->second.get(); }
  double      get(doubleConst param) const { return doubleConstValues.find(param)->second.get(); }

  void        set(intParam param, int value) { intParamValues.find(param)->second.set(value); }
  void        set(doubleParam param, double value) { doubleParamValues.find(param)->second.set(value); }
  void        set(stringParam param, std::string value) { stringParamValues.find(param)->second.set(value); }

  /** Set parameters by name */
  bool        set(std::string tmpname, const std::string tmp) {
    std::string name = upperCaseStringNoUnderscore(tmpname);
//    for (unsigned i = 0; i < intParam::NUM_INT_PARAMS; i++) {
    for (unsigned i = 0; i < intParam::NUM_INT_PARAMS; i++) {
      intParam param = static_cast<intParam>(i);
      std::string name = intParamValues.find(param)->second.name();
      std::string str1 = upperCaseStringNoUnderscore(name);
      if (str1.compare(name) == 0) {
        int value;
        if (tmp.compare("+inf") == 0) {
          value = std::numeric_limits<int>::max();
        } else if (tmp.compare("-inf") == 0) {
          value = std::numeric_limits<int>::min();
        } else {
          value = parseInt(tmp.c_str(), value);
        }
        set(param, value);
        return true;
      }
    }
    for (unsigned i = 0; i < doubleParam::NUM_DOUBLE_PARAMS; i++) {
      doubleParam param = static_cast<doubleParam>(i);
      std::string name = doubleParamValues.find(param)->second.name();
      std::string str1 = upperCaseStringNoUnderscore(name);
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
    for (unsigned i = 0; i < stringParam::NUM_STRING_PARAMS; i++) {
      stringParam param = static_cast<stringParam>(i);
      std::string name = stringParamValues.find(param)->second.name();
      std::string str1 = upperCaseStringNoUnderscore(name);
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
