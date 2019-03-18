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
#include <type_traits> // is_arithmetic
#include <unordered_map>
//#include <functional> // hash

#include "utility.hpp" // parseInt/Double, stringValue, lowerCaseString, and overloading << for vectors

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
  TEMP, // useful for various temporary parameter changes; see corresponding enum
  // Objective options
  USE_ALL_ONES, // 0: do not use, 1: use
  USE_DISJ_LB, // 0: do not use, 1: use
  USE_ITER_BILINEAR, // 0: do not use, 1+: number of iterations to do
  USE_TIGHT_POINTS, // 0: do not use, 1+: num points to try
  USE_TIGHT_RAYS, // 0: do not use, 1+: num rays to try, -1: try the first sqrt(# rays) rays
  USE_UNIT_VECTORS, // 0: do not use, 1+: num to try, <0: abs(val) * sqrt(n)
  // Other options
  VERBOSITY,
  // BB options
  RANDOM_SEED,
  BB_RUNS, // number of times to run b&b (if negative, also test with VPCs+GMICs if GMICs have been generated)
  //  off = 0,
  //  cbc = 2,
  //  cplex = 4,
  //  gurobi = 8,
  //  user_cuts = 16,
  //  all_cuts_off = 32,
  //  all_cuts_on = 64,
  //  gmics_off = 128,
  //  gmics_on = 256,
  //  presolve_off = 512,
  //  presolve_on = 1024,
  //  heuristics_off = 2048,
  //  heuristics_on = 4096,
  //  use_best_bound = 8192,
  //  strong_branching_on = 16384
  BB_STRATEGY, // bit vector; sum of above bits
  BB_MODE, // 111: each bit represents whether to branch with gmics, vpcs, and no cuts (from largest to smallest bit)
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

enum class TempOptions {
  NONE = 0,
  PREPROCESS = 1,
  PREPROCESS_CUSTOM = 2,
  CHECK_CUTS_AGAINST_BB_OPT = 3,
  CALC_NUM_GOMORY_ROUNDS_TO_MATCH = 4,
  // Options for generating tikz string
  GEN_TIKZ_STRING_WITH_VPCS = 10,
  GEN_TIKZ_STRING_WITH_GMICS = 11,
  GEN_TIKZ_STRING_WITH_VPCS_AND_GMICS = 12,
  GEN_TIKZ_STRING_NO_CUTS = 13,
  GEN_TIKZ_STRING_AND_RETURN = 14,
  GEN_TIKZ_STRING_AND_EXIT = 15,
};

/********** DEFINITIONS **********/
//template <class T> class Parameter;
//namespace std {
//  template<class T>
//  struct hash<Parameter<T>> {
//  public:
//    size_t operator()(const Parameter<T> &param) const { return std::hash<std::string>{}(param.name()); };
//  };
//}
/**
 * EnumClassHash is used in order to hash for unordered_map
 */
struct EnumClassHash {
  template<typename T>
  std::size_t operator()(const T& t) const {
    return static_cast<std::size_t>(t);
  }
};

template <class T>
class Parameter {
public:
  Parameter(std::string name, const T& val) : param_name(name), val(val) {}
  virtual ~Parameter() {}
  virtual std::string to_string(const int amountToPrint = 0, const char* fmt = NULL) const = 0;

  virtual const T get() const final { return this->val; }
  virtual std::string name() const final { return this->param_name; }
  virtual bool set(const T& val) final { this->val = val; return check(); }
  virtual bool check() const { return true; }

  bool operator<(const Parameter& other) const {
    return this->param_name.compare(other.param_name) < 0;
  }
  bool operator==(const Parameter& other) const {
    return this->param_name == other.param_name;
  }
protected:
  std::string param_name;
  T val;
}; /* Parameter */

template <class T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
//template <class T>
class NumericParameter : public Parameter<T> {
public:
    using Parameter<T>::Parameter;
    NumericParameter(const std::string& name, const T& val, const T& min_val, const T& max_val)
    : Parameter<T>(name, val), min_val(min_val), max_val(max_val) { check(); }
    NumericParameter(const std::string& name, const T& val, const std::vector<T>& allowed_vals)
    : Parameter<T>(name, val), allowed_vals(allowed_vals) {
      this->min_val = *std::min_element(allowed_vals.begin(), allowed_vals.end());
      this->max_val = *std::max_element(allowed_vals.begin(), allowed_vals.end());
      check();
    }
    virtual ~NumericParameter() {}

    /**
     * to_string
     * amountToPrint:
     *  0 = both names and values,
     *  1 = only names,
     *  2 = only values
     */
    virtual std::string to_string(const int amountToPrint = 2, const char* fmt = NULL) const {
      std::string retval = "";
      switch (amountToPrint) {
        case 1: {
          retval = lowerCaseString(this->name());
          break;
        }
        case 2: {
          retval = stringValue(this->val, fmt);
          break;
        }
        default: {
          retval = lowerCaseString(this->name()) + "," + stringValue(this->val, fmt);
        }
      }
      return retval;
    } /* to_string */

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
          << ": Val = " << stringValue(this->val) << ".";
      if (this->allowed_vals.size() > 0 && this->allowed_vals.size() < 10)
        std::cerr << " Allowed values: " << allowed_vals << ".";
      else
        std::cerr << " Min = " << stringValue(this->min_val)
            << ". Max = " << stringValue(this->max_val) << ".";
      std::cerr << std::endl;
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
  IntParameter(intParam num, const std::string& name, const int& val, const int& min_val, const int& max_val)
  : NumericParameter<int>(name, val, min_val, max_val), param_num(num) {}
  IntParameter(intParam num, const std::string& name, const int& val, const std::vector<int>& allowed_vals)
  : NumericParameter<int>(name, val, allowed_vals), param_num(num) {}
  virtual std::string to_string(const int amountToPrint = 2, const char* fmt = NULL) const {
    if (fmt) {
      return NumericParameter<int>::to_string(amountToPrint, fmt);
    } else {
      return NumericParameter<int>::to_string(amountToPrint, "%d");
    }
  }
protected:
  intParam param_num;
}; /* IntParameter */

class DoubleParameter : public NumericParameter<double> {
public:
  using NumericParameter<double>::NumericParameter;
  DoubleParameter(doubleParam num, const std::string& name, const double& val, const double& min_val, const double& max_val)
  : NumericParameter<double>(name, val, min_val, max_val), param_num(num) {}
  DoubleParameter(doubleParam num, const std::string& name, const double& val, const std::vector<double>& allowed_vals)
  : NumericParameter<double>(name, val, allowed_vals), param_num(num) {}
  virtual std::string to_string(const int amountToPrint = 2, const char* fmt = NULL) const {
    if (fmt) {
      return NumericParameter<double>::to_string(amountToPrint, fmt);
    } else {
      return NumericParameter<double>::to_string(amountToPrint, "%.3e");
    }
  }
protected:
  doubleParam param_num;
}; /* DoubleParameter */

class StringParameter : public Parameter<std::string> {
public:
  using Parameter<std::string>::Parameter;
  StringParameter(stringParam num, const std::string& name, const std::string& val)
  : Parameter<std::string>(name, val), param_num(num) {}

  /**
   * to_string
   * amountToPrint:
   * 0 = all (also adds a newline after each param/constant),
   * 1 = only names,
   * 2 = only values
   */
  virtual std::string to_string(const int amountToPrint = 2, const char* fmt = NULL) const {
    std::string retval = "";
    switch (amountToPrint) {
      case 1: {
        retval = lowerCaseString(this->name());
        break;
      }
      case 2: {
        retval = stringValue(this->val, fmt);
        break;
      }
      default: {
        retval = lowerCaseString(this->name()) + "," + stringValue(this->val, fmt);
      }
    }
    return retval;
  } /* to_string */
protected:
  stringParam param_num;
}; /* StringParameter */

/********** VPC PARAMETERS STRUCT **********/
struct VPCParameters {
  FILE* logfile = NULL; // NB: right now this is a shallow copy if this struct gets copied

  /** unordered_map gets printed in reverse order; advantage over map is constant access time on average */
  std::unordered_map<intParam, IntParameter, EnumClassHash> intParamValues {
    {intParam::BB_MODE, IntParameter("BB_MODE", 10, 0, 111)}, // 010 = branch with vpcs only
    {intParam::BB_STRATEGY, IntParameter("BB_STRATEGY", 10776, std::numeric_limits<int>::min(), std::numeric_limits<int>::max())}, // see BBHelper.hpp; 10776 = 010101000011000 => gurobi: 1, user_cuts: 1, presolve_off: 1, heuristics_off: 1, use_best_bound: 1
    {intParam::BB_RUNS, IntParameter("BB_RUNS", 0, std::numeric_limits<int>::min(), std::numeric_limits<int>::max())}, // see BBHelper.hpp; 10776 = 010101000011000 => gurobi: 1, user_cuts: 1, presolve_off: 1, heuristics_off: 1, use_best_bound: 1
    {intParam::RANDOM_SEED, IntParameter("RANDOM_SEED", 628, -1, std::numeric_limits<int>::max())},
#ifdef TRACE
    {intParam::VERBOSITY, IntParameter("VERBOSITY", 1, 0, 2)},
#else
    {intParam::VERBOSITY, IntParameter("VERBOSITY", 0, 0, 2)},
#endif
    {intParam::USE_UNIT_VECTORS, IntParameter("USE_UNIT_VECTORS", 0, std::numeric_limits<int>::min(), std::numeric_limits<int>::max())},
    {intParam::USE_TIGHT_RAYS, IntParameter("USE_TIGHT_RAYS", 0, std::numeric_limits<int>::min(), std::numeric_limits<int>::max())},
    {intParam::USE_TIGHT_POINTS, IntParameter("USE_TIGHT_POINTS", 0, std::numeric_limits<int>::min(), std::numeric_limits<int>::max())},
    {intParam::USE_ITER_BILINEAR, IntParameter("USE_ITER_BILINEAR", 1, 0, std::numeric_limits<int>::max())},
    {intParam::USE_DISJ_LB, IntParameter("USE_DISJ_LB", 1, 0, 1)},
    {intParam::USE_ALL_ONES, IntParameter("USE_ALL_ONES", 1, 0, 1)},
    {intParam::TEMP, IntParameter("TEMP", 0, std::numeric_limits<int>::min(), std::numeric_limits<int>::max())},
    {intParam::STRENGTHEN, IntParameter("STRENGTHEN", 1, 0, 2)},
    {intParam::ROUNDS, IntParameter("ROUNDS", 1, 0, std::numeric_limits<int>::max())},
    {intParam::PRLP_FLIP_BETA, IntParameter("PRLP_FLIP_BETA", 0, -1, 1)},
    {intParam::PARTIAL_BB_NUM_STRONG, IntParameter("PARTIAL_BB_NUM_STRONG", 5, std::numeric_limits<int>::min(), std::numeric_limits<int>::max())},
    {intParam::PARTIAL_BB_STRATEGY, IntParameter("PARTIAL_BB_STRATEGY", 4, std::numeric_limits<int>::min(), std::numeric_limits<int>::max())},
    {intParam::MODE, IntParameter("MODE", 0, {0, 1, 3})},
    {intParam::DISJ_TERMS, IntParameter("DISJ_TERMS", 0, std::numeric_limits<int>::min(), std::numeric_limits<int>::max())},
    {intParam::CUTLIMIT, IntParameter("CUTLIMIT", -1, std::numeric_limits<int>::min(), std::numeric_limits<int>::max())},
  }; /* intParamValues */
  std::unordered_map<doubleParam, DoubleParameter, EnumClassHash> doubleParamValues {
    {doubleParam::TIMELIMIT, DoubleParameter("TIMELIMIT", 60, 0., std::numeric_limits<double>::max())},
    {doubleParam::PRLP_TIMELIMIT, DoubleParameter("PRLP_TIMELIMIT", -1, -1., std::numeric_limits<double>::max())},
    {doubleParam::PARTIAL_BB_TIMELIMIT, DoubleParameter("PARTIAL_BB_TIMELIMIT", 3600, 0., std::numeric_limits<double>::max())},
    {doubleParam::MIN_ORTHOGONALITY, DoubleParameter("MIN_ORTHOGONALITY", 0., 0., 1.)},
    {doubleParam::EPS, DoubleParameter("EPS", 1e-7, 0., 1.)},
  }; /* doubleParamValues */
  std::unordered_map<stringParam, StringParameter, EnumClassHash> stringParamValues {
    {stringParam::OPTFILE, StringParameter("OPTFILE", "")},
    {stringParam::LOGFILE, StringParameter("LOGFILE", "")},
    {stringParam::FILENAME, StringParameter("FILENAME", "")},
  }; /* stringParamValues */

  // Constants
  std::unordered_map<intConst, IntParameter, EnumClassHash> intConstValues {
    {intConst::PRLP_PRESOLVE, IntParameter("PRLP_PRESOLVE", 2, 2, 2)}, // use presolve when solving PRLP (either initial or resolve)
    {intConst::NB_SPACE, IntParameter("NB_SPACE", 1, 1, 1)}, // currently only works with true
    {intConst::NUM_OBJ_PER_POINT, IntParameter("NUM_OBJ_PER_POINT", -2, -2, -2)}, // sqrt(n)
    {intConst::MODE_OBJ_PER_POINT, IntParameter("MODE_OBJ_PER_POINT", 121, 121, 121)}, // one ray at a time, points+rays+variables, large to small angle (descending)
    {intConst::MAX_SUPPORT_ABS, IntParameter("MAX_SUPPORT_ABS", std::numeric_limits<int>::max(), std::numeric_limits<int>::max(), std::numeric_limits<int>::max())},
    {intConst::LUB, IntParameter("LUB", 1e3, 1e3, 1e3)},
    {intConst::CHECK_DUPLICATES, IntParameter("CHECK_DUPLICATES", 1, 1, 1)},
  }; /* intConstValues */
  std::unordered_map<doubleConst, DoubleParameter, EnumClassHash> doubleConstValues {
    {doubleConst::MAX_SUPPORT_REL, DoubleParameter("MAX_SUPPORT_REL", 0.9, 0.9, 0.9)},
    {doubleConst::MAX_DYN_LUB, DoubleParameter("MAX_DYN_LUB", 1e13, 1e13, 1e13)},
    {doubleConst::MAX_DYN, DoubleParameter("MAX_DYN", 1e8, 1e8, 1e8)},
    {doubleConst::MIN_VIOL_REL, DoubleParameter("MIN_VIOL_REL", 1e-7, 1e-7, 1e-7)},
    {doubleConst::MIN_VIOL_ABS, DoubleParameter("MIN_VIOL_ABS", 1e-7, 1e-7, 1e-7)},
    {doubleConst::EPS_COEFF_LUB, DoubleParameter("EPS_COEFF_LUB", 1e-13, 1e-13, 1e-13)},
    {doubleConst::EPS_COEFF, DoubleParameter("EPS_COEFF", 1e-5, 1e-5, 1e-5)},
    {doubleConst::MIN_PRLP_TIMELIMIT, DoubleParameter("MIN_PRLP_TIMELIMIT", 5., 5., 5.)},
    {doubleConst::BB_TIMELIMIT, DoubleParameter("BB_TIMELIMIT", 3600., 3600., 3600.)},
    {doubleConst::RAYEPS, DoubleParameter("RAYEPS", 1e-7, 1e-7, 1e-7)},
    {doubleConst::INF, DoubleParameter("INF", std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max())},
    {doubleConst::DIFFEPS, DoubleParameter("DIFFEPS", 1e-3, 1e-3, 1e-3)}, // to check whether something is different enough to throw an error
    {doubleConst::AWAY, DoubleParameter("AWAY", 1e-3, 1e-3, 1e-3)},
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

  void        set(intParam param, const int value) { intParamValues.find(param)->second.set(value); }
  void        set(doubleParam param, const double value) { doubleParamValues.find(param)->second.set(value); }
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
 * amountToPrint:
 *  0 = all param name/values except string params (newline-separated)
 *  1 = only param names (comma-separated)
 *  2 = only param values (comma-separated)
 *  3 = string param name/values (newline-sep)
 *  4 = only string names (comma-sep)
 *  5 = only string values (comma-sep)
 *  10 = all const name/values (newline-separated)
 *  11 = only const names (comma-separated)
 *  12 = only const values (comma-separated)
 */
inline void printParams(const VPCParameters& params, FILE* logfile = stdout, const int amountToPrint = 0) {
  if (!logfile)
    return;

//  !(std::find(allowed_vals.begin(), allowed_vals.end(), test_val) == allowed_vals.end())
  const bool printIntParams = amountToPrint == 0 || amountToPrint == 1 || amountToPrint == 2;
  const bool printDoubleParams = amountToPrint == 0 || amountToPrint == 1 || amountToPrint == 2;
  const bool printStringParams = amountToPrint == 3 || amountToPrint == 4 || amountToPrint == 5;
  const bool printIntConsts = amountToPrint == 10 || amountToPrint == 11 || amountToPrint == 12;
  const bool printDoubleConsts = amountToPrint == 10 || amountToPrint == 11 || amountToPrint == 12;

  if (printIntParams)
  for (auto param : params.intParamValues) {
    if (amountToPrint == 0) {
      fprintf(logfile, "%s\n", param.second.to_string(amountToPrint).c_str());
    } else {
      fprintf(logfile, "%s,", param.second.to_string(amountToPrint).c_str());
    }
  }
  if (printDoubleParams)
  for (auto param : params.doubleParamValues) {
    if (amountToPrint == 0) {
      fprintf(logfile, "%s\n", param.second.to_string(amountToPrint).c_str());
    } else {
      fprintf(logfile, "%s,", param.second.to_string(amountToPrint).c_str());
    }
  }
  if (printStringParams)
  for (auto param : params.stringParamValues) {
    if (amountToPrint == 0) {
      fprintf(logfile, "%s\n", param.second.to_string(amountToPrint).c_str());
    } else {
      fprintf(logfile, "%s,", param.second.to_string(amountToPrint).c_str());
    }
  }
  if (printIntConsts)
  for (auto param : params.intConstValues) {
    if (amountToPrint == 0) {
      fprintf(logfile, "%s\n", param.second.to_string(amountToPrint).c_str());
    } else {
      fprintf(logfile, "%s,", param.second.to_string(amountToPrint).c_str());
    }
  }
  if (printDoubleConsts)
  for (auto param : params.doubleConstValues) {
    if (amountToPrint == 0) {
      fprintf(logfile, "%s\n", param.second.to_string(amountToPrint).c_str());
    } else {
      fprintf(logfile, "%s,", param.second.to_string(amountToPrint).c_str());
    }
  }
  fflush(logfile);
} /* printParams */
