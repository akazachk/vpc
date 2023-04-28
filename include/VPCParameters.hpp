/**
 * @file VPCParameters.hpp
 * @author A. M. Kazachkov
 * @date 2019-Feb-20
 */
#pragma once

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

/**
 * @brief Namespace to contain all VPC parameters (and to enable interaction with other codes using generic names like #intParam)
 * @details This file contains the parameters/constants that are used in the CglVPC class
 * (also in the PRLP and Disjunction classes)
 *
 * To add a new parameter/constant:
 * 1. Add it to the relevant enum (#intParam, #doubleParam, #stringParam, #intConst, #doubleConst)
 * 2. In the struct VPCParameters, add the parameter using the name/default/min/max constructor of the parameter
 * 3. Optionally, add a way to set the parameter in the option handling part of the code (in the main file)
 **/
namespace VPCParametersNamespace {
#define ENUM_OPTION_0 1
#define ENUM_OPTION_1 2
#define ENUM_OPTION_2 4
#define ENUM_OPTION_3 8
#define ENUM_OPTION_4 16
#define ENUM_OPTION_5 32
#define ENUM_OPTION_6 64
#define ENUM_OPTION_7 128
#define ENUM_OPTION_8 256
#define ENUM_OPTION_9 512
#define ENUM_OPTION_10 1024
#define ENUM_OPTION_11 2048
#define ENUM_OPTION_12 4096
#define ENUM_OPTION_13 8192
#define ENUM_OPTION_14 16384
#define ENUM_OPTION_15 32768

/********** PARAMETERS **********/
/// Integer-valued parameters
enum intParam {
  CUTLIMIT, ///< max number of cuts generated; 0 = none, -k = k * # fractional variables at root
  DISJ_TERMS, ///< number of disjunctive terms or number of disjunctions, depending on ::MODE
  GOMORY, ///< Gomory cut mode, 0: none, +/-1: use CglGMI class to generate cuts (-1: do not add them to LP before generating VPCs; 1: do add them)
  MODE, ///< 0: partial b&b tree, 1: splits, 2: crosses (not implemented), 3: custom
  /// Total is used to decide the choose variable/branch/node decision:
  /// \li variable decision => hundreds digit: 0: default, 1: default+second criterion, 2: max min change+second (max max change), 3: second-best default, 4: second-best max-min change, -x: -1 * (1+x);
  /// \li branch decision => tens digit: 0: default, 1: dynamic, 2: strong, 3: none;
  /// \li node comparison decision => ones digit: 0: default: 1: bfs, 2: depth, 3: estimate, 4: objective
  PARTIAL_BB_STRATEGY,
  PARTIAL_BB_NUM_STRONG, ///< -1: num cols, -2: sqrt(num cols), >= 0: that many
  PREPROCESS, ///< 0: off, 1: on (with solver), 2: on (solver + custom cleaning process)
  PRLP_FLIP_BETA, ///< controls rhs in nb space, -1: do not cut away LP opt, 0: cut away LP opt, 1: both
  ROUNDS, ///< number of VPC rounds to do
  STRENGTHEN, ///< 0: no, 1: yes, when possible, 2: same as 1 plus add GMICs to strengthen each disjunctive term
  TEMP, ///< useful for various temporary parameter changes; see corresponding enum
  // Objective options
  USE_ALL_ONES, ///< 0: do not use, 1: use
  USE_DISJ_LB, ///< 0: do not use, 1: use
  USE_ITER_BILINEAR, ///< 0: do not use, 1+: number of iterations to do
  USE_TIGHT_POINTS, ///< 0: do not use, 1+: num points to try
  USE_TIGHT_RAYS, ///< 0: do not use, 1+: num rays to try, -1: try the first sqrt(# rays) rays
  USE_UNIT_VECTORS, ///< 0: do not use, 1+: num to try, <0: abs(val) * sqrt(n)
  // Other options
  VERBOSITY, ///< how much to print
  // BB options
  RANDOM_SEED, ///< random seed to use in b&b experiments (and elsewhere); this will be multiplied by i = 1 to BB_RUNS if multiple runs are used
  BB_RUNS, ///< number of times to run b&b (if negative, also test with VPCs+GMICs if GMICs have been generated)
  /// bit vector; sum of bits, default: 536 = gurobi (8) + user_cuts (16) + presolve_off (512).
  /// \li  off = 0,
  /// \li  cbc = 2,
  /// \li  cplex = 4,
  /// \li  gurobi = 8,
  /// \li  user_cuts = 16,
  /// \li  all_cuts_off = 32,
  /// \li  all_cuts_on = 64,
  /// \li  gmics_off = 128,
  /// \li  gmics_on = 256,
  /// \li  presolve_off = 512,
  /// \li  presolve_on = 1024,
  /// \li  heuristics_off = 2048,
  /// \li  heuristics_on = 4096,
  /// \li  use_best_bound = 8192,
  /// \li  strong_branching_on = 16384
  BB_STRATEGY, 
  BB_MODE, ///< 111: each bit represents whether to branch with gmics, vpcs, and no cuts (from largest to smallest bit)
  NUM_INT_PARAMS ///< number of integer params
}; /* intParam */
/// Double-valued parameters
enum doubleParam {
  BB_TIMELIMIT, ///< time limit for doing branch-and-bound
  EPS, ///< global epsilon (may be further refined based on instance-specific data)
  INF, ///< infinity (INFINITY keyword is reserved as a macro from math header)
  IP_OBJ, ///< way to give just the objective for this instance rather than reading it from a file
  MIN_ORTHOGONALITY, ///< minimum orthogonality between cuts added to the collection
  PARTIAL_BB_TIMELIMIT, ///< time allotted for generating partial b&b tree
  PRLP_TIMELIMIT, ///< -1 = unlimited time for initial solve, and then half that time subsequently; -2 = unlimited time always
  TIMELIMIT, ///< time limit used for overall cut generation process
  NUM_DOUBLE_PARAMS ///< number of double params
}; /* doubleParam */
/// String-valued parameters
enum stringParam {
  DISJ_OPTIONS,
  FILENAME,
  LOGFILE,
  OPTFILE,
  SOLFILE,
  NUM_STRING_PARAMS
}; /* stringParam */

/********** CONSTANTS **********/
/// Integer parameters that we do not let the user change
enum class intConst {
  CHECK_DUPLICATES, ///< do not add duplicate cuts
  LUB, ///< value for var upper bound considered "large"
  MAX_SUPPORT_ABS, ///< max absolute number of nonzero coefficients in any cut we generate
  /// VPC objective-related constants
  /// For each point, MODE_OBJ_PER_POINT says which objectives to try
  /// \li 0 = all of the rows that are not tight, and subtract as they get tight
  /// \li 1 = one point/ray at time
  /// \li 2 = keep trying even if the point/ray has become tight in the process?
  /// \li 0 is more expensive at each step and not clear that it would be better
  /// \li 0x: only rays
  /// \li 1x: points+rays
  /// \li 2x: points+rays+variables
  /// \li 0xx: small to large angle with obj (ascending)
  /// \li 1xx: large to small angle with obj (descending)
  /// \li 2xx: small to large slack (ascending)
  /// \li 3xx: large to small slack (descending)
  MODE_OBJ_PER_POINT,
  NUM_OBJ_PER_POINT, ///< number cuts to try to generate for the strong branching lb point (and others)
  NB_SPACE, ///< whether to generate cuts in the nonbasic space (currently must be set to true)
  PRLP_PRESOLVE, ///< 0: no presolve, 1: only initial solve, 2: both initial solve and resolve
  NUM_INT_CONST ///< number of integer constants
}; /* intConst */
/// Double parameters that we do not let the user change
enum class doubleConst {
  AWAY, ///< min fractionality before a variable is considered not integer-valued
  DIFFEPS, ///< to check whether something is different enough to throw an error
  RAYEPS, ///< value for which a ray coefficient will be treated as zero
  // Time limits
  MIN_PRLP_TIMELIMIT, ///< minimum amount of time allotted for solving/resolving PRLP
  // Safety related constants:
  EPS_COEFF, ///< any cut coefficient smaller than this will be replaced by zero
  EPS_COEFF_LUB, ///< for variables with large upper bound, any cut coefficient smaller than this will be replaced by zero
  MIN_VIOL_ABS, ///< minimum amount a cut should remove the LP opt
  MIN_VIOL_REL, ///< minimum amount a cut should remove the LP opt (using ratio)
  MAX_DYN, ///< |alpha_max| / |alpha_min| upper bound; Maximum ratio between largest and smallest non zero coefficients in a cut
  MAX_DYN_LUB, ///< Same as doubleConst::MAX_DYN but when some of the variables involved in the cut have a large upper bound; should be >= doubleConst::MAX_DYN logically
  MAX_SUPPORT_REL, ///< max relative (to number of variables) number of nonzero coefficients in any cut we generate
  NUM_DOUBLE_CONST ///< number of double constants
}; /* doubleConst */

/// @brief Temporary options that might be invoked during testing or running of the code
enum class TempOptions {
  NONE = 0, ///< default
  CHECK_CUTS_AGAINST_BB_OPT = ENUM_OPTION_3, ///< if integer optimal solution is available, check if it violates any cuts
  CALC_NUM_GOMORY_ROUNDS_TO_MATCH = ENUM_OPTION_4, ///< calculate number of rounds of GMICs needed to get same bound as obtained by our cuts
  SAVE_IP_OPT = ENUM_OPTION_5, ///< save IP optimum solution
  // Options for generating tikz string (if negative, then print nodes without the LP value at the node)
  GEN_TIKZ_STRING = ENUM_OPTION_10, ///< generate code for printing (partial) b&b tree in tikz format
  GEN_TIKZ_STRING_WITH_VPCS = ENUM_OPTION_11, ///< generate code for printing b&b tree in tikz format, after our cuts have been added
  GEN_TIKZ_STRING_WITH_GMICS = ENUM_OPTION_12, ///< generate code for printing b&b tree in tikz format, after GMICs have been added
//  GEN_TIKZ_STRING_WITH_VPCS_AND_GMICS = (ENUM_OPTION_10 & ENUM_OPTION_11),
  GEN_TIKZ_STRING_AND_RETURN = ENUM_OPTION_14, ///< after generating tikz string, stop cut generation process, but do not quit
  GEN_TIKZ_STRING_AND_EXIT = ENUM_OPTION_15, ///< quit after generating tikz string
};

/// @brief Shortcut for checking if a bit is enabled
inline bool use_temp_option(const int strategy, const TempOptions option) {
  return strategy & static_cast<int>(option);
}
/// @brief Shortcut for bitwise enabling an option
inline int enable_temp_option(const int strategy, const TempOptions option) {
  return strategy | static_cast<int>(option);
}
/// @brief Shortcut for bitwise disabling an option
inline int disable_temp_option(const int strategy, const TempOptions option) {
  return strategy & ~static_cast<int>(option);
}

/// @brief Options for the parameter BB_STRATEGY
enum class BB_Strategy_Options {
  off                 = 0,              ///< do not do any b&b tests
  cbc                 = ENUM_OPTION_1,  ///< use Cbc as the branch-and-bound solver
  cplex               = ENUM_OPTION_2,  ///< use CPLEX as the branch-and-bound solver
  gurobi              = ENUM_OPTION_3,  ///< use Gurobi as the branch-and-bound solver
  user_cuts           = ENUM_OPTION_4,  ///< tell solver that user cuts may be used
  all_cuts_off        = ENUM_OPTION_5,  ///< in b&b, do not use any cuts
  all_cuts_on         = ENUM_OPTION_6,  ///< in b&b, use all possible (default) cuts
  gmics_off           = ENUM_OPTION_7,  ///< in b&b, turn off gomory cuts
  gmics_on            = ENUM_OPTION_8,  ///< in b&b, turn on gomory cuts
  presolve_off        = ENUM_OPTION_9,  ///< in b&b, turn off presolve
  presolve_on         = ENUM_OPTION_10, ///< in b&b, turn on presolve
  heuristics_off      = ENUM_OPTION_11, ///< in b&b, turn off heuristics for finding primal-feasible solutions
  heuristics_on       = ENUM_OPTION_12, ///< in b&b, turn on heuristics for finding primal-feasible solutions
  use_best_bound      = ENUM_OPTION_13, ///< tell b&b solver to use best known obj value bound to prune subtrees
  strong_branching_on = ENUM_OPTION_14, ///< tell b&b solver to use strong branching
}; /* BB_Strategy_Options */

/// @brief Shortcut for checking if a bit is enabled
inline bool use_bb_option(const int strategy, const BB_Strategy_Options option) {
  return strategy & static_cast<int>(option);
}
/// @brief Shortcut for bitwise enabling an option
inline int enable_bb_option(const int strategy, const BB_Strategy_Options option) {
  return strategy | static_cast<int>(option);
}
/// @brief Shortcut for bitwise disabling an option
inline int disable_bb_option(const int strategy, const BB_Strategy_Options option) {
  return strategy & ~static_cast<int>(option);
}
/// @brief Get int value of a set of #BB_Strategy_Options
inline int get_bb_option_value(const std::vector<BB_Strategy_Options>& options) {
  int val = 0;
  for (const auto option : options) {
    val += static_cast<int>(option);
  }
  return val;
} /* get bb_option_value */

/********** DEFINITIONS **********/
//template <class T> class Parameter;
//namespace std {
//  template<class T>
//  struct hash<Parameter<T>> {
//  public:
//    size_t operator()(const Parameter<T> &param) const { return std::hash<std::string>{}(param.name()); };
//  };
//}

/// EnumClassHash is used in order to hash for unordered_map
struct EnumClassHash {
  /// Get size_t of \p t
  template<typename T>
  std::size_t operator()(const T& t) const {
    return static_cast<std::size_t>(t);
  }
};

/// Generic parameter class, with a name and value, ability to check allowable values, and sorting rule
template <class T>
class Parameter {
public:
  /// Set #param_name and #val
  Parameter(std::string name, const T& val) : param_name(name), val(val) {}
  /// Destructor
  virtual ~Parameter() {}
  /// Print parameter (depending on \p amountToPrint)
  virtual std::string to_string(const int amountToPrint = 0, const char* fmt = NULL) const = 0;

  /// Return the #val of the parameter
  virtual const T get() const final { return this->val; }
  /// Return the #name of the parameter
  virtual std::string name() const final { return this->param_name; }
  /// Set parameter to a value
  virtual bool set(const T& val) final { this->val = val; return check(); }
  /// Does the parameter take a valid value? Default: true
  virtual bool check() const { return true; }

  /// Check if one parameter has value coming before the other
  bool operator<(const Parameter& other) const {
    return this->param_name.compare(other.param_name) < 0;
  }
  /// Check if two parameters are equal
  bool operator==(const Parameter& other) const {
    return this->param_name == other.param_name;
  }
protected:
  std::string param_name; ///< name of the parameter
  T val; ///< value of the parameter
}; /* Parameter */

/// Generic parameter of any arithmetic type
template <class T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
//template <class T>
class NumericParameter : public Parameter<T> {
public:
    using Parameter<T>::Parameter;
    /// Constructor setting #param_name, #val, and allowable range of values
    NumericParameter(const std::string& name, const T& val, const T& min_val, const T& max_val)
    : Parameter<T>(name, val), min_val(min_val), max_val(max_val) { check(); }
    /// Constructor setting #param_name, #val, and allowable list of values
    NumericParameter(const std::string& name, const T& val, const std::vector<T>& allowed_vals)
    : Parameter<T>(name, val), allowed_vals(allowed_vals) {
      this->min_val = *std::min_element(allowed_vals.begin(), allowed_vals.end());
      this->max_val = *std::max_element(allowed_vals.begin(), allowed_vals.end());
      check();
    }
    /// Destructor
    virtual ~NumericParameter() {}

    /**
     * @brief Convert parameter to string
     *
     * @details Depending on \p amountToPrint, we can print names, values, or both (comma-separated)
     */
    virtual std::string to_string(
        /// 0 = both names and values,
        /// 1 = only names,
        /// 2 = only values
        const int amountToPrint = 2,
        /// Passed to stringValue()
        const char* fmt = NULL) const {
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

    /// Return #min_val
    virtual const T get_min() const { return this->min_val; }
    /// Return #max_val
    virtual const T get_max() const { return this->max_val; }
    /// Check presence of \p test_val in #allowed_vals (if that is set), or compare it to #min_val and #max_val
    virtual bool val_is_allowed(const T test_val) const {
      if (!this->allowed_vals.empty()) {
        return !(std::find(allowed_vals.begin(), allowed_vals.end(), test_val) == allowed_vals.end());
      } else {
        return !(lessThanVal(this->val, this->min_val) || greaterThanVal(this->val, this->max_val));
      }
    }
    /// Return #allowed_vals
    virtual const std::vector<T>& get_allowed_vals() const { return this->allowed_vals; }
    /// Check if val_is_allowed() and otherwise print error
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
    T min_val; ///< min allowable value
    T max_val; ///< max allowable value
    std::vector<T> allowed_vals; ///< list of allowable values (overwrites #min_val, #max_val)
}; /* NumericParameter */

/// Specialization of NumericParameter to #intParam values
class IntParameter : public NumericParameter<int> {
public:
  using NumericParameter<int>::NumericParameter;
  /// Constructor setting #param_name, #val, allowable range of values, and #param_id
  IntParameter(intParam num, const std::string& name, const int& val, const int& min_val, const int& max_val)
  : NumericParameter<int>(name, val, min_val, max_val), param_id(num) {}
  /// Constructor setting #param_name, #val, allowable list of values, and #param_id
  IntParameter(intParam num, const std::string& name, const int& val, const std::vector<int>& allowed_vals)
  : NumericParameter<int>(name, val, allowed_vals), param_id(num) {}

  virtual std::string to_string(const int amountToPrint = 2, const char* fmt = NULL) const {
    if (fmt) {
      return NumericParameter<int>::to_string(amountToPrint, fmt);
    } else {
      return NumericParameter<int>::to_string(amountToPrint, "%d");
    }
  }

  /// Return #param_id
  virtual const intParam id() const { return this->param_id; }
protected:
  intParam param_id; ///< reference back to the enum intParam
}; /* IntParameter */

/// Specialization of NumericParameter to #doubleParam values
class DoubleParameter : public NumericParameter<double> {
public:
  using NumericParameter<double>::NumericParameter;
  /// Constructor setting #param_name, #val, allowable range of values, and #param_id
  DoubleParameter(doubleParam num, const std::string& name, const double& val, const double& min_val, const double& max_val)
  : NumericParameter<double>(name, val, min_val, max_val), param_id(num) {}
  /// Constructor setting #param_name, #val, allowable list of values, and #param_id
  DoubleParameter(doubleParam num, const std::string& name, const double& val, const std::vector<double>& allowed_vals)
  : NumericParameter<double>(name, val, allowed_vals), param_id(num) {}

  virtual std::string to_string(const int amountToPrint = 2, const char* fmt = NULL) const {
    if (fmt) {
      return NumericParameter<double>::to_string(amountToPrint, fmt);
    } else {
      return NumericParameter<double>::to_string(amountToPrint, "%.3e");
    }
  }

  /// Return #param_id
  virtual const doubleParam id() const { return this->param_id; }
protected:
  doubleParam param_id; ///< reference back to the enum doubleParam
}; /* DoubleParameter */

/// Specialization of NumericParameter to #stringParam values
class StringParameter : public Parameter<std::string> {
public:
  using Parameter<std::string>::Parameter;

  /// Constructor setting #param_name, #val, and #param_id
  StringParameter(stringParam num, const std::string& name, const std::string& val)
  : Parameter<std::string>(name, val), param_id(num) {}

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

  /// Return #param_id
  virtual const stringParam id() const { return this->param_id; }
protected:
  stringParam param_id; ///< reference back to the enum stringParam
}; /* StringParameter */

/********** VPC PARAMETERS STRUCT **********/
/// Class for parameters and carrying the logfile FILE pointer
struct VPCParameters {
  FILE* logfile = NULL; ///< NB: right now this is a shallow copy if this struct gets copied

  /// @brief int parameter values
  /// unordered_map gets printed in reverse order; advantage over map is constant access time on average
  std::unordered_map<intParam, IntParameter, EnumClassHash> intParamValues {
    /// BB_MODE: 010 = branch with vpcs only
    {intParam::BB_MODE,
        IntParameter(intParam::BB_MODE, "BB_MODE",
            10, 0, 111)},
    /// BB_STRATEGY: see BBHelper.hpp; 536 = 000001000011000 => gurobi: 8, user_cuts: 16, presolve_off: 512
    ///previously default was 10776 = 010101000011000 => gurobi: 8, user_cuts: 16, presolve_off: 512, heuristics_off: 2048, use_best_bound: 8192
    {intParam::BB_STRATEGY,
        IntParameter(intParam::BB_STRATEGY, "BB_STRATEGY",
            get_bb_option_value({
              BB_Strategy_Options::gurobi,
              BB_Strategy_Options::user_cuts,
              BB_Strategy_Options::presolve_off
            }),
            std::numeric_limits<int>::min(), std::numeric_limits<int>::max())},
    {intParam::BB_RUNS,
        IntParameter(intParam::BB_RUNS, "BB_RUNS",
            0, std::numeric_limits<int>::min(), std::numeric_limits<int>::max())},
    {intParam::RANDOM_SEED,
        IntParameter(intParam::RANDOM_SEED, "RANDOM_SEED",
            628, -1, std::numeric_limits<int>::max())},
#ifdef TRACE
    {intParam::VERBOSITY,
        IntParameter(intParam::VERBOSITY, "VERBOSITY",
            1, 0, 2)},
#else
    {intParam::VERBOSITY,
        IntParameter(intParam::VERBOSITY, "VERBOSITY",
            0, 0, 2)},
#endif
    {intParam::USE_UNIT_VECTORS,
        IntParameter(intParam::USE_UNIT_VECTORS, "USE_UNIT_VECTORS",
            0, std::numeric_limits<int>::min(), std::numeric_limits<int>::max())},
    {intParam::USE_TIGHT_RAYS,
        IntParameter(intParam::USE_TIGHT_RAYS, "USE_TIGHT_RAYS",
            0, std::numeric_limits<int>::min(), std::numeric_limits<int>::max())},
    {intParam::USE_TIGHT_POINTS,
        IntParameter(intParam::USE_TIGHT_POINTS, "USE_TIGHT_POINTS",
            0, std::numeric_limits<int>::min(), std::numeric_limits<int>::max())},
    {intParam::USE_ITER_BILINEAR,
        IntParameter(intParam::USE_ITER_BILINEAR, "USE_ITER_BILINEAR",
            1, 0, std::numeric_limits<int>::max())},
    {intParam::USE_DISJ_LB,
        IntParameter(intParam::USE_DISJ_LB, "USE_DISJ_LB",
            1, 0, 1)},
    {intParam::USE_ALL_ONES,
        IntParameter(intParam::USE_ALL_ONES, "USE_ALL_ONES",
            1, 0, 1)},
    {intParam::TEMP,
        IntParameter(intParam::TEMP, "TEMP",
            0, std::numeric_limits<int>::min(), std::numeric_limits<int>::max())},
    {intParam::STRENGTHEN,
        IntParameter(intParam::STRENGTHEN, "STRENGTHEN",
            1, 0, 2)},
    {intParam::ROUNDS,
        IntParameter(intParam::ROUNDS, "ROUNDS",
            1, 0, std::numeric_limits<int>::max())},
    {intParam::PRLP_FLIP_BETA,
        IntParameter(intParam::PRLP_FLIP_BETA, "PRLP_FLIP_BETA",
            0, -1, 1)},
    {intParam::PREPROCESS,
            IntParameter(intParam::PREPROCESS, "PREPROCESS",
                0, 0, 2)},
    {intParam::PARTIAL_BB_NUM_STRONG,
        IntParameter(intParam::PARTIAL_BB_NUM_STRONG, "PARTIAL_BB_NUM_STRONG",
            5, std::numeric_limits<int>::min(), std::numeric_limits<int>::max())},
    /// PARTIAL_BB_STRATEGY: 004 => default variable decision, default branch decision, objective-based node comparison
    {intParam::PARTIAL_BB_STRATEGY,
        IntParameter(intParam::PARTIAL_BB_STRATEGY, "PARTIAL_BB_STRATEGY",
            4, std::numeric_limits<int>::min(), std::numeric_limits<int>::max())},
    {intParam::MODE,
        IntParameter(intParam::MODE, "MODE",
            0, {0, 1, 3, 4})},
    {intParam::GOMORY,
        IntParameter(intParam::GOMORY, "GOMORY",
            0, -1, 1)},
    {intParam::DISJ_TERMS,
        IntParameter(intParam::DISJ_TERMS, "DISJ_TERMS",
            0, std::numeric_limits<int>::min(), std::numeric_limits<int>::max())},
    {intParam::CUTLIMIT,
        IntParameter(intParam::CUTLIMIT, "CUTLIMIT",
            -1, std::numeric_limits<int>::min(), std::numeric_limits<int>::max())},
  }; /* intParamValues */

  /// @brief double parameter values
  std::unordered_map<doubleParam, DoubleParameter, EnumClassHash> doubleParamValues {
    {doubleParam::TIMELIMIT,
      DoubleParameter(doubleParam::TIMELIMIT, "TIMELIMIT",
          60, 0., std::numeric_limits<double>::max())},
    {doubleParam::PRLP_TIMELIMIT,
        DoubleParameter(doubleParam::PRLP_TIMELIMIT, "PRLP_TIMELIMIT",
            -1, -1., std::numeric_limits<double>::max())},
    {doubleParam::PARTIAL_BB_TIMELIMIT,
        DoubleParameter(doubleParam::PARTIAL_BB_TIMELIMIT, "PARTIAL_BB_TIMELIMIT",
            3600, 0., std::numeric_limits<double>::max())},
    {doubleParam::MIN_ORTHOGONALITY,
        DoubleParameter(doubleParam::MIN_ORTHOGONALITY, "MIN_ORTHOGONALITY",
            0., 0., 1.)},
    {doubleParam::IP_OBJ,
        DoubleParameter(doubleParam::IP_OBJ, "IP_OBJ",
            std::numeric_limits<double>::max(), std::numeric_limits<double>::lowest(), std::numeric_limits<double>::max())},
    {doubleParam::INF,
        DoubleParameter(doubleParam::EPS, "INF",
            std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max())},
    {doubleParam::EPS,
        DoubleParameter(doubleParam::EPS, "EPS",
            1e-7, 0., 1.)},
    {doubleParam::BB_TIMELIMIT,
        DoubleParameter(doubleParam::BB_TIMELIMIT, "BB_TIMELIMIT",
            3600., 0., std::numeric_limits<double>::max())},
  }; /* doubleParamValues */

  /// @brief string parameter values
  std::unordered_map<stringParam, StringParameter, EnumClassHash> stringParamValues {
    {stringParam::DISJ_OPTIONS, StringParameter("DISJ_OPTIONS", "")},
    {stringParam::SOLFILE,
      StringParameter(stringParam::SOLFILE, "SOLFILE", "")},
    {stringParam::OPTFILE,
      StringParameter(stringParam::OPTFILE, "OPTFILE", "")},
    {stringParam::LOGFILE,
        StringParameter(stringParam::LOGFILE, "LOGFILE", "")},
    {stringParam::FILENAME,
        StringParameter(stringParam::FILENAME, "FILENAME", "")},
  }; /* stringParamValues */

  /// @brief Integer constants
  std::unordered_map<intConst, IntParameter, EnumClassHash> intConstValues {
    {intConst::PRLP_PRESOLVE, IntParameter("PRLP_PRESOLVE", 2, 2, 2)}, // use presolve when solving PRLP (either initial or resolve)
    {intConst::NB_SPACE, IntParameter("NB_SPACE", 1, 1, 1)}, // currently only works with true
    {intConst::NUM_OBJ_PER_POINT, IntParameter("NUM_OBJ_PER_POINT", -2, -2, -2)}, // sqrt(n)
    {intConst::MODE_OBJ_PER_POINT, IntParameter("MODE_OBJ_PER_POINT", 121, 121, 121)}, // one ray at a time, points+rays+variables, large to small angle (descending)
    {intConst::MAX_SUPPORT_ABS, IntParameter("MAX_SUPPORT_ABS", std::numeric_limits<int>::max(), std::numeric_limits<int>::max(), std::numeric_limits<int>::max())},
    {intConst::LUB, IntParameter("LUB", 1e3, 1e3, 1e3)},
    {intConst::CHECK_DUPLICATES, IntParameter("CHECK_DUPLICATES", 1, 1, 1)},
  }; /* intConstValues */

  /// @brief Double constants
  std::unordered_map<doubleConst, DoubleParameter, EnumClassHash> doubleConstValues {
    {doubleConst::MAX_SUPPORT_REL, DoubleParameter("MAX_SUPPORT_REL", 0.9, 0.9, 0.9)},
    {doubleConst::MAX_DYN_LUB, DoubleParameter("MAX_DYN_LUB", 1e13, 1e13, 1e13)},
    {doubleConst::MAX_DYN, DoubleParameter("MAX_DYN", 1e8, 1e8, 1e8)},
    {doubleConst::MIN_VIOL_REL, DoubleParameter("MIN_VIOL_REL", 1e-7, 1e-7, 1e-7)},
    {doubleConst::MIN_VIOL_ABS, DoubleParameter("MIN_VIOL_ABS", 1e-7, 1e-7, 1e-7)},
    {doubleConst::EPS_COEFF_LUB, DoubleParameter("EPS_COEFF_LUB", 1e-13, 1e-13, 1e-13)},
    {doubleConst::EPS_COEFF, DoubleParameter("EPS_COEFF", 1e-5, 1e-5, 1e-5)},
    {doubleConst::MIN_PRLP_TIMELIMIT, DoubleParameter("MIN_PRLP_TIMELIMIT", 5., 5., 5.)},
    {doubleConst::RAYEPS, DoubleParameter("RAYEPS", 1e-7, 1e-7, 1e-7)},
    {doubleConst::DIFFEPS, DoubleParameter("DIFFEPS", 1e-3, 1e-3, 1e-3)}, // to check whether something is different enough to throw an error
    {doubleConst::AWAY, DoubleParameter("AWAY", 1e-3, 1e-3, 1e-3)},
  }; /* doubleConstValues */

  ///@{
  /// @name Constructors

  /// @brief Default constructor
  VPCParameters() { }

  /// @brief Copy constructor
  VPCParameters(const VPCParameters& source) {
    this->logfile = source.logfile;
    this->intParamValues = source.intParamValues;
    this->doubleParamValues = source.doubleParamValues;
    this->stringParamValues = source.stringParamValues;
    this->intConstValues = source.intConstValues;
    this->doubleConstValues = source.doubleConstValues;
  } /* copy constructor */
  ///@}

  ///@{
  /// @name Name methods

  /// @brief Return name of intParam \p param
  std::string name(intParam param) const { return intParamValues.find(param)->second.name(); }
  /// @brief Return name of doubleParam \p param
  std::string name(doubleParam param) const { return doubleParamValues.find(param)->second.name(); }
  /// @brief Return name of stringParam \p param
  std::string name(stringParam param) const { return stringParamValues.find(param)->second.name(); }
  /// @brief Return name of intConst \p param
  std::string name(intConst param) const { return intConstValues.find(param)->second.name(); }
  /// @brief Return name of doubleConst \p param
  std::string name(doubleConst param) const { return doubleConstValues.find(param)->second.name(); }
  ///@}

  ///@{
  /// @name Get / set methods

  /// @brief Get value of intParam \p param
  int         get(intParam param) const { return intParamValues.find(param)->second.get(); }
  /// @brief Get value of doubleParam \p param
  double      get(doubleParam param) const { return doubleParamValues.find(param)->second.get(); }
  /// @brief Get value of stringParam \p param
  std::string get(stringParam param) const { return stringParamValues.find(param)->second.get(); }
  /// @brief Get value of intConst \p param
  int         get(intConst param) const { return intConstValues.find(param)->second.get(); }
  /// @brief Get value of doubleConst \p param
  double      get(doubleConst param) const { return doubleConstValues.find(param)->second.get(); }

  /// @brief Set intParam \p param to \p value
  void        set(intParam param, const int value) { intParamValues.find(param)->second.set(value); }
  /// @brief Set doubleParam \p param to \p value
  void        set(doubleParam param, const double value) { doubleParamValues.find(param)->second.set(value); }
  /// @brief Set stringParam \p param to \p value
  void        set(stringParam param, std::string value) { stringParamValues.find(param)->second.set(value); }

  /// Set parameters by name
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
  ///@}
}; /* struct VPCParameters */

/// @brief Read parameters and constants from a given file with the format that each line has param_name,param_value (comma-separated)
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

/// @brief Print parameters and constants
inline void printParams(
    /// parameters to print
    const VPCParameters& params, 
    /// where to print
    FILE* logfile = stdout, 
    /// \li 0 = numeric param name/values except string params (newline-separated)
    /// \li 1 = only numeric param names (comma-separated)
    /// \li 2 = only numeric param values (comma-separated)
    /// \li 3 = string param name/values (newline-separated)
    /// \li 4 = only string names (comma-separated)
    /// \li 5 = only string values (comma-separated)
    /// \li 6 = numeric+string param name/values (newline-separated)
    /// \li 7 = only numeric+string names (comma-separated)
    /// \li 8 = only numeric+string values (comma-separated)
    /// \li 9 = all const name/values (newline-separated)
    /// \li 10 = only const names (comma-separated)
    /// \li 11 = only const values (comma-separated)
    const int whatToPrint = 0) {
  if (!logfile)
    return;

//  !(std::find(allowed_vals.begin(), allowed_vals.end(), test_val) == allowed_vals.end())
  const bool printIntParams = whatToPrint == 0 || whatToPrint == 1 || whatToPrint == 2 || whatToPrint == 6 || whatToPrint == 7 || whatToPrint == 8;
  const bool printDoubleParams = whatToPrint == 0 || whatToPrint == 1 || whatToPrint == 2 || whatToPrint == 6 || whatToPrint == 7 || whatToPrint == 8;
  const bool printStringParams = whatToPrint == 3 || whatToPrint == 4 || whatToPrint == 5 || whatToPrint == 6 || whatToPrint == 7 || whatToPrint == 8;
  const bool printIntConsts = whatToPrint == 9 || whatToPrint == 10 || whatToPrint == 11;
  const bool printDoubleConsts = whatToPrint == 9 || whatToPrint == 10 || whatToPrint == 11;
  //const bool printNames = whatToPrint == 0 || whatToPrint == 1 || whatToPrint == 3 || whatToPrint == 4 || whatToPrint == 6 || whatToPrint == 7 || whatToPrint == 9 || whatToPrint == 10;
  //const bool printValues = whatToPrint == 0 || whatToPrint == 2 || whatToPrint == 3 || whatToPrint == 5 || whatToPrint == 6 || whatToPrint == 8 || whatToPrint == 9 || whatToPrint == 11;
  const char SEP = (whatToPrint == 0 || whatToPrint == 3 || whatToPrint == 6 || whatToPrint == 10) ? '\n' : ',';
  const int amountToPrint = whatToPrint % 3;

  if (printIntParams)
  for (auto param : params.intParamValues) {
    if (param.second.id() == intParam::DISJ_TERMS) {
      // Special handling for disjunctive terms to sort well in the logfile
      fprintf(logfile, "%s%c", param.second.to_string(amountToPrint, "%05d").c_str(), SEP);
    } else {
      fprintf(logfile, "%s%c", param.second.to_string(amountToPrint).c_str(), SEP);
    }
  }
  if (printDoubleParams)
  for (auto param : params.doubleParamValues) {
    fprintf(logfile, "%s%c", param.second.to_string(amountToPrint).c_str(), SEP);
  }
  if (printStringParams)
  for (auto param : params.stringParamValues) {
    fprintf(logfile, "%s%c", param.second.to_string(amountToPrint).c_str(), SEP);
  }
  if (printIntConsts)
  for (auto param : params.intConstValues) {
    fprintf(logfile, "%s%c", param.second.to_string(amountToPrint).c_str(), SEP);
  }
  if (printDoubleConsts)
  for (auto param : params.doubleConstValues) {
    fprintf(logfile, "%s%c", param.second.to_string(amountToPrint).c_str(), SEP);
  }
  fflush(logfile);
} /* printParams */
} /* namespace VPCParameters */
