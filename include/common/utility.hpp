/**
 * @file utility.hpp
 * @brief Utility functions
 *
 * The header is going to be included in many places,
 * so ideally it will not have too much included that is unnecessary
 *
 * @author A. M. Kazachkov
 * @date 2018-Dec-24
 */
#pragma once

// standard library
#include <cstdio>
#include <iostream> // cerr
#include <memory> // shared_ptr
#include <string>
#include <vector>
#include <cmath> // abs
#include <limits> // numeric_limits

// coin-or
#include "SolverInterface.hpp" // SolverInterface

// VPC
#include <OsiSolverInterface.hpp>

class CoinPackedVectorBase;
class CoinPackedVector;
class CoinPackedMatrix;

#define macro_to_string(s) #s
#define x_macro_to_string(s) macro_to_string(s)

/// Templated format string
template <typename T> inline constexpr const char *outFormat = "%g"; // pick a default.
template <> inline constexpr const char *outFormat<int> = "%d";
template <> inline constexpr const char *outFormat<float> = "%g";
template <> inline constexpr const char *outFormat<double> = "%lf";
template <> inline constexpr const char *outFormat<long double> = "%Lf";
template <> inline constexpr const char *outFormat<char> = "%c";

/// @brief overload << operator for vector of ints
inline std::ostream& operator<<(std::ostream& os, const std::vector<int> &input) {
  os << "{";
  for (unsigned tmp_i = 0; tmp_i < input.size(); tmp_i++) {
    int i = input[tmp_i];
    os << i << ((tmp_i < input.size() - 1) ? ", " : "");
  }
  os << "}";
  return os;
} /* overload << operator for vector of ints */

/// @brief overload << operator for vector of doubles
inline std::ostream& operator<<(std::ostream& os, const std::vector<double> &input) {
  os << "{";
  for (unsigned tmp_i = 0; tmp_i < input.size(); tmp_i++) {
    double i = input[tmp_i];
    os << i << ((tmp_i < input.size() - 1) ? ", " : "");
  }
  os << "}";
  return os;
} /* overload << operator for vector of doubles */

///@{
/// @name Sorting functions
// @defgroup Sorting Sorting functions
/// @details NOTE: the below sorting functions will only work when the index list given is {0,...,size}
/// (not necessarily in that order)
/// because the index list is used to access elements of arr
///
/// If you want to, say, sort [10,1,2] by [0.5,20.5,10.0], the below would not work
/// Typically you need to use a struct then combining the data together

// @addtogroup Sorting
/// @brief In descending order
template<class T> struct index_cmp_dsc {

  const T arr; ///< array used to be sorted

  /// Constructor
  index_cmp_dsc(const T _arr) :
      arr(_arr) {
  }

  /// Comparator
  bool operator()(const size_t a, const size_t b) const {
    //return arr[a] > arr[b] + param.getEPS();
    return arr[a] > arr[b];
  }
};

// @addtogroup Sorting
/// @brief In ascending order
template<class T> struct index_cmp_asc {

  const T arr; ///< array used to be sorted

  /// Constructor
  index_cmp_asc(const T _arr) :
      arr(_arr) {
  }

  /// Comparator
  bool operator()(const size_t a, const size_t b) const {
    //return arr[a] < arr[b] - param.getEPS();
    return arr[a] < arr[b];
  }
};
///@}

/* #define error_msg(str, fmt) \
  char str[500]; \
  snprintf(str, sizeof(str) / sizeof(char), "*** ERROR: %s:%d: " fmt, __FILE__, __LINE__); \
  std::cerr << str

#define warning_msg(str, fmt) \
  char str[500]; \
  snprintf(str, sizeof(str) / sizeof(char), "*** WARNING: %s:%d: " fmt, __FILE__, __LINE__); \
  std::cerr << str */

#ifdef VPC_VERSION
#define error_msg(str, fmt, ...) \
  char str[2028]; \
  snprintf(str, sizeof(str) / sizeof(char), "*** ERROR (version %.8s): %s:%d: " fmt, x_macro_to_string(VPC_VERSION), __FILE__, __LINE__, ##__VA_ARGS__); \
  std::cerr << str
#else
#define error_msg(str, fmt, ...) \
  char str[2028]; \
  snprintf(str, sizeof(str) / sizeof(char), "*** ERROR: %s:%d: " fmt, __FILE__, __LINE__, ##__VA_ARGS__); \
  std::cerr << str
#endif

#ifdef VPC_VERSION
#define warning_msg(str, fmt, ...) \
  char str[2028]; \
  snprintf(str, sizeof(str) / sizeof(char), "*** WARNING (version %.8s): %s:%d: " fmt, x_macro_to_string(VPC_VERSION), __FILE__, __LINE__, ##__VA_ARGS__); \
  std::cerr << str
#else
#define warning_msg(str, fmt, ...) \
  char str[2028]; \
  snprintf(str, sizeof(str) / sizeof(char), "*** WARNING: %s:%d: " fmt, __FILE__, __LINE__, ##__VA_ARGS__); \
  std::cerr << str
#endif

/**
 * @brief Writes error and closes file myfile.
 */
inline void writeErrorToLog(std::string text, FILE *myfile) {
  if (myfile == NULL || myfile == stdout || myfile == stderr)
    return;
  fprintf(myfile, "||%c%s", ',', text.c_str());
  fclose(myfile);
}

/** @brief Create temporary filename */
void createTmpFilename(std::string& f_name, const std::string add_ext = "");

/** @brief Separate filename into the directory, instance name, and extension */
int parseFilename(std::string& dir, std::string& instname, std::string& in_file_ext, const std::string& fullfilename, FILE* logfile);

/// @brief Get objective value from file \p opt_filename where each line is "instance,value"
double getObjValueFromFile(std::string opt_filename, std::string fullfilename, FILE* logfile);

/// @brief Retrieve solution from a file \p filename.
void getSolFromFile(const char* filename, std::vector<double>& sol);

/// @brief Check if a file exists
bool fexists(const char* filename);

/// @brief Parses int from string using strtol
bool parseInt(const char *str, int &val);

/// @brief Parses long int from string using strtol
bool parseLong(const char *str, long &val);

/// @brief Parses double from string using strtod
bool parseDouble(const char *str, double &val);

/// @brief Converts string to lower case
std::string lowerCaseString(const std::string& data);
/// @brief Convert vector of strings to lower case
std::vector<std::string> lowerCaseStringVector(const std::vector<std::string>& strVec);
/// @brief Converts string to upper case
std::string upperCaseString(const std::string& tmpData);
/// @brief Converts string to upper case and remove underscores (for processing parameters)
std::string upperCaseStringNoUnderscore(const std::string& tmpData);

/// @brief Check equality within tolerance \p eps
inline bool isVal(const double val1, const double val2, const double eps = 1e-7) {
  return (std::abs((val1) - (val2)) <= eps);
} /* isVal */

/// @brief Check equality to zero within tolerance \p eps; calls #isVal
inline bool isZero(const double val, const double eps = 1e-7) {
  return isVal(val, 0.0, eps);
} /* isZero */

/// @brief Check whether \p val1 is at least \p eps less than \p val2
inline bool lessThanVal(const double val1, const double val2, const double eps = 1e-7) {
  return (val1 < (val2 - eps));
} /* lessThanVal */

/// @brief Check whether \p val1 is at least \p eps greater than \p val2
inline bool greaterThanVal(const double val1, const double val2, const double eps = 1e-7) {
  return (val1 > (val2 + eps));
} /* greaterThanVal */

/// @brief Global infinity
inline double getInfinity() {
  return std::numeric_limits<double>::max();
}

/// @brief Check whether value is at infinity using #lessThanVal
inline bool isInfinity(const double val, const double infinity = __DBL_MAX__, const double eps = 1e-7) {
  return !lessThanVal(val, infinity, eps);
} /* isInfinity */

/// @brief Check whether value is at negative infinity using #greaterThanVal
inline bool isNegInfinity(const double val, const double neg_infinity = __DBL_MIN__, const double eps = 1e-7) {
  return !greaterThanVal(val, neg_infinity, eps);
} /* isInfinity */

///// @brief Covert an integer into a string, accounting for possible infinite values
//inline const std::string stringValue(const int value,
//    const char* format = "%d",
//    const int MIN = std::numeric_limits<int>::min(),
//    const int MAX = std::numeric_limits<int>::max()) {
//  if (!lessThanVal(value, MAX)) {
//    const std::string infty = "\'inf\'";
//    return infty;
//  } else if (!greaterThanVal(value, MIN)) {
//    const std::string neg_infty = "\'-inf\'";
//    return neg_infty;
//  } else {
//    char temp[500];
//    snprintf(temp, sizeof(temp) / sizeof(char), format, value);
//    std::string tmp(temp);
//    return tmp;
//  }
//} /* stringValue (int) */
//
///// @brief Covert a long into a string, accounting for possible infinite values
//inline const std::string stringValue(const long value,
//    const char* format = "%ld",
//    const long MIN = std::numeric_limits<long>::min(),
//    const long MAX = std::numeric_limits<long>::max()) {
//  if (!lessThanVal(value, MAX)) {
//    const std::string infty = "\'inf\'";
//    return infty;
//  } else if (!greaterThanVal(value, MIN)) {
//    const std::string neg_infty = "\'-inf\'";
//    return neg_infty;
//  } else {
//    char temp[500];
//    snprintf(temp, sizeof(temp) / sizeof(char), format, value);
//    std::string tmp(temp);
//    return tmp;
//  }
//} /* stringValue (long) */

/// @brief Covert a numeric type into a string, accounting for possible infinite values
template <class T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
inline const std::string stringValue(
    /// value to be converted
    const T value,
    /// formatting to be used (note that outFormat<T> needs to be defined for the specific type T)
    const char* format = outFormat<T>,
    /// value of infinity to be used
    const double INF = std::numeric_limits<T>::max(),
    /// how many digits before the decimal to print (-1 default implies do not limit)
    const int NUM_DIGITS_BEFORE_DEC = -1,
    /// how many digits after the decimal to print (-1 default implies do not limit)
    const int NUM_DIGITS_AFTER_DEC = -1) {
  char temp[500];
  if (!lessThanVal(value, INF)) {
    if (NUM_DIGITS_BEFORE_DEC == -1) {
      snprintf(temp, sizeof(temp) / sizeof(char), "%s", "\'inf\'");
    } else {
      snprintf(temp, sizeof(temp) / sizeof(char), "%-*s", NUM_DIGITS_BEFORE_DEC, "\'inf\'");
    }
  } else if (!greaterThanVal(value, -INF)) {
    if (NUM_DIGITS_BEFORE_DEC == -1) {
      snprintf(temp, sizeof(temp) / sizeof(char), "%s", "\'-inf\'");
    } else {
      snprintf(temp, sizeof(temp) / sizeof(char), "%-*s", NUM_DIGITS_BEFORE_DEC, "\'-inf\'");
    }
  } else {
    if (NUM_DIGITS_BEFORE_DEC == -1 && NUM_DIGITS_AFTER_DEC == -1) {
      snprintf(temp, sizeof(temp) / sizeof(char), format, value);
    } else if (NUM_DIGITS_BEFORE_DEC >= 0 && NUM_DIGITS_AFTER_DEC >= 0) {
      snprintf(temp, sizeof(temp) / sizeof(char), format, NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC, value);
    } else if (NUM_DIGITS_BEFORE_DEC >= 0) {
      snprintf(temp, sizeof(temp) / sizeof(char), format, NUM_DIGITS_BEFORE_DEC, value);
    } else {
      snprintf(temp, sizeof(temp) / sizeof(char), format, NUM_DIGITS_AFTER_DEC, value);
    }
  }
  std::string tmp(temp);
  return tmp;
} /* stringValue (numeric) */

/// @brief Currently returns \p value; mostly for ease of use with params; in the future we may way to format strings though
inline const std::string stringValue(const std::string value, const char* format = "%s") {
  return value;
} /* stringValue (string) */

/// @brief Calls #dotProductNoCompensation(int, const int*, const double*, const double*)
double dotProductNoCompensation(const CoinPackedVector& vec1, const double* vec2);
/// @brief Calls #dotProductNoCompensation(int, const int*, const double*, const int*, const double*)
double dotProductNoCompensation(const CoinPackedVector& vec1, const CoinPackedVector& vec2);

/// @brief Compute dot product using compensated summation to have small numerical error.
/// First version: dense vectors.
double dotProductNoCompensation(const double* a, const double* b, int dimension);

/// @brief Second version: first vector is sparse, second one is dense.
double dotProductNoCompensation(int sizea, const int* indexa, const double* a,
    const double* b);

/// @brief Third version: sparse vectors.
double dotProductNoCompensation(int sizea, const int* indexa, const double* a, int sizeb,
    const int* indexb, const double* b);

/**********************************************************
 * The following code is for dot product with compensation,
 * taken from the source cited below.
 **********************************************************/
// Name:     common_definitions.hpp
// Author:   Giacomo Nannicini
//           Singapore University of Technology and Design, Singapore
//           email: nannicini@sutd.edu.sg
// Date:     04/09/11
// Last edit: 03/27/12
//-----------------------------------------------------------------------------
// Copyright (C) 2012, Giacomo Nannicini.  All Rights Reserved.

/// @brief Dot product between sparse and dense vec
double dotProduct(const CoinPackedVector& vec1, const double* vec2);

/// @brief Dot product between two sparse vectors
double dotProduct(const CoinPackedVector& vec1, const CoinPackedVector& vec2);

/// @brief Compute dot product using compensated summation to have small
/// numerical error. First version: dense vectors
double dotProduct(const double* a, const double* b, int dimension);

/// @brief Second version: first vector is sparse, second one is dense
double dotProduct(int sizea, const int* indexa, const double* a,
    const double* b);

/// @brief Third version: sparse vectors
double dotProduct(int sizea, const int* indexa, const double* a, int sizeb,
    const int* indexb, const double* b);

double getRowTwoNorm(const int row, const CoinPackedMatrix* const mat);
void packedSortedVectorSum(CoinPackedVector& sum, const double mult1,
    const CoinPackedVectorBase& vec1, const double mult2,
    const CoinPackedVectorBase& vec2, const double eps);

bool variableBoundsContained(const OsiSolverInterface* const solver1,
                             const OsiSolverInterface* const solver2);

/// @brief create a mutable solver interface
std::shared_ptr<SolverInterface> getSolver(const OsiSolverInterface* const si, FILE* logfile = NULL);
