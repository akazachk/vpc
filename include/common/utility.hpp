/**
 * @file utility.hpp
 * @author A. M. Kazachkov
 * @date 2018-Dec-24
 * @brief Utility functions
 * @details The header is going to be included in many places,
 * so ideally it will not have too much included that is unnecessary
 */
#pragma once

#include <cstdio>
#include <iostream> // cerr
#include <string>
#include <vector>
#include <cmath> // abs
#include <limits> // numeric_limits

class CoinPackedVectorBase;
class CoinPackedVector;
class CoinPackedMatrix;

#define macro_to_string(s) #s
#define x_macro_to_string(s) macro_to_string(s)

inline std::ostream& operator<<(std::ostream& os, const std::vector<int> &input) {
  os << "{";
  for (unsigned tmp_i = 0; tmp_i < input.size(); tmp_i++) {
    int i = input[tmp_i];
    os << i << ((tmp_i < input.size() - 1) ? ", " : "");
  }
  os << "}";
  return os;
} /* overload << operator for vector of ints */

inline std::ostream& operator<<(std::ostream& os, const std::vector<double> &input) {
  os << "{";
  for (unsigned tmp_i = 0; tmp_i < input.size(); tmp_i++) {
    double i = input[tmp_i];
    os << i << ((tmp_i < input.size() - 1) ? ", " : "");
  }
  os << "}";
  return os;
} /* overload << operator for vector of doubles */

/**
 * @brief For sorting
 *
 * @details NOTE: the below sorting functions will only work when the index list given is {0,...,size}
 * (not necessarily in that order)
 * because the index list is used to access elements of arr
 *
 * If you want to, say, sort [10,1,2] by [0.5,20.5,10.0], the below would not work
 * Typically you need to use a struct then combining the data together
 */
template<class T> struct index_cmp_dsc {

  const T arr; // array used to be sorted

  // Constructor
  index_cmp_dsc(const T _arr) :
      arr(_arr) {
  }

  // Comparator
  bool operator()(const size_t a, const size_t b) const {
    //return arr[a] > arr[b] + param.getEPS();
    return arr[a] > arr[b];
  }
};

template<class T> struct index_cmp_asc {

  const T arr; // array used to be sorted

  // Constructor
  index_cmp_asc(const T _arr) :
      arr(_arr) {
  }

  // Comparator
  bool operator()(const size_t a, const size_t b) const {
    //return arr[a] < arr[b] - param.getEPS();
    return arr[a] < arr[b];
  }
};

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

double getObjValueFromFile(std::string opt_filename, std::string fullfilename, FILE* logfile);

/**
  * @brief Check if a file exists
  */
bool fexists(const char* filename);

/**
 * @brief Parses int from string using strtol
 */
bool parseInt(const char *str, int &val);

/**
 * @brief Parses long int from string using strtol
 */
bool parseLong(const char *str, long &val);

/**
 * @brief Parses double from string using strtod
 */
bool parseDouble(const char *str, double &val);

std::string lowerCaseString(const std::string& data);
std::vector<std::string> lowerCaseStringVector(const std::vector<std::string>& strVec);
std::string upperCaseString(const std::string& tmpData);
std::string upperCaseStringNoUnderscore(const std::string& tmpData);

inline bool isVal(const double val1, const double val2, const double eps = 1e-7) {
  return (std::abs((val1) - (val2)) <= eps);
} /* isVal */

inline bool isZero(const double val, const double eps = 1e-7) {
  return isVal(val, 0.0, eps);
} /* isZero */

inline bool lessThanVal(const double val1, const double val2, const double eps = 1e-7) {
  return (val1 < (val2 - eps));
} /* lessThanVal */

inline bool greaterThanVal(const double val1, const double val2, const double eps = 1e-7) {
  return (val1 > (val2 + eps));
} /* greaterThanVal */

inline double getInfinity() {
  return std::numeric_limits<double>::max();
}

inline bool isInfinity(const double val, const double infinity = __DBL_MAX__, const double eps = 1e-7) {
  return !lessThanVal(val, infinity, eps);
} /* isInfinity */

inline bool isNegInfinity(const double val, const double neg_infinity = __DBL_MIN__, const double eps = 1e-7) {
  return !greaterThanVal(val, neg_infinity, eps);
} /* isInfinity */

inline const std::string stringValue(const int value,
    const char* format = "%d") {
  if (!lessThanVal(value, std::numeric_limits<int>::max())) {
    const std::string infty = "\'inf\'";
    return infty;
  } else if (!greaterThanVal(value, std::numeric_limits<int>::min())) {
    const std::string neg_infty = "\'-inf\'";
    return neg_infty;
  } else {
    char temp[500];
    snprintf(temp, sizeof(temp) / sizeof(char), format, value);
    std::string tmp(temp);
    return tmp;
  }
} /* stringValue (int) */

inline const std::string stringValue(const long value,
    const char* format = "%ld") {
  if (!lessThanVal(value, std::numeric_limits<long>::max())) {
    const std::string infty = "\'inf\'";
    return infty;
  } else if (!greaterThanVal(value, std::numeric_limits<long>::min())) {
    const std::string neg_infty = "\'-inf\'";
    return neg_infty;
  } else {
    char temp[500];
    snprintf(temp, sizeof(temp) / sizeof(char), format, value);
    std::string tmp(temp);
    return tmp;
  }
} /* stringValue (long) */

inline const std::string stringValue(const double value, const char* format = "%f",
    const int NUM_DIGITS_BEFORE_DEC = -1, const int NUM_DIGITS_AFTER_DEC = -1) {
  if (!lessThanVal(value, std::numeric_limits<double>::max())) {
    const std::string infty = "\'inf\'";
    return infty;
  } else if (!greaterThanVal(value, std::numeric_limits<double>::lowest())) {
    const std::string neg_infty = "\'-inf\'";
    return neg_infty;
  } else {
    char temp[500];
    if (NUM_DIGITS_BEFORE_DEC == -1 && NUM_DIGITS_AFTER_DEC == -1) {
      snprintf(temp, sizeof(temp) / sizeof(char), format, value);
    } else if (NUM_DIGITS_BEFORE_DEC >= 0 && NUM_DIGITS_AFTER_DEC >= 0) {
      snprintf(temp, sizeof(temp) / sizeof(char), format, NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC, value);
    } else if (NUM_DIGITS_BEFORE_DEC >= 0) {
      snprintf(temp, sizeof(temp) / sizeof(char), format, NUM_DIGITS_BEFORE_DEC, value);
    } else {
      snprintf(temp, sizeof(temp) / sizeof(char), format, NUM_DIGITS_AFTER_DEC, value);
    }
    std::string tmp(temp);
    return tmp;
  }
} /* stringValue (double) */

/** @brief The below is mostly for ease of use with params; in the future we may way to format strings though */
inline const std::string stringValue(const std::string value, const char* format = "%s") {
  return value;
} /* stringValue (string) */

double dotProductNoCompensation(const CoinPackedVector& vec1, const double* vec2);
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
