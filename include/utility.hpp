#pragma once

#include <cstdio>
#include <iostream> // cerr
#include <string>
#include <vector>
#include <cmath> // abs
#include <limits> // numeric_limits

// COIN-OR includes
#include <CoinPackedVector.hpp>
#include <OsiSolverInterface.hpp>

// VPC includes
#include "VPCParameters.hpp"

/* #define error_msg(str, fmt) \
  char str[500]; \
  snprintf(str, sizeof(str) / sizeof(char), "*** ERROR: %s:%d: " fmt, __FILE__, __LINE__); \
  std::cerr << str

#define warning_msg(str, fmt) \
  char str[500]; \
  snprintf(str, sizeof(str) / sizeof(char), "*** WARNING: %s:%d: " fmt, __FILE__, __LINE__); \
  std::cerr << str */

#define error_msg(str, fmt, ...) \
  char str[500]; \
  snprintf(str, sizeof(str) / sizeof(char), "*** ERROR: %s:%d: " fmt, __FILE__, __LINE__, ##__VA_ARGS__); \
  std::cerr << str

#define warning_msg(str, fmt, ...) \
  char str[500]; \
  snprintf(str, sizeof(str) / sizeof(char), "*** WARNING: %s:%d: " fmt, __FILE__, __LINE__, ##__VA_ARGS__); \
  std::cerr << str

/***********************************************************************/
/**
 * @brief Writes error and closes file myfile.
 */
inline void writeErrorToLog(std::string text, FILE *myfile) {
  if (myfile == NULL || myfile == stdout)
    return;
  fprintf(myfile, "||%c%s", ',', text.c_str());
  fclose(myfile);
}

double getObjValueFromFile(const VPCParameters& params);

void setCompNBCoorPoint(CoinPackedVector& vec, double& objViolation,
    const VPCParameters& params, const OsiSolverInterface* const tmpSolver,
    const OsiSolverInterface* const origSolver,
    const std::vector<int>& nonBasicVarIndex,
    const std::vector<double>& nonBasicReducedCost, const int deletedVar = -1);
void setCompNBCoorRay(CoinPackedVector& vec, const double* ray, double& objViolation, double& scale,
    const VPCParameters& params, const OsiSolverInterface* const tmpSolver,
    const OsiSolverInterface* const origSolver,
    const std::vector<int>& rowOfOrigNBVar,
    const std::vector<int>& nonBasicVarIndex,
    const std::vector<double>& nonBasicReducedCost,
    const int tmpNBVar, const int deletedVar = -1,
    const bool rayNeedsCalculation = false);
void setCompNBCoor(CoinPackedVector& vec, double& objViolation,
    const VPCParameters& params, const double* const currColValue,
    const double* const currSlackValue,
    const OsiSolverInterface* const origSolver,
    const std::vector<int>& nonBasicVarIndex,
    const std::vector<double>& nonBasicReducedCost, const int deletedVar = -1);
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

inline bool isInfinity(const double val, const double infinity = __DBL_MAX__, const double eps = 1e-7) {
  return (val >= (infinity - eps));
} /* isInfinity */

inline const std::string stringValue(const int value,
    const char* format = "%d") {
  if (!isInfinity(std::abs(value), std::numeric_limits<int>::max())) {
    char temp[500];
    snprintf(temp, sizeof(temp) / sizeof(char), format, value);
    std::string tmp(temp);
    return tmp;
  } else {
    if (greaterThanVal(value, 0.0)) {
      const std::string infty = "\'inf\'";
      return infty;
    } else {
      const std::string neg_infty = "\'-inf\'";
      return neg_infty;
    }
  }
} /* stringValue (int) */

inline const std::string stringValue(const long value,
    const char* format = "%ld") {
  if (!isInfinity(std::abs(value), std::numeric_limits<long>::max())) {
    char temp[500];
    snprintf(temp, sizeof(temp) / sizeof(char), format, value);
    std::string tmp(temp);
    return tmp;
  } else {
    if (greaterThanVal(value, 0.0)) {
      const std::string infty = "\'inf\'";
      return infty;
    } else {
      const std::string neg_infty = "\'-inf\'";
      return neg_infty;
    }
  }
} /* stringValue (long) */

inline const std::string stringValue(const double value,
    const char* format = "%f") {
  if (!isInfinity(std::abs(value))) {
    char temp[500];
    snprintf(temp, sizeof(temp) / sizeof(char), format, value);
    std::string tmp(temp);
    return tmp;
  } else {
    if (greaterThanVal(value, 0.0)) {
      const std::string infty = "\'+inf\'";
      return infty;
    } else {
      const std::string neg_infty = "\'-inf\'";
      return neg_infty;
    }
  }
} /* stringValue (double) */

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

double dotProduct(const CoinPackedVector& vec1, const double* vec2);
double dotProduct(const CoinPackedVector& vec1, const CoinPackedVector& vec2);

// Compute dot product using compensated summation to have small
// numerical error. First version: dense vectors
double dotProduct(const double* a, const double* b, int dimension);

// Second version: first vector is sparse, second one is dense
double dotProduct(int sizea, const int* indexa, const double* a,
    const double* b);

// Third version: sparse vectors
double dotProduct(int sizea, const int* indexa, const double* a, int sizeb,
    const int* indexb, const double* b);
