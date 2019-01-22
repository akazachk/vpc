#pragma once

#include <cstdio>
#include <iostream> // cerr
#include <string>
#include <cmath> // abs
#include <limits> // numeric_limits

#include <CbcModel.hpp>
#include <OsiClpSolverInterface.hpp>

#include "params.hpp"

//#define error_msg(str, fmt) \
//  char str[500]; \
//  snprintf(str, sizeof(str) / sizeof(char), "*** ERROR: %s:%d: " fmt, __FILE__, __LINE__); \
//  std::cerr << str
//
//#define warning_msg(str, fmt) \
//  char str[500]; \
//  snprintf(str, sizeof(str) / sizeof(char), "*** WARNING: %s:%d: " fmt, __FILE__, __LINE__); \
//  std::cerr << str

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

/**
 * Checks whether a solver is optimal
 * Something that can go wrong (e.g., bc1 -64 sb5 tmp_ind = 14):
 * 	Solver is declared optimal for the scaled problem but there are primal or dual infeasibilities in the unscaled problem
 * 	In this case, secondary status will be 2 or 3
 * Sometimes cleanup helps
 * Sometimes things break after enabling factorization (e.g., secondary status = 0, but sum infeasibilities > 0)
 */
bool checkSolverOptimality(OsiSolverInterface* const solver,
    const bool exitOnDualInfeas,
    const double timeLimit = std::numeric_limits<double>::max(),
    const int maxNumResolves = 2);

/**
 * Enable factorization, and check whether some cleanup needs to be done
 */
bool enableFactorization(OsiSolverInterface* const solver, const double EPS, const int resolveFlag = 2);

void setLPSolverParameters(OsiSolverInterface* const solver,
    const double max_time = std::numeric_limits<double>::max());
void setIPSolverParameters(CbcModel* const cbc_model);
double getObjValueFromFile(const VPCParameters& params);

void setCompNBCoor(CoinPackedVector& vec, double& objViolation,
    const VPCParameters& params, const OsiSolverInterface* const tmpSolver,
    const OsiSolverInterface* const origSolver,
    const std::vector<int>& nonBasicVarIndex,
    const std::vector<double>& nonBasicReducedCost, const int deletedVar = -1);
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

/*****************************
 * Decide status of variables
 *****************************/
/**
 * Status of structural variables.
 * The only tricky thing is fixed variables.
 * This is decided by reduced cost, as for slacks.
 * This is a minimization problem.
 * The reduced cost, at optimality, is
 *   >= 0 for a variable at its lower bound,
 *     (increasing the variable would increase the objective value)
 *   <= 0 for a variable at its upper bound
 *     (decreasing the variable would increase the objective value).
 * Thus, if a reduced cost is < 0, the variable will be treated as being at its ub.
 */
inline bool isBasicVar(const ClpSimplex::Status stat) {
  return (stat == ClpSimplex::Status::basic);
}

inline bool isNonBasicFreeVar(const ClpSimplex::Status stat) {
  return (stat == ClpSimplex::Status::isFree);
}

// This may be wrong for fixed variables that end up being counted as at ub
inline bool isNonBasicUBCol(const ClpSimplex::Status stat) {
  return (stat == ClpSimplex::Status::atUpperBound);
}

// This may be wrong for fixed variables that end up being counted as at lb
inline bool isNonBasicLBCol(const ClpSimplex::Status stat) {
  return (stat == ClpSimplex::Status::atLowerBound);
}

inline bool isNonBasicFixedVar(const ClpSimplex::Status stat) {
  return (stat == ClpSimplex::Status::isFixed);
}

// Aliases from cstat/rstat
// Flipping we are worried about below has already occurred
inline bool isBasicCol(const std::vector<int> &cstat, const int col) {
  return isBasicVar(static_cast<ClpSimplex::Status>(cstat[col]));
}

inline bool isNonBasicUBCol(const std::vector<int> &cstat, const int col) {
  return isNonBasicUBCol(static_cast<ClpSimplex::Status>(cstat[col]));
}

inline bool isNonBasicLBCol(const std::vector<int> &cstat, const int col) {
  return isNonBasicLBCol(static_cast<ClpSimplex::Status>(cstat[col]));
}

inline bool isBasicSlack(const std::vector<int> &rstat, const int row) {
  return isBasicVar(static_cast<ClpSimplex::Status>(rstat[row]));
}

inline bool isNonBasicUBSlack(const std::vector<int> &rstat, const int row) {
  return isNonBasicUBCol(static_cast<ClpSimplex::Status>(rstat[row]));
}

inline bool isNonBasicLBSlack(const std::vector<int> &rstat, const int row) {
  return isNonBasicLBCol(static_cast<ClpSimplex::Status>(rstat[row]));
}

// Aliases from solver
inline bool isBasicCol(const OsiClpSolverInterface* const solver, const int col) {
  return isBasicVar(solver->getModelPtr()->getColumnStatus(col));
}

inline bool isNonBasicFreeCol(const OsiClpSolverInterface* const solver,
    const int col) {
  return isNonBasicFreeVar(solver->getModelPtr()->getColumnStatus(col));
}

inline bool isNonBasicUBCol(const OsiClpSolverInterface* const solver,
    const int col) {
  const ClpSimplex::Status status = solver->getModelPtr()->getColumnStatus(col);
  return ((status == ClpSimplex::Status::atUpperBound)
      || (status == ClpSimplex::Status::isFixed
          && lessThanVal(solver->getReducedCost()[col], 0.)));
}

inline bool isNonBasicLBCol(const OsiClpSolverInterface* const solver,
    const int col) {
  const ClpSimplex::Status status = solver->getModelPtr()->getColumnStatus(col);
  return ((status == ClpSimplex::Status::atLowerBound)
        || (status == ClpSimplex::Status::isFixed
            && !lessThanVal(solver->getReducedCost()[col], 0.)));
}

inline bool isNonBasicFixedCol(const OsiClpSolverInterface* const solver,
    const int col) {
  return isNonBasicFixedVar(solver->getModelPtr()->getColumnStatus(col));
}

/**
 * Status of slack variables.
 * The only tricky aspect is that the model pointer row status for slack variables
 * at ub is actually that they are at lb, and vice versa.
 * Basically, when Clp says a row is atLowerBound,
 * then it means it is a >= row, so that rowActivity == rowLower.
 * In Osi, the status of the corresponding slack is at its UPPER bound,
 * because Osi keeps all slacks with a +1 coefficient in the row.
 * When it says a row is atUpperBound,
 * then the row is a <= row and rowActivity == rowUpper.
 * In Osi, the status of the corresponding slack is at its LOWER bound,
 * since this is the regular slack we know and love.
 *
 * One more issue: fixed slack (artificial) variables.
 * The status in Osi is decided by the reduced cost.
 * This is a minimization problem; the solution is optimal when the reduced cost is
 * (as above)
 *   >= 0 for variables at their lower bound,
 *   <= 0 for variables at their upper bound.
 * Thus, if a reduced cost is < 0, the variable will be treated as being at its ub.
 * How do we get the reduced cost for slacks?
 * It is actually just the negative of the dual variable value,
 * which is stored in rowPrice.
 */
inline bool isNonBasicFreeSlack(const OsiClpSolverInterface* const solver,
    const int row) {
  return isNonBasicFreeVar(solver->getModelPtr()->getRowStatus(row));
}

inline bool isBasicSlack(const OsiClpSolverInterface* const solver, const int row) {
  return isBasicVar(solver->getModelPtr()->getRowStatus(row));
}

// NOTE: We do at *lower* bound on purpose because of comment above:
// the underlying model flips how it holds these with respect to Clp
inline bool isNonBasicUBSlack(const OsiClpSolverInterface* const solver,
    const int row) {
  const ClpSimplex::Status status = solver->getModelPtr()->getRowStatus(row);
  return ((status == ClpSimplex::Status::atLowerBound)
      || (status == ClpSimplex::Status::isFixed
          && greaterThanVal(solver->getRowPrice()[row], 0.)));
}

// NOTE: We do at *upper* bound on purpose because of comment above:
// the underlying model flips how it holds these with respect to Clp
inline bool isNonBasicLBSlack(const OsiClpSolverInterface* const solver,
    const int row) {
  const ClpSimplex::Status status = solver->getModelPtr()->getRowStatus(row);
  return ((status == ClpSimplex::Status::atUpperBound)
      || (status == ClpSimplex::Status::isFixed
          && !greaterThanVal(solver->getRowPrice()[row], 0.)));
}

inline bool isNonBasicFixedSlack(const OsiClpSolverInterface* const solver,
    const int row) {
  return isNonBasicFixedVar(solver->getModelPtr()->getRowStatus(row));
}

inline bool isBasicVar(const OsiClpSolverInterface* const solver, const int col) {
  return isBasicVar(solver->getModelPtr()->getStatus(col));
}

/**
 * The variables that lead to rays are non-basic structural variables
 * and non-basic slacks not coming from equality constraints.
 */
bool isRayVar(const OsiClpSolverInterface* const solver, const int var);

inline bool isNonBasicLBVar(const OsiClpSolverInterface* const solver,
    const int var) {
  return
      (var < solver->getNumCols()) ?
          isNonBasicLBCol(solver, var) :
          isNonBasicLBSlack(solver, var - solver->getNumCols());
}

inline bool isNonBasicUBVar(const OsiClpSolverInterface* const solver,
    const int var) {
  return
      (var < solver->getNumCols()) ?
          isNonBasicUBCol(solver, var) :
          isNonBasicUBSlack(solver, var - solver->getNumCols());
}

inline bool isNonBasicFreeVar(const OsiClpSolverInterface* const solver,
    const int var) {
  return
      (var < solver->getNumCols()) ?
          isNonBasicFreeCol(solver, var) :
          isNonBasicFreeSlack(solver, var - solver->getNumCols());
}

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
