/** 
 * @file SolverHelper.hpp
 * @brief Various helpful functions for solver-agnostic queries
 * @author A. M. Kazachkov
 */
#pragma once

#include <limits> // numeric_limits

#ifdef USE_CLP
  #include <OsiClpSolverInterface.hpp>
#else
  #include <OsiSolverInterface.hpp>
#endif /* USE_CLP */

#ifdef USE_CBC
class CbcModel;
/// @brief Set log level and print frequency for \p cbc_model
void setIPSolverParameters(CbcModel* const cbc_model,
    const int verbosity =
#ifdef TRACE
        1
#else
        0
#endif
    );
#endif /* USE_CBC */

class OsiCuts;

/// @brief Set log level and max time
void setLPSolverParameters(OsiSolverInterface* const solver,
    const int verbosity =
#ifdef TRACE
        1,
#else
        0,
#endif
    const double max_time = std::numeric_limits<double>::max());

/// @brief Set solver time limit
void setTimeLimit(OsiSolverInterface* const solver, const double timeLimit);
/// @brief Check if time limit is hit
bool hitTimeLimit(const OsiSolverInterface* const solver);

#ifdef USE_CLP
/// @brief Set options for Clp to use strong branching
void setupClpForStrongBranching(OsiClpSolverInterface* const solver,
    const int hot_start_iter_limit = std::numeric_limits<int>::max());
#endif

/// @brief Sets message handler and special options when using solver as part of B&B
void setupClpForCbc(OsiSolverInterface* const solver,
    const int verbosity =
#ifdef TRACE
        1,
#else
        0,
#endif
    const double max_time = std::numeric_limits<double>::max(),
    const int hot_start_iter_limit = std::numeric_limits<int>::max());

/// @brief Set the cut solver objective coefficients based on a packed vector
void addToObjectiveFromPackedVector(OsiSolverInterface* const solver,
    const CoinPackedVectorBase* vec, const bool zeroOut = false,
    const double mult = 1., const std::vector<int>* const nonZeroColIndices = 0,
    const bool SHOULD_SCALE = true);

/// @brief Set all objective coefficients to \p val
void setConstantObjectiveFromPackedVector(OsiSolverInterface* const solver,
    const double val = 0., const int numIndices = 0, const int* indices = 0);

/// @brief Set solution (set to zero if \p sol is NULL)
void setSolverSolution(OsiSolverInterface* const solver, const double* const sol = NULL);

/// @brief Overload solve from hot start because of issues
bool solveFromHotStart(OsiSolverInterface* const solver, const int col,
    const bool isChangedUB, const double origBound, const double newBound);

/// @brief Checks whether a solver is optimal
bool checkSolverOptimality(OsiSolverInterface* const solver,
    const bool exitOnDualInfeas,
    const double timeLimit = std::numeric_limits<double>::max(),
    const int maxNumResolves = 2);

/// @brief Enable factorization, and check whether some cleanup needs to be done
bool enableFactorization(OsiSolverInterface* const solver, const double EPS, const int resolveFlag = 2);

/// @brief Generic way to apply a set of cuts to a solver
double applyCutsCustom(OsiSolverInterface* const solver, const OsiCuts& cs,
    FILE* logfile, const int startCutIndex, const int numCutsOverload, double compare_obj);
/// @brief Generic way to apply a set of cuts to a solver
double applyCutsCustom(OsiSolverInterface* const solver, const OsiCuts& cs,
    FILE* logfile = NULL,
    const int numCutsOverload = -1,
    double compare_obj = std::numeric_limits<double>::lowest());

/// @brief Get objective value from solver with cuts added
double getObjValue(OsiSolverInterface* const solver, const OsiCuts* const cuts);

#ifdef USE_CLP
///@{ 
/// @name Decide status of variables using ClpSimplex status

/**
 * @brief Status of structural variables.
 * @details The only tricky thing is fixed variables.
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

/// Check if \p stat is ClpSimplex::Status::isFree
inline bool isNonBasicFreeVar(const ClpSimplex::Status stat) {
  return (stat == ClpSimplex::Status::isFree);
}

/// This may be wrong for fixed variables that end up being counted as at ub
inline bool isNonBasicUBVar(const ClpSimplex::Status stat) {
  return (stat == ClpSimplex::Status::atUpperBound);
}

/// This may be wrong for fixed variables that end up being counted as at lb
inline bool isNonBasicLBVar(const ClpSimplex::Status stat) {
  return (stat == ClpSimplex::Status::atLowerBound);
}

/// Check if \p stat is ClpSimplex::Status::isFixed
inline bool isNonBasicFixedVar(const ClpSimplex::Status stat) {
  return (stat == ClpSimplex::Status::isFixed);
}
///@}
#endif // USE_CLP

// <!---------------- Aliases from cstat/rstat ---------------->
///@{
/// @name Aliases from cstat/rstat

/// @brief Check if column (structural variable) is basic
/// @details Reverts to isBasicVar(const int) if USE_CLP is defined
inline bool isBasicCol(const std::vector<int> &cstat, const int col) {
#ifdef USE_CLP
  return isBasicVar(static_cast<ClpSimplex::Status>(cstat[col]));
#else
  return cstat[col] == 1;
#endif
}

/// @brief Check if column (structural variable) is nonbasic at its upper bound
/// @details Reverts to isNonBasicUBVar(const int) if USE_CLP is defined
inline bool isNonBasicUBCol(const std::vector<int> &cstat, const int col) {
#ifdef USE_CLP
  return isNonBasicUBVar(static_cast<ClpSimplex::Status>(cstat[col]));
#else
  return cstat[col] == 2;
#endif
}

/// @brief Check if column (structural variable) is nonbasic at its lower bound
/// @details Reverts to isNonBasicLBVar(const int) if USE_CLP is defined
inline bool isNonBasicLBCol(const std::vector<int> &cstat, const int col) {
#ifdef USE_CLP
  return isNonBasicLBVar(static_cast<ClpSimplex::Status>(cstat[col]));
#else
  return cstat[col] == 3;
#endif
}

/// @brief Check if row (slack variable) is basic
/// @details Flipping we are worried about in isBasicSlack(const OsiSolverInterface* const, const int) has already occurred
/// Reverts to isBasicVar(const int) if USE_CLP is defined
inline bool isBasicSlack(const std::vector<int> &rstat, const int row) {
#ifdef USE_CLP
  return isBasicVar(static_cast<ClpSimplex::Status>(rstat[row]));
#else
  return rstat[row] == 1;
#endif
}

/// @brief Check if row (slack variable) is nonbasic at its upper bound
/// @details Reverts to isNonBasicUBVar(const int) if USE_CLP is defined
inline bool isNonBasicUBSlack(const std::vector<int> &rstat, const int row) {
#ifdef USE_CLP
  return isNonBasicUBVar(static_cast<ClpSimplex::Status>(rstat[row]));
#else
  return rstat[row] == 2;
#endif
}

/// @brief Check if row (slack variable) is nonbasic at its lower bound
/// @details Reverts to isNonBasicLBVar(const int) if USE_CLP is defined
inline bool isNonBasicLBSlack(const std::vector<int> &rstat, const int row) {
#ifdef USE_CLP
  return isNonBasicLBVar(static_cast<ClpSimplex::Status>(rstat[row]));
#else
  return rstat[row] == 3;
#endif
}
///@} // aliases using cstat/rstat

// <!---------------- Aliases from solver ---------------->
/** @name Check basic/nonbasic status of structural/slack variables given OsiSolverInterface pointer
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
///@{

/// @brief Check if structural column is basic
bool isBasicCol(const OsiSolverInterface* const solver, const int col);

/// @brief Check if column (structural variable) is nonbasic and free
bool isNonBasicFreeCol(const OsiSolverInterface* const solver, const int col);

/// @brief Check if column (structural variable) is nonbasic at its upper bound
bool isNonBasicUBCol(const OsiSolverInterface* const solver, const int col);

/// @brief Check if column (structural variable) is nonbasic at its lower bound
bool isNonBasicLBCol(const OsiSolverInterface* const solver, const int col);

/// @brief Check if column (structural variable) is fixed
bool isNonBasicFixedCol(const OsiSolverInterface* const solver, const int col);

/// @brief Check if slack variable is basic
bool isBasicSlack(const OsiSolverInterface* const solver, const int row);

/// @brief Check if slack variable is nonbasic and free
bool isNonBasicFreeSlack(const OsiSolverInterface* const solver, const int row);

/// @brief Check if slack variable is nonbasic and at its upper bound
bool isNonBasicUBSlack(const OsiSolverInterface* const solver, const int row);

/// @brief Check if slack variable is nonbasic and at its lower bound
bool isNonBasicLBSlack(const OsiSolverInterface* const solver, const int row);

/// @brief Check if slack variable is nonbasic and fixed
bool isNonBasicFixedSlack(const OsiSolverInterface* const solver, const int row);
///@} // check basic / nb status given OsiSolverInterface

///@{
/// @name Aliases checking if var is structural or slack
/// @brief Check if variable is nonbasic at its lower bound
inline bool isBasicVar(const OsiSolverInterface* const solver, const int var) {
#ifdef USE_CLP
  try {
    return isBasicVar(dynamic_cast<const OsiClpSolverInterface* const>(solver)->getModelPtr()->getStatus(var));
  } catch (std::exception& e) {

  }
#endif
  return
      (var < solver->getNumCols()) ?
          isBasicCol(solver, var) :
          isBasicSlack(solver, var - solver->getNumCols());
}

/** 
 * @brief Check if variable is nonbasic at its lower bound
 *
 * @details Implementation is not simply using getStatus because of fixed variables 
 */
inline bool isNonBasicLBVar(const OsiSolverInterface* const solver, const int var) {
  return
      (var < solver->getNumCols()) ?
          isNonBasicLBCol(solver, var) :
          isNonBasicLBSlack(solver, var - solver->getNumCols());
}

/** 
 * @brief Check if variable is nonbasic at its upper bound
 *
 * @details Implementation is not simply using getStatus because of fixed variables 
 */
inline bool isNonBasicUBVar(const OsiSolverInterface* const solver, const int var) {
  return
      (var < solver->getNumCols()) ?
          isNonBasicUBCol(solver, var) :
          isNonBasicUBSlack(solver, var - solver->getNumCols());
}

/** 
 * @brief Check if variable is nonbasic and free
 *
 * @details Implementation is not simply using getStatus because of fixed variables 
 */
inline bool isNonBasicFreeVar(const OsiSolverInterface* const solver, const int var) {
  return
      (var < solver->getNumCols()) ?
          isNonBasicFreeCol(solver, var) :
          isNonBasicFreeSlack(solver, var - solver->getNumCols());
}

/**
 * @brief Check if a variable is nonbasic and fixed
 */
inline bool isNonBasicFixedVar(const OsiSolverInterface* const solver, const int var) {
  return
      (var < solver->getNumCols()) ?
          isNonBasicFixedCol(solver, var) :
          isNonBasicFixedSlack(solver, var - solver->getNumCols());
}

/**
 * @brief Check if a variable corresponds to a ray of the basis cone (equivalently, it is not a basic variable)
 *
 * @details The variables that lead to rays are nonbasic structural variables
 * and nonbasic slacks not coming from equality constraints.
 */
inline bool isRayVar(const OsiSolverInterface* const solver, const int var) {
  return !isBasicVar(solver, var);
}
///@} // status of variables (checking if it is structural or slack)
