#pragma once
#include "utility.hpp"

#include <limits> // numeric_limits

#ifdef USE_CLP
  #include <OsiClpSolverInterface.hpp>
#else
  #include <OsiSolverInterface.hpp>
#endif

#ifdef USE_CBC
#include <CbcModel.hpp>
void setIPSolverParameters(CbcModel* const cbc_model);
#endif

void setLPSolverParameters(OsiSolverInterface* const solver,
    const double max_time = std::numeric_limits<double>::max());

/**
 * Checks whether a solver is optimal
 * Something that can go wrong (e.g., bc1 -64 sb5 tmp_ind = 14):
 *  Solver is declared optimal for the scaled problem but there are primal or dual infeasibilities in the unscaled problem
 *  In this case, secondary status will be 2 or 3
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

/*****************************
 * Decide status of variables
 *****************************/
#ifdef USE_CLP
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
inline bool isNonBasicUBVar(const ClpSimplex::Status stat) {
  return (stat == ClpSimplex::Status::atUpperBound);
}

// This may be wrong for fixed variables that end up being counted as at lb
inline bool isNonBasicLBVar(const ClpSimplex::Status stat) {
  return (stat == ClpSimplex::Status::atLowerBound);
}

inline bool isNonBasicFixedVar(const ClpSimplex::Status stat) {
  return (stat == ClpSimplex::Status::isFixed);
}
#endif // USE_CLP

// Aliases from cstat/rstat
// Flipping we are worried about below has already occurred
inline bool isBasicCol(const std::vector<int> &cstat, const int col) {
#ifdef USE_CLP
  return isBasicVar(static_cast<ClpSimplex::Status>(cstat[col]));
#else
  return cstat[col] == 1;
#endif
}

inline bool isNonBasicUBCol(const std::vector<int> &cstat, const int col) {
#ifdef USE_CLP
  return isNonBasicUBVar(static_cast<ClpSimplex::Status>(cstat[col]));
#else
  return cstat[col] == 2;
#endif
}

inline bool isNonBasicLBCol(const std::vector<int> &cstat, const int col) {
#ifdef USE_CLP
  return isNonBasicLBVar(static_cast<ClpSimplex::Status>(cstat[col]));
#else
  return cstat[col] == 3;
#endif
}

inline bool isBasicSlack(const std::vector<int> &rstat, const int row) {
#ifdef USE_CLP
  return isBasicVar(static_cast<ClpSimplex::Status>(rstat[row]));
#else
  return rstat[row] == 1;
#endif
}

inline bool isNonBasicUBSlack(const std::vector<int> &rstat, const int row) {
#ifdef USE_CLP
  return isNonBasicUBVar(static_cast<ClpSimplex::Status>(rstat[row]));
#else
  return rstat[row] == 2;
#endif
}

inline bool isNonBasicLBSlack(const std::vector<int> &rstat, const int row) {
#ifdef USE_CLP
  return isNonBasicLBVar(static_cast<ClpSimplex::Status>(rstat[row]));
#else
  return rstat[row] == 3;
#endif
}

bool isBasicCol(const OsiSolverInterface* const solver, const int col);

bool isNonBasicFreeCol(const OsiSolverInterface* const solver, const int col);

bool isNonBasicUBCol(const OsiSolverInterface* const solver, const int col);

bool isNonBasicLBCol(const OsiSolverInterface* const solver, const int col);

bool isNonBasicFixedCol(const OsiSolverInterface* const solver, const int col);

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
bool isBasicSlack(const OsiSolverInterface* const solver, const int row);

bool isNonBasicFreeSlack(const OsiSolverInterface* const solver, const int row);

// NOTE: We do at *lower* bound on purpose because of comment above:
// the underlying model flips how it holds these with respect to Clp
bool isNonBasicUBSlack(const OsiSolverInterface* const solver, const int row);

// NOTE: We do at *upper* bound on purpose because of comment above:
// the underlying model flips how it holds these with respect to Clp
bool isNonBasicLBSlack(const OsiSolverInterface* const solver, const int row);

bool isNonBasicFixedSlack(const OsiSolverInterface* const solver, const int row);

// Aliases
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

/** Not simply calling using getStatus because of fixed variables */
inline bool isNonBasicLBVar(const OsiSolverInterface* const solver, const int var) {
  return
      (var < solver->getNumCols()) ?
          isNonBasicLBCol(solver, var) :
          isNonBasicLBSlack(solver, var - solver->getNumCols());
}

/** Not simply calling using getStatus because of fixed variables */
inline bool isNonBasicUBVar(const OsiSolverInterface* const solver, const int var) {
  return
      (var < solver->getNumCols()) ?
          isNonBasicUBCol(solver, var) :
          isNonBasicUBSlack(solver, var - solver->getNumCols());
}

/** Not simply calling using getStatus because of fixed variables */
inline bool isNonBasicFreeVar(const OsiSolverInterface* const solver, const int var) {
  return
      (var < solver->getNumCols()) ?
          isNonBasicFreeCol(solver, var) :
          isNonBasicFreeSlack(solver, var - solver->getNumCols());
}

/**
 * The variables that lead to rays are nonbasic structural variables
 * and nonbasic slacks not coming from equality constraints.
 */
inline bool isRayVar(const OsiSolverInterface* const solver, const int var) {
  return !isBasicVar(solver, var);
}
