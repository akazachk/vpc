#include "SolverHelper.hpp"

#ifdef USE_CBC
void setIPSolverParameters(CbcModel* const cbc_model) {
#ifdef TRACE
  cbc_model->setLogLevel(3);
  cbc_model->messagesPointer()->setDetailMessages(10, 10000, (int *) NULL);
#else
  cbc_model->setLogLevel(0);
  cbc_model->messagesPointer()->setDetailMessages(10,5,5000);
#endif
  if (cbc_model->solver()) {
    cbc_model->solver()->setHintParam(OsiDoReducePrint, true, OsiHintTry);
  }
  cbc_model->setPrintFrequency(1);
} /* setIPSolverParameters (Cbc) */
#endif

void setLPSolverParameters(OsiSolverInterface* const solver,
    const double max_time) {
#ifndef TRACE
  solver->messageHandler()->setLogLevel(0);
#ifdef USE_CLP
  try {
    dynamic_cast<OsiClpSolverInterface*>(solver)->getModelPtr()->messageHandler()->setLogLevel(0);
  } catch (std::exception& e) {

  }
#endif
#endif
#ifdef USE_CLP
  try {
    dynamic_cast<OsiClpSolverInterface*>(solver)->getModelPtr()->setMaximumSeconds(max_time);
  } catch (std::exception& e) {

  }
#endif
  // Try turning on scaling with enableFactorization
//  solver->setSpecialOptions(solver->specialOptions() | 512);
} /* setLPSolverParameters (OsiClp) */

/**
 * @brief Checks whether a solver is optimal
 *
 * Something that can go wrong (e.g., bc1 -64 sb5 tmp_ind = 14):
 *  Solver is declared optimal for the scaled problem but there are primal or dual infeasibilities in the unscaled problem
 *  In this case, secondary status will be 2 or 3
 * Sometimes cleanup helps
 * Sometimes things break after enabling factorization (e.g., secondary status = 0, but sum infeasibilities > 0)
 */
bool checkSolverOptimality(OsiSolverInterface* const solver,
    const bool exitOnDualInfeas, const double timeLimit,
    const int maxNumResolves) {
#ifdef USE_CLP
  OsiClpSolverInterface* clpsolver = NULL;
  try {
    clpsolver = dynamic_cast<OsiClpSolverInterface*>(solver);
  } catch (std::exception& e) {
    // Disregard / catch below
  }
  if (!clpsolver) {
    return solver->isProvenOptimal();
  }

  // If it is supposedly proven primal infeasible, might be good to do a resolve first
  // This is, e.g., something that arises with miplib2003/pp08aCUTS_presolved -32
  int resolve_count = 0;
  if (clpsolver->isProvenPrimalInfeasible()) {
    clpsolver->getModelPtr()->setNumberIterations(0);
    clpsolver->getModelPtr()->setMaximumSeconds(timeLimit);
    clpsolver->resolve();
    resolve_count++;
  }

  // First clean
  if (resolve_count < maxNumResolves && clpsolver) {
    int status = clpsolver->getModelPtr()->secondaryStatus();
    bool is_cleaned = false;
    if (status == 2 || !isZero(clpsolver->getModelPtr()->sumPrimalInfeasibilities())) {
      clpsolver->getModelPtr()->cleanup(1);
      is_cleaned = true;
    } else if (status == 3 || status == 9 || !isZero(clpsolver->getModelPtr()->sumDualInfeasibilities())) {
      clpsolver->getModelPtr()->cleanup(2);
      is_cleaned = true;
    } else if (status == 4) {
      clpsolver->getModelPtr()->cleanup(1);
      is_cleaned = true;
    }
    if (is_cleaned) {
      resolve_count++;
      return checkSolverOptimality(solver, exitOnDualInfeas, timeLimit, 0);
    }

    // Do some resolves, but first save whether the initial status is dual infeasible
    // The reason is that for neos15 w/str=2, we were getting a dual infeasible problem
    // turn into a primal infeasible problem after resolves
    // (probably the strengthening causes numerical issues)
    const bool oldStatusDualInfeasible = (clpsolver->isProvenDualInfeasible());
    const bool oldStatusOptimal = (clpsolver->isProvenOptimal());
    bool resolve = true;
    double infeas = std::numeric_limits<double>::max();
    while (resolve && resolve_count < maxNumResolves) {
      const double curr_infeas =
          clpsolver->getModelPtr()->sumPrimalInfeasibilities()
              + clpsolver->getModelPtr()->sumDualInfeasibilities();
      resolve = (resolve_count == 0)
          || (!isZero(curr_infeas) && (curr_infeas < infeas));
      if (resolve) {
        clpsolver->getModelPtr()->setNumberIterations(0);
        clpsolver->getModelPtr()->setMaximumSeconds(timeLimit);
        clpsolver->resolve();
        resolve_count++;
      }
      infeas = curr_infeas;
    }

    // If dual infeas -> primal infeas, do a dual resolve, which should fix the issue
    if (oldStatusDualInfeasible && clpsolver->isProvenPrimalInfeasible()) {
      clpsolver->getModelPtr()->dual(2,0); // just do values pass
      if (!clpsolver->isProvenDualInfeasible() && !clpsolver->isProvenOptimal()) {
        clpsolver->getModelPtr()->setProblemStatus(2);
        resolve_count = maxNumResolves; // stop resolving
      }
    }

    if (oldStatusOptimal && clpsolver->isProvenPrimalInfeasible()) {
      clpsolver->enableFactorization();
      clpsolver->getModelPtr()->setNumberIterations(0);
      clpsolver->getModelPtr()->setMaximumSeconds(timeLimit);
      clpsolver->resolve();
      resolve_count++;
      clpsolver->disableFactorization();
    }

    // Clean once more if needed and possible
    status = clpsolver->getModelPtr()->secondaryStatus();
    is_cleaned = false;
    if (status == 2
        || !isZero(clpsolver->getModelPtr()->sumPrimalInfeasibilities())) {
      clpsolver->getModelPtr()->cleanup(1);
      is_cleaned = true;
    } else if (status == 3 || status == 9
        || !isZero(clpsolver->getModelPtr()->sumDualInfeasibilities())) {
      clpsolver->getModelPtr()->cleanup(2);
      is_cleaned = true;
    } else if (status == 4) {
      clpsolver->getModelPtr()->cleanup(1);
      is_cleaned = true;
    }
    if (is_cleaned) {
      resolve_count++;
      return checkSolverOptimality(solver, exitOnDualInfeas, timeLimit, 0);
    }
  }
#endif // USE_CLP

  if (solver->isProvenPrimalInfeasible()) {
    return false;
  } else if (!(solver->isProvenOptimal())) {
    // Sometimes need to resolve once more to get the correct status
    if (resolve_count < maxNumResolves) {
      solver->resolve();
    }
    if (solver->isProvenPrimalInfeasible()) {
      return false;
    } else if (solver->isProvenDualInfeasible()) {
      if (exitOnDualInfeas) {
        error_msg(errstr,
            "Solver is dual infeasible. Check why this happened!\n");
        exit(1);
      } else {
        return false;
      }
    } else if (!(solver->isProvenOptimal())) {
      return false;
    }
  }

  return true;
} /* checkSolverOptimality */

/**
 * @brief Enable factorization, and check whether some cleanup needs to be done
 */
bool enableFactorization(OsiSolverInterface* const solver, const double EPS, const int resolveFlag) {
  solver->enableFactorization();

#ifdef USE_CLP
  if (resolveFlag > 0) {
    try {
      OsiClpSolverInterface* clpsolver =
          dynamic_cast<OsiClpSolverInterface*>(solver);

      // After enabling factorization, things sometimes look different and it is worth resolving (e.g., bc1 with num_strong = 5, row price on row 458)
      // This is motivated by seymour-disj-10; after adding one round of GMICs, row 5077 had negative row price
      if ((clpsolver->getModelPtr()->sumPrimalInfeasibilities() > EPS)
          || (clpsolver->getModelPtr()->sumDualInfeasibilities() > EPS)) {
        if (resolveFlag == 1) {
          clpsolver->initialSolve(); // hopefully this fixes things rather than breaks things; actually it breaks things, e.g., for rd2 SICs for coral/neos17, the variable basic in row 648 changes from 164 to 23
        } else {
          clpsolver->resolve();
        }
        clpsolver->disableFactorization();
        clpsolver->enableFactorization(); // this will hopefully solve the neos17 problem
      }
    } catch (std::exception& e) {
      // Disregard
    }
  }
#endif // USE_CLP

  return solver->isProvenOptimal();
} /* enableFactorization */

// Aliases from solver
bool isBasicCol(const OsiSolverInterface* const solver, const int col) {
#ifdef USE_CLP
  try {
    const OsiClpSolverInterface* const clpsolver = dynamic_cast<const OsiClpSolverInterface* const>(solver);
    return isBasicVar(clpsolver->getModelPtr()->getColumnStatus(col));
  } catch (std::exception& e) {
  }
#endif
  // If not using clp, we have to get basis status
  std::vector<int> cstat(solver->getNumCols()), rstat(solver->getNumRows());
  solver->getBasisStatus(&cstat[0], &rstat[0]);
  return isBasicCol(cstat, col);
}

bool isNonBasicFreeCol(const OsiSolverInterface* const solver, const int col) {
#ifdef USE_CLP
  try {
    const OsiClpSolverInterface* const clpsolver = dynamic_cast<const OsiClpSolverInterface* const>(solver);
    return isNonBasicFreeVar(clpsolver->getModelPtr()->getColumnStatus(col));
  } catch (std::exception& e) {
  }
#endif
  // If not using clp, we have to get basis status
  std::vector<int> cstat(solver->getNumCols()), rstat(solver->getNumRows());
  solver->getBasisStatus(&cstat[0], &rstat[0]);
  return cstat[col] == 0;
}

bool isNonBasicUBCol(const OsiSolverInterface* const solver, const int col) {
#ifdef USE_CLP
  try {
    const OsiClpSolverInterface* const clpsolver = dynamic_cast<const OsiClpSolverInterface* const>(solver);
    const ClpSimplex::Status status = clpsolver->getModelPtr()->getColumnStatus(col);
    return ((status == ClpSimplex::Status::atUpperBound)
        || (status == ClpSimplex::Status::isFixed
            && lessThanVal(clpsolver->getReducedCost()[col], 0.)));
  } catch (std::exception& e) {
  }
#endif
  // If not using clp, we have to get basis status
  std::vector<int> cstat(solver->getNumCols()), rstat(solver->getNumRows());
  solver->getBasisStatus(&cstat[0], &rstat[0]);
  return isNonBasicUBCol(cstat, col);
}

bool isNonBasicLBCol(const OsiSolverInterface* const solver, const int col) {
#ifdef USE_CLP
  try {
    const OsiClpSolverInterface* const clpsolver = dynamic_cast<const OsiClpSolverInterface* const>(solver);
    const ClpSimplex::Status status = clpsolver->getModelPtr()->getColumnStatus(col);
    return ((status == ClpSimplex::Status::atLowerBound)
        || (status == ClpSimplex::Status::isFixed
            && !lessThanVal(clpsolver->getReducedCost()[col], 0.)));
  } catch (std::exception& e) {
  }
#endif
  // If not using clp, we have to get basis status
  std::vector<int> cstat(solver->getNumCols()), rstat(solver->getNumRows());
  solver->getBasisStatus(&cstat[0], &rstat[0]);
  return isNonBasicLBCol(cstat, col);
}

bool isNonBasicFixedCol(const OsiSolverInterface* const solver, const int col) {
#ifdef USE_CLP
  try {
    const OsiClpSolverInterface* const clpsolver = dynamic_cast<const OsiClpSolverInterface* const>(solver);
    return isNonBasicFixedVar(clpsolver->getModelPtr()->getColumnStatus(col));
  } catch (std::exception& e) {
  }
#endif
  // If not using clp, we have to get basis status
  std::vector<int> cstat(solver->getNumCols()), rstat(solver->getNumRows());
  solver->getBasisStatus(&cstat[0], &rstat[0]);
  return cstat[col] == 5;
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
bool isBasicSlack(const OsiSolverInterface* const solver, const int row) {
#ifdef USE_CLP
  try {
    const OsiClpSolverInterface* const clpsolver = dynamic_cast<const OsiClpSolverInterface* const>(solver);
    return isBasicVar(clpsolver->getModelPtr()->getRowStatus(row));
  } catch (std::exception& e) {
  }
#endif
  // If not using clp, we have to get basis status
  std::vector<int> cstat(solver->getNumCols()), rstat(solver->getNumRows());
  solver->getBasisStatus(&cstat[0], &rstat[0]);
  return isBasicSlack(rstat, row);
}

bool isNonBasicFreeSlack(const OsiSolverInterface* const solver, const int row) {
#ifdef USE_CLP
  try {
    const OsiClpSolverInterface* const clpsolver = dynamic_cast<const OsiClpSolverInterface* const>(solver);
    return isNonBasicFreeVar(clpsolver->getModelPtr()->getRowStatus(row));
  } catch (std::exception& e) {
  }
#endif
  // If not using clp, we have to get basis status
  std::vector<int> cstat(solver->getNumCols()), rstat(solver->getNumRows());
  solver->getBasisStatus(&cstat[0], &rstat[0]);
  return rstat[row] == 0;
}

// NOTE: We do at *lower* bound on purpose because of comment above:
// the underlying model flips how it holds these with respect to Clp
bool isNonBasicUBSlack(const OsiSolverInterface* const solver, const int row) {
#ifdef USE_CLP
  try {
    const OsiClpSolverInterface* const clpsolver = dynamic_cast<const OsiClpSolverInterface* const>(solver);
    const ClpSimplex::Status status = clpsolver->getModelPtr()->getRowStatus(row);
    return ((status == ClpSimplex::Status::atLowerBound)
        || (status == ClpSimplex::Status::isFixed
            && greaterThanVal(clpsolver->getRowPrice()[row], 0.)));
  } catch (std::exception& e) {
  }
#endif
  // If not using clp, we have to get basis status
  std::vector<int> cstat(solver->getNumCols()), rstat(solver->getNumRows());
  solver->getBasisStatus(&cstat[0], &rstat[0]);
  return isNonBasicUBSlack(rstat, row); // need to check if flipping is done here too
}

// NOTE: We do at *upper* bound on purpose because of comment above:
// the underlying model flips how it holds these with respect to Clp
bool isNonBasicLBSlack(const OsiSolverInterface* const solver, const int row) {
#ifdef USE_CLP
  try {
    const OsiClpSolverInterface* const clpsolver = dynamic_cast<const OsiClpSolverInterface* const>(solver);
    const ClpSimplex::Status status = clpsolver->getModelPtr()->getRowStatus(row);
    return ((status == ClpSimplex::Status::atUpperBound)
        || (status == ClpSimplex::Status::isFixed
            && !greaterThanVal(clpsolver->getRowPrice()[row], 0.)));
  } catch (std::exception& e) {
  }
#endif
  // If not using clp, we have to get basis status
  std::vector<int> cstat(solver->getNumCols()), rstat(solver->getNumRows());
  solver->getBasisStatus(&cstat[0], &rstat[0]);
  return isNonBasicLBSlack(rstat, row); // need to check if flipping is done here too
}

bool isNonBasicFixedSlack(const OsiSolverInterface* const solver, const int row) {
#ifdef USE_CLP
  try {
    const OsiClpSolverInterface* const clpsolver = dynamic_cast<const OsiClpSolverInterface* const>(solver);
    return isNonBasicFixedVar(clpsolver->getModelPtr()->getRowStatus(row));
  } catch (std::exception& e) {
  }
#endif
  // If not using clp, we have to get basis status
  std::vector<int> cstat(solver->getNumCols()), rstat(solver->getNumRows());
  solver->getBasisStatus(&cstat[0], &rstat[0]);
  return rstat[row] == 5;
}
