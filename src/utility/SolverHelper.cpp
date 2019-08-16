#include "SolverHelper.hpp"

#include "utility.hpp"

#include "OsiCuts.hpp"

#ifdef USE_CBC
#include <CbcModel.hpp>

void setIPSolverParameters(CbcModel* const cbc_model, const int verbosity) {
  if (verbosity > 0) {
    cbc_model->setLogLevel(3);
    cbc_model->messagesPointer()->setDetailMessages(10, 10000, (int *) NULL);
  } else {
    cbc_model->setLogLevel(0);
    cbc_model->messagesPointer()->setDetailMessages(10,5,5000);
  }
  if (cbc_model->solver()) { // catches when we are doing partial b&b
    cbc_model->solver()->setHintParam(OsiDoReducePrint, true, OsiHintTry);
  }
  cbc_model->setPrintFrequency(1);
} /* setIPSolverParameters (Cbc) */
#endif /* USE_CBC */

void setLPSolverParameters(OsiSolverInterface* const solver,
    const int verbosity,
    const double max_time) {
  if (verbosity == 0) {
    solver->messageHandler()->setLogLevel(0);
#ifdef USE_CLP
    try {
      dynamic_cast<OsiClpSolverInterface*>(solver)->getModelPtr()->messageHandler()->setLogLevel(0);
    } catch (std::exception& e) {
    }
#endif
  }
  setTimeLimit(solver, max_time);
  // Try turning on scaling with enableFactorization
//  solver->setSpecialOptions(solver->specialOptions() | 512);
} /* setLPSolverParameters */

void setTimeLimit(OsiSolverInterface* const solver, const double timeLimit) {
#ifdef USE_CLP
  try {
    dynamic_cast<OsiClpSolverInterface*>(solver)->getModelPtr()->setMaximumSeconds(timeLimit);
  } catch (std::exception& e) {

  }
#endif
} /* setTimeLimit */

bool hitTimeLimit(const OsiSolverInterface* const solver) {
#ifdef USE_CLP
  try {
    if (dynamic_cast<const OsiClpSolverInterface* const>(solver)->getModelPtr()->hitMaximumIterations())
      return true;
  } catch (std::exception& e) {

  }
#endif
  return false;
} /* hitTimeLimit */

#ifdef USE_CLP
/**
 * We need to be careful with the strong branching options;
 * sometimes the Clp strong branching fails, such as with arki001, branching down on variable 924
 */
void setupClpForStrongBranching(OsiClpSolverInterface* const solver, const int hot_start_iter_limit) {
  solver->setIntParam(OsiMaxNumIterationHotStart, hot_start_iter_limit);
  solver->setSpecialOptions(16); // use standard strong branching rather than clp's
} /* setupClpForStrongBranching */

/**
 * Sets message handler and special options when using solver as part of B&B
 * (in which we want to run full strong branching and enable the fixing of variables)
 */
void setupClpForCbc(OsiSolverInterface* const solver,
    const int hot_start_iter_limit) {
  setLPSolverParameters(solver);
  solver->setHintParam(OsiDoPresolveInInitial, false);
  solver->setHintParam(OsiDoPresolveInResolve, false);
  try {
    setupClpForStrongBranching(dynamic_cast<OsiClpSolverInterface*>(solver), hot_start_iter_limit);
    // Do not switch from dual to primal, or something to this effect;
    // This allows infeasible branches to be fixed during strong branching
    dynamic_cast<OsiClpSolverInterface*>(solver)->getModelPtr()->setMoreSpecialOptions(
        dynamic_cast<OsiClpSolverInterface*>(solver)->getModelPtr()->moreSpecialOptions() + 256);
  } catch (std::exception& e) {
    std::cerr << "Unable to cast solver as OsiClpSolverInterface." << std::endl;
    exit(1);
  }
} /* setupClpForCbc */
#endif // USE_CLP

/** Overload solve from hot start because of issues */
bool solveFromHotStart(OsiSolverInterface* const solver, const int col,
    const bool isChangedUB, const double origBound, const double newBound) {
  solver->solveFromHotStart();
  if (solver->isIterationLimitReached()) {
    // This sometimes happens, e.g., with arki001
    solver->unmarkHotStart();
    solver->resolve();
    if (isChangedUB) {
      solver->setColUpper(col, origBound);
    } else {
      solver->setColLower(col, origBound);
    }
    solver->resolve();
    solver->markHotStart();
    if (isChangedUB) {
      solver->setColUpper(col, newBound);
    } else {
      solver->setColLower(col, newBound);
    }
    solver->solveFromHotStart();
  }
  return (solver->isProvenOptimal());
} /* solveFromHotStart overload */

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
  int resolve_count = 0;
#ifdef USE_CLP
  OsiClpSolverInterface* clpsolver = NULL;
  try {
    clpsolver = dynamic_cast<OsiClpSolverInterface*>(solver);
  } catch (std::exception& e) {
    // Disregard / catch below
  }
  if (clpsolver && maxNumResolves > 0) {
    // If it is supposedly proven primal infeasible, might be good to do a resolve first
    // This is, e.g., something that arises with miplib2003/pp08aCUTS_presolved -32
    if (clpsolver->isProvenPrimalInfeasible()) {
      clpsolver->getModelPtr()->setNumberIterations(0);
      clpsolver->getModelPtr()->setMaximumSeconds(timeLimit);
      clpsolver->resolve();
      resolve_count++;
    }

    // First clean
    if (resolve_count < maxNumResolves) {
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
        resolve = (!isZero(curr_infeas) && (resolve_count == 0 || curr_infeas < infeas))
            || (resolve_count == 0);
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
    } // check resolve_count < max number of resolves
  } // check clpsolver exists and some strictly positive number of resolves allowed
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

/**
 * Generic way to apply a set of cuts to a solver
 * @return Solver value after adding the cuts, if they were added successfully
 */
double applyCutsCustom(OsiSolverInterface* const solver, const OsiCuts& cs,
    const int startCutIndex, const int numCutsOverload, double compare_obj) {
  const int numCuts = (numCutsOverload >= 0) ? numCutsOverload : cs.sizeCuts();
  if (solver && !solver->isProvenOptimal()) {
    solver->resolve();
    checkSolverOptimality(solver, false);
  }
  if (solver == NULL || !solver->isProvenOptimal()) {
    return solver->getInfinity();
  }

  const double initObjValue = solver->getObjValue();
  if (isNegInfinity(compare_obj)) {
    compare_obj = initObjValue;
  }

  // If the number of cuts is 0, we simply get the value without these cuts applied
  if (numCuts > 0) {
#ifdef TRACE
    OsiSolverInterface::ApplyCutsReturnCode code;
    int num_applied = 0;
#endif
    if ((startCutIndex == 0) && (numCuts == cs.sizeCuts())) {
#ifdef TRACE
      code = solver->applyCuts(cs);
      num_applied = code.getNumApplied();
#else
      solver->applyCuts(cs);
#endif
    } else {
#ifdef TRACE
      const int old_num_rows = solver->getNumRows();
#endif
      OsiRowCut* cuts;
      if (numCuts > 0) {
        cuts = new OsiRowCut[numCuts];

        for (int i = startCutIndex; i < startCutIndex + numCuts; i++) {
          cuts[i - startCutIndex] = cs.rowCut(i);
        }
        solver->applyRowCuts(numCuts, cuts);
        delete[] cuts;
      }
#ifdef TRACE
      num_applied = (solver->getNumRows() - old_num_rows);
#endif
    }

#ifdef TRACE
    printf("Cuts applied: %d (of the %d possible).\n",
        num_applied, numCuts);
#endif
  } /* check if num cuts > 0 */

  const bool useResolve = solver->isProvenOptimal();
  if (useResolve) {
    solver->resolve();
  } else {
    solver->initialSolve();
  }
  checkSolverOptimality(solver, false);

  // At this point it is useful to check for integer feasibility quickly
  // NB: If we reintroduce the subspace notion, then this could actually happen
  bool inSubspace = false;
  if (!inSubspace && !solver->isProvenOptimal() && !hitTimeLimit(solver) && !solver->isIterationLimitReached()) {
    try {
      error_msg(errstr,
          "Not in subspace and solver not optimal after adding cuts, so we may have an integer infeasible problem. Exit status: %d.\n",
          dynamic_cast<OsiClpSolverInterface*>(solver)->getModelPtr()->status());
    } catch (std::exception& e) {
      error_msg(errstr,
          "Not in subspace and solver not optimal after adding cuts, so we may have an integer infeasible problem.\n");
    }
#ifdef SAVE_INFEASIBLE_LP
    std::string infeas_f_name = "infeasible_lp";
    solver->writeMps(infeas_f_name.c_str(), "mps", 1);
#endif
    exit(1);
  }

  // Sometimes there are slight inaccuracies at the end that we can get rid of
  if (solver->isProvenOptimal() && solver->getObjValue() < compare_obj) {
    solver->resolve();
    checkSolverOptimality(solver, false);
  }

  if (solver->isProvenOptimal()) {
    return solver->getObjValue();
  } else {
    return solver->getInfinity();  // Minimization problem, so this means infeas
  }
} /* applyCutsCustom */

/**
 * Generic way to apply a set of cuts to a solver
 * @return Solver value after adding the cuts, if they were added successfully
 */
double applyCutsCustom(OsiSolverInterface* const solver, const OsiCuts& cs,
    const int numCutsOverload, double compare_obj) {
  const int numCuts = (numCutsOverload >= 0) ? numCutsOverload : cs.sizeCuts();
  return applyCutsCustom(solver, cs, 0, numCuts, compare_obj);
} /* applyCutsCustom */

/**
 * Get objective value from solver with cuts added
 */
double getObjValue(OsiSolverInterface* const solver, const OsiCuts* const cuts) {
  if (cuts != NULL && cuts->sizeCuts() > 0) {
    const int numRows = solver->getNumRows();
    const double obj = applyCutsCustom(solver, *cuts);
    solver->restoreBaseModel(numRows);
    return obj;
  } else {
    if (solver->isProvenOptimal()) {
      return solver->getObjValue();
    } else {
      if (solver->isProvenPrimalInfeasible()) {
        return solver->getInfinity();
      } else {
        return -1 * solver->getInfinity();
      }
    }
  }
} /* getObjValue */

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
