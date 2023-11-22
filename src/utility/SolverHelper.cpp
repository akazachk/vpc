#include "SolverHelper.hpp"

#include "SolverInterface.hpp"
#include "utility.hpp"

#include <OsiCuts.hpp>

void initializeSolver(
    OsiSolverInterface* &solver,
    const std::string FILENAME,
    const int VERBOSITY,
    const double TIMELIMIT,
    FILE* logfile) {
  solver = new SolverInterface;
  setLPSolverParameters(solver, VERBOSITY, TIMELIMIT);

  std::string dir, instname, in_file_ext;
  if (parseFilename(dir, instname, in_file_ext, FILENAME, logfile) != 0) {
    error_msg(errorstring,
        "Unable to parse filename: %s. Found: dir=\"%s\", instname=\"%s\",ext=\"%s\".\n",
        FILENAME.c_str(), dir.c_str(),
        instname.c_str(), in_file_ext.c_str());
    exit(1);
  }

  int status = 0;
  if (in_file_ext.compare("lp") == 0) {
#ifdef TRACE
    printf("\n## Reading LP file. ##\n");
#endif
    status = solver->readLp(FILENAME.c_str());
  } else {
    if (in_file_ext.compare("mps") == 0) {
#ifdef TRACE
      printf("\n## Reading MPS file. ##\n");
#endif
      status = solver->readMps(FILENAME.c_str());
    } else {
      try {
#ifdef TRACE
        printf("\n## Reading MPS file. ##\n");
#endif
        status = solver->readMps(FILENAME.c_str());
      } catch (std::exception& e) {
        error_msg(errorstring, "Unrecognized extension: %s.\n",
            in_file_ext.c_str());
        writeErrorToLog(errorstring, logfile);
        exit(1);
      }
    }
  } // read file
  if (status < 0) {
    error_msg(errorstring, "Unable to read in file %s.\n",
        FILENAME.c_str());
    writeErrorToLog(errorstring, logfile);
    exit(1);
  }

  // Make sure we are doing a minimization problem; this is just to make later
  // comparisons simpler (i.e., a higher LP obj after adding the cut is better).
  if (solver->getObjSense() < 1e-3) {
    printf(
        "\n## Detected maximization problem. Negating objective function to make it minimization. ##\n");
    solver->setObjSense(1.0);
    const double* obj = solver->getObjCoefficients();
    for (int col = 0; col < solver->getNumCols(); col++) {
      solver->setObjCoeff(col, -1. * obj[col]);
    }
    double objOffset = 0.;
    solver->getDblParam(OsiDblParam::OsiObjOffset, objOffset);
    if (objOffset != 0.) {
      solver->setDblParam(OsiDblParam::OsiObjOffset, -1. * objOffset);
    }
  }
} /* initializeSolver */

double getObjOffset(const OsiSolverInterface* const solver) {
  double offset = 0.;
  solver->getDblParam(OsiDblParam::OsiObjOffset, offset);
  return offset;
} /* getObjOffset */

#ifdef USE_CBC
#include <CbcModel.hpp>

void setIPSolverParameters(CbcModel* const cbc_model, const int verbosity) {
  if (verbosity > 0) {
    cbc_model->setLogLevel(verbosity);
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
 * @details Sets OsiMaxNumIterationHotStart to \p hot_start_iter_limit and special options to 16 (use strong branching rather than Clp's).
 * We need to be careful with the strong branching options;
 * sometimes the Clp strong branching fails, such as with arki001, branching down on variable 924.
 */
void setupClpForStrongBranching(OsiClpSolverInterface* const solver, const int hot_start_iter_limit) {
  solver->setIntParam(OsiMaxNumIterationHotStart, hot_start_iter_limit);
  solver->setSpecialOptions(16); // use standard strong branching rather than clp's
} /* setupClpForStrongBranching */

/**
 * @details We want to run full strong branching and enable the fixing of variables.
 */
void setupClpForCbc(OsiSolverInterface* const solver,
    const int verbosity, const double max_time,
    const int hot_start_iter_limit) {
  setLPSolverParameters(solver, verbosity, max_time);
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

/**
 * @details We do not zero out the other coefficients unless requested
 */
void addToObjectiveFromPackedVector(
    /// [in,out] solver object to be modified
    OsiSolverInterface* const solver,
    /// [in] add values in vec to the objective
    const CoinPackedVectorBase* vec,
    /// [in] set all obj coeff to 0 before proceeding
    const bool zeroOut,
    /// [in] a constant multiple to scale the coefficients by (for example, we can set \p SHOULD_SCALE to false and pass the inverse of the 2-norm of \p vec for \p mult)
    const double mult,
    /// [in] assumed sorted in increasing order, in case we want to change a specific set of indices
    const std::vector<int>* const nonZeroColIndices,
    /// [in] add a scaled version of \p vec, by its two norm
    const bool SHOULD_SCALE) {
  if (zeroOut) {
    for (int i = 0; i < solver->getNumCols(); i++) {
      solver->setObjCoeff(i, 0.);
    }
  }

  if (vec) {
    double twoNorm = 1.;
    if (!isZero(mult) && SHOULD_SCALE) {
      const double realTwoNorm = vec->twoNorm() / solver->getNumCols();
      if (greaterThanVal(realTwoNorm, 0.)) {
        twoNorm = realTwoNorm;
      }
    }
    if (nonZeroColIndices) {
      int ind = 0; // assuming nonZeroColIndices are sorted in increasing order
      const int numNZ = (*nonZeroColIndices).size();
      for (int i = 0; i < (*vec).getNumElements(); i++) {
        const int col = (*vec).getIndices()[i];

        while ((ind < numNZ) && (col > (*nonZeroColIndices)[ind])) {
          ind++;
        }
        if (ind >= numNZ)
          break;
        if (col < (*nonZeroColIndices)[ind])
          continue;

        const double elem = (*vec).getElements()[i];
        const double coeff = solver->getObjCoefficients()[ind];
        solver->setObjCoeff(ind, isZero(mult) ? 0. : coeff + elem * mult / twoNorm);
      }
    } else {
      for (int i = 0; i < (*vec).getNumElements(); i++) {
        const int col = (*vec).getIndices()[i];
        const double elem = (*vec).getElements()[i];
        const double coeff = solver->getObjCoefficients()[col];
        solver->setObjCoeff(col, isZero(mult) ? 0. : coeff + elem * mult / twoNorm);
      }
    }
  }
} /* addToObjectiveFromPackedVector */

/**
 * @details Only use specified indices if both \p numIndices > 0 and \p indices != NULL
 */
void setConstantObjectiveFromPackedVector(
    /// [in,out] solver to be modified
    OsiSolverInterface* const solver, 
    /// [in] value to set the coefficients to (e.g., 0)
    const double val,
    /// [in] if specified, only the indices in \p indices will be set
    const int numIndices,
    /// [in] if specified, only these indices will be set
    const int* indices) {
  if (numIndices > 0 && indices) {
    for (int i = 0; i < numIndices; i++) {
      solver->setObjCoeff(indices[i], val);
    }
  } else {
    for (int i = 0; i < solver->getNumCols(); i++) {
      solver->setObjCoeff(i, val);
    }
  }
} /* setConstantObjectiveFromPackedVector */

/// @details Set solver solution to specified value (if \p sol is NULL, set all to to zero)
void setSolverSolution(OsiSolverInterface* const solver, const double* const sol) {
  if (sol) {
    solver->setColSolution(sol);
  } else {
    std::vector<double> zerosol(solver->getNumCols(), 0.);
    solver->setColSolution(zerosol.data());
  }
} /* setSolverSolution */

/// @details Some problems with normal OsiSolverInterface::solverFromHotStart() can happen,
/// such as with arki001, in which we can have OsiSolverInterface::isIterationLimitReached().
/// The plan then is to try one more time with a hot start, and otherwise switch to resolving
/// (depending on the parameter \p MAX_NUM_HOT_START_VIOLS).
///
/// @return OsiSolverInterface::isProvenOptimal()
bool solveFromHotStart(
    /// [in,out] Original solver
    OsiSolverInterface* const solver,
    /// [in] Which variable is being modified
    const int col,
    /// [in] Is the upper bound changing?
    const bool isChangedUB,
    /// [in] Old value of the variable bound
    const double origBound,
    /// [in] New value of the variable bound
    const double newBound,
    /// [in,out] If in a prior call, hot start seemed to not work, we might disable it, in which case we need to use normal resolving
    int& numHotStartViolations,
    /// [in] Maximum number of \p numHotStartViolations before switching to resolve
    const int MAX_NUM_HOT_START_VIOLS) {
  if (numHotStartViolations < MAX_NUM_HOT_START_VIOLS) {
    solver->solveFromHotStart();
    if (solver->isIterationLimitReached()) {
      numHotStartViolations++;
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
      if (solver->isIterationLimitReached()) {
        numHotStartViolations++;
      }
    }
  } // check if hot start is not disabled

  if (numHotStartViolations >= MAX_NUM_HOT_START_VIOLS) {
    // Else, hot start is disabled, revert to resolve
    solver->resolve();
  }
  return (solver->isProvenOptimal());
} /* solveFromHotStart overload */

/**
 * @details
 * Resolve if there are inaccuracies detected
 * (e.g., the secondary status of the solver is not optimal,
 * or if the sum of the primal or dual infeasibilities is not zero)
 *
 * Something that can go wrong (e.g., bc1 -64 sb5 tmp_ind = 14):
 *  \p solver is declared optimal for the scaled problem but there are primal or dual infeasibilities in the unscaled problem.
 *  In this case, secondary status will be 2 or 3.
 * Sometimes cleanup helps.
 * Sometimes things break after enabling factorization (e.g., secondary status = 0, but sum infeasibilities > 0).
 *
 * @return OsiSolverInterface::isProvenOptimal()
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
 * @details If USE_CLP macro is defined and \p resolveFlag > 0, then we sometimes might want to resolve,
 * because occasionally the act of enabling factorization creates a bad (primal or dual) solution
 *
 * @return OsiSolverInterface::isProvenOptimal()
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
 * @return Solver value after adding the cuts, if they were added successfully
 */
double applyCutsCustom(OsiSolverInterface* const solver, const OsiCuts& cs,
    FILE* logfile, const int startCutIndex, const int numCutsOverload, double compare_obj) {
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
      writeErrorToLog(errstr, logfile);
    } catch (std::exception& e) {
      error_msg(errstr,
          "Not in subspace and solver not optimal after adding cuts, so we may have an integer infeasible problem.\n");
      writeErrorToLog(errstr, logfile);
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
 * @details Calls #applyCutsCustom(OsiSolverInterface* const, const OsiCuts*, FILE*, const int, const int, double)
 * with \p numCuts set to size of cut collection in \p cs, unless \p numCutsOverload takes a nonnegative value
 * @return Solver value after adding the cuts, if they were added successfully
 */
double applyCutsCustom(OsiSolverInterface* const solver, const OsiCuts& cs,
    FILE* logfile, const int numCutsOverload, double compare_obj) {
  const int numCuts = (numCutsOverload >= 0) ? numCutsOverload : cs.sizeCuts();
  return applyCutsCustom(solver, cs, logfile, 0, numCuts, compare_obj);
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

// <!---------------- Aliases from solver ---------------->
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
} /* isBasicCol */

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
} /* isNonBasicFreeCol */

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
} /* isNonBasicUBCol */

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
} /* isNonBasicLBCol */

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
} /* isNonBasicFixedCol */

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
} /* isBasicSlack */

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
} /* isNonBasicFreeSlack */

/// @details NOTE: We do at *lower* bound on purpose because of comment above:
/// the underlying model flips how it holds these with respect to Clp
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
} /* isNonBasicUBSlack */

/// @details NOTE: We do at *upper* bound on purpose because of comment above:
/// the underlying model flips how it holds these with respect to Clp
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
} /* isNonBasicLBSlack */

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
} /* isNonBasicFixedSlack */
