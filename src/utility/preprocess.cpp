/**
 * @file preprocess.cpp
 * @author A. M. Kazachkov
 * @date 2019-03-16
 */
#include "preprocess.hpp"

// Project files
#include "analysis.hpp"
#include "BBHelper.hpp"
#ifdef USE_GUROBI
#include "GurobiHelper.hpp"
#endif
#ifdef USE_CPLEX
#include "CplexHelper.hpp"
#endif
#include "SolverHelper.hpp"
#include "SolverInterface.hpp"
#include "VPCParameters.hpp"
using namespace VPCParametersNamespace;

// COIN-OR
#include <OsiSolverInterface.hpp>
#ifdef USE_CLP
#include <OsiClpSolverInterface.hpp>
#endif

const int countBoundInfoEntries = 2;
const int countSummaryBBInfoEntries = 4 * 2;
const int countOrigProbEntries = 13;
const int countCleanedProbEntries = 15;
const int countFullBBInfoEntries = static_cast<int>(BB_INFO_CONTENTS.size()) * 4 * 2;

/**
 * Perform preprocessing and get statistics
 */
void performCleaning(const VPCParametersNamespace::VPCParameters& params,
    OsiSolverInterface* const solver, const std::string& filename,
    const double ip_obj, const int CLEANING_MODE_OPTION, const char SEP) {
  FILE* logfile = params.logfile;
  if (logfile == NULL) {
    std::cerr << "Cleaning operations currently will not be performed when logfile == NULL.\n";
    return;
  }

  // Save strategy
  const int strategy = params.get(BB_STRATEGY);
  SummaryBBInfo orig_info;

  // Original instance info
//  const bool DO_STRONG_BRANCHING = false;
  const int origNumRows = solver->getNumRows();
  const int origNumCols = solver->getNumCols();
  const int origNumNonZero = solver->getMatrixByRow()->getNumElements();
  const double origLPOpt = solver->getObjValue();

  // Get original instance branching info
#ifdef TRACE
    printf("Performing branch-and-bound on original instance.\n");
#endif
  runBBTests(params, &orig_info, NULL, NULL,
      params.get(stringParam::FILENAME),
      solver, ip_obj, NULL, NULL);

  /********** Now we do the cleaning **********/
  SolverInterface* cleanedSolver = new SolverInterface;
  std::string presolved_name_stub =
      (CLEANING_MODE_OPTION <= 1) ? filename + "_presolved" : "";
  const std::string cleaned_name =
      (CLEANING_MODE_OPTION <= 1) ?
          presolved_name_stub : filename + "_cleaned";

  // First get presolved opt using commercial solver of choice
  double presolvedLPOpt;
  int numBoundsChanged = 0;
#ifdef USE_GUROBI
  if (use_bb_option(strategy, BB_Strategy_Options::gurobi)) {
#ifdef TRACE
    printf("\n## Presolve model with Gurobi. ##\n");
#endif
    presolveModelWithGurobi(params, strategy, params.get(stringParam::FILENAME).c_str(),
        presolvedLPOpt, presolved_name_stub, ip_obj); // returns mps.gz
    cleanedSolver->readMps(presolved_name_stub.c_str());

    // Check if any bounds were changed, but only if no rows/cols deleted
    // NB: this is being done to avoid redoing BB in case presolve did nothing
    // If the # of rows and cols are the same and all their bounds are the same...
    // it seems safe to assume that the instance has not changed
    if (solver->getNumRows() == cleanedSolver->getNumRows()) {
      for (int row = 0; row < cleanedSolver->getNumRows(); row++) {
        if (!isVal(solver->getRightHandSide()[row],
            cleanedSolver->getRightHandSide()[row])) {
          numBoundsChanged++;
        }
      }
    }
    if (solver->getNumCols() == cleanedSolver->getNumCols()) {
      for (int col = 0; col < cleanedSolver->getNumCols(); col++) {
        if (!isVal(solver->getColLower()[col],
            cleanedSolver->getColLower()[col])) {
          numBoundsChanged++;
        }
        if (!isVal(solver->getColUpper()[col],
            cleanedSolver->getColUpper()[col])) {
          numBoundsChanged++;
        }
      }
    }

    // Make sure we are doing a minimization problem; this is just to make later
    // comparisons simpler (i.e., a higher LP obj after adding the cut is better).
    if (cleanedSolver->getObjSense() < 1e-3) {
      printf(
          "\n## Detected maximization problem. Negating objective function to make it minimization. ##\n");
      cleanedSolver->setObjSense(1.0);
      const double* obj = cleanedSolver->getObjCoefficients();
      for (int col = 0; col < cleanedSolver->getNumCols(); col++) {
        cleanedSolver->setObjCoeff(col, -1. * obj[col]);
      }
      double objOffset = 0.;
      cleanedSolver->getDblParam(OsiDblParam::OsiObjOffset, objOffset);
      if (objOffset != 0.) {
        cleanedSolver->setDblParam(OsiDblParam::OsiObjOffset, -1. * objOffset);
      }
    }

    // Perform initial solve
    cleanedSolver->initialSolve();
    if (!checkSolverOptimality(cleanedSolver, false)) {
      error_msg(errorstring, "After initial solve, solver is not optimal.\n");
      writeErrorToLog(errorstring, params.logfile);
      exit(1);
    }

    if (CLEANING_MODE_OPTION > 1) {
      remove(presolved_name_stub.c_str()); // remove the temporary file
    }
  } // gurobi
#endif // use_gurobi
#ifdef USE_CPLEX
  if (use_bb_option(strategy, BB_Strategy_Options::cplex)) {
#ifdef TRACE
    printf("\n## Presolve model with CPLEX. ##\n");
#endif
    presolveModelWithCplexCallable(params, strategy, params.get(stringParam::FILENAME).c_str(),
        presolvedLPOpt, presolved_name_stub, ip_obj); // returns mps.gz
    cleanedSolver->readMps(presolved_name_stub.c_str());

    // Check if any bounds were changed, but only if no rows/cols deleted
    // NB: this is being done to avoid redoing BB in case presolve did nothing
    // If the # of rows and cols are the same and all their bounds are the same...
    // it seems safe to assume that the instance has not changed
    if (solver->getNumRows() == cleanedSolver->getNumRows()) {
      for (int row = 0; row < cleanedSolver->getNumRows(); row++) {
        if (!isVal(solver->getRightHandSide()[row],
            cleanedSolver->getRightHandSide()[row])) {
          numBoundsChanged++;
        }
      }
    }
    if (solver->getNumCols() == cleanedSolver->getNumCols()) {
      for (int col = 0; col < cleanedSolver->getNumCols(); col++) {
        if (!isVal(solver->getColLower()[col],
            cleanedSolver->getColLower()[col])) {
          numBoundsChanged++;
        }
        if (!isVal(solver->getColUpper()[col],
            cleanedSolver->getColUpper()[col])) {
          numBoundsChanged++;
        }
      }
    }

    // Make sure we are doing a minimization problem; this is just to make later
    // comparisons simpler (i.e., a higher LP obj after adding the cut is better).
    if (cleanedSolver->getObjSense() < 1e-3) {
      printf(
          "\n## Detected maximization problem. Negating objective function to make it minimization. ##\n");
      cleanedSolver->setObjSense(1.0);
      const double* obj = cleanedSolver->getObjCoefficients();
      for (int col = 0; col < cleanedSolver->getNumCols(); col++) {
        cleanedSolver->setObjCoeff(col, -1. * obj[col]);
      }
      double objOffset = 0.;
      cleanedSolver->getDblParam(OsiDblParam::OsiObjOffset, objOffset);
      if (objOffset != 0.) {
        cleanedSolver->setDblParam(OsiDblParam::OsiObjOffset, -1. * objOffset);
      }
    }

    // Perform initial solve
    cleanedSolver->initialSolve();
    if (!checkSolverOptimality(cleanedSolver, false)) {
      error_msg(errorstring, "After initial solve, solver is not optimal.\n");
      writeErrorToLog(errorstring, params.logfile);
      exit(1);
    }

    if (CLEANING_MODE_OPTION > 1) {
      remove(presolved_name_stub.c_str()); // remove the temporary file
    }
  } // cplex
#endif // use_cplex

  // Now (perhaps) clean using our own methods
  int numSBFixed = 0;
  if (CLEANING_MODE_OPTION > 1) {
    bool is_clean = false; // 2 (or higher) = do our own cleaning
    const int MAX_ITER = 5;
    int iter = 0;
    while (!is_clean && iter < MAX_ITER) {
#ifdef TRACE
      printf("Clean model with custom method (remove one-sided split disjunctions and tighten bounds, iteratively).\n");
#endif
      is_clean = cleanProblem(params, cleanedSolver, numBoundsChanged, numSBFixed);
      iter++;
    }
    cleanedSolver->writeMps(cleaned_name.c_str(), "mps", cleanedSolver->getObjSense());

    if (numBoundsChanged > 0 || numSBFixed > 0) {
      cleanedSolver->resolve();
      if (!checkSolverOptimality(cleanedSolver, false)) {
        error_msg(errorstring, "After specialized cleaning, cleanedSolver is not optimal.\n");
        writeErrorToLog(errorstring, params.logfile);
        exit(1);
      }
    }
  } // check cleaning mode > 1

  // Get new solver info, including primal and dual degeneracy values
  const int cleanedNumCols = cleanedSolver->getNumCols();
  const int cleanedNumRows = cleanedSolver->getNumRows();
  const int cleanedNumNonZero = cleanedSolver->getMatrixByRow()->getNumElements();
  const double cleanedLPOpt = (cleanedNumNonZero > 0) ? cleanedSolver->getObjValue() : ip_obj;

  // Set up new BBInfos
  SummaryBBInfo cleaned_info;
  const bool was_cleaned = (numSBFixed > 0) || (numBoundsChanged > 0)
      || (cleanedNumCols != origNumCols) || (cleanedNumRows != origNumRows)
      || (cleanedNumNonZero != origNumNonZero);
  if (cleanedNumNonZero > 0 && was_cleaned) {
#ifdef TRACE
    printf("\n## Performing branch-and-bound on cleaned instance. ##\n");
#endif
    runBBTests(params, &cleaned_info, NULL, NULL, cleaned_name + ".mps", cleanedSolver,
        ip_obj, NULL, NULL);
  } /* check if any cleaning was performed */
  else if (cleanedNumNonZero == 0) {
    presolvedLPOpt = ip_obj;
//    initializeBBInfo(cleanedBBInfo, ip_obj);
  } /* check if cleaning yielded an empty problem */
  else {
    cleaned_info = orig_info;
  } /* strong branching did nothing so do not repeat the experiments */

  // Now save everything to the logfile; instance name is already saved
  fprintf(logfile, "%d%c", strategy, SEP);
  { // BOUND INFO
    int count = 0;
    fprintf(logfile, "%s%c", stringValue(origLPOpt, "%2.20f").c_str(), SEP); count++;
    fprintf(logfile, "%s%c", stringValue(cleanedLPOpt, "%2.20f").c_str(), SEP); count++;
    assert( count == countBoundInfoEntries );
  }
  printSummaryBBInfo({orig_info, cleaned_info}, logfile);
  printOrigProbInfo(solver, logfile, SEP);
  printOrigProbInfo(cleanedSolver, logfile, SEP);
  fprintf(logfile, "%d%c", numBoundsChanged, SEP);
  fprintf(logfile, "%d%c", numSBFixed, SEP);
  printFullBBInfo({orig_info, cleaned_info}, logfile);

  if (cleanedSolver) { delete cleanedSolver; }
} /* performCleaning */

/**
 * Makes sure no variable bounds can be tightened,
 * including via strong branching
 */
bool cleanProblem(const VPCParametersNamespace::VPCParameters& params, OsiSolverInterface* solver,
    int& numBoundsChanged, int& numSBFixed) {
  bool is_clean = true;

  const int numCols = solver->getNumCols();
//  std::vector<double> origObjCoeff(solver->getObjCoefficients(),
//      solver->getObjCoefficients() + numCols);
  std::vector<double> origColSolution(solver->getColSolution(),
        solver->getColSolution() + numCols);
  std::vector<double> obj(numCols, 0.);

  // Set up solver for checking bounds
  OsiSolverInterface* boundSolver = solver->clone();
  boundSolver->setObjective(obj.data());
  boundSolver->setObjSense(solver->getObjSense());
  double objOffset = 0.;
  solver->getDblParam(OsiDblParam::OsiObjOffset, objOffset);
  boundSolver->setDblParam(OsiDblParam::OsiObjOffset, objOffset);

  // For strong branching
  try {
#ifdef USE_CLP
    setupClpForStrongBranching(dynamic_cast<OsiClpSolverInterface*>(solver));
#endif
  } catch (std::exception& e) {
    // continue anyway
  }
  solver->enableFactorization();
  solver->markHotStart();
  for (int col = 0; col < numCols; col++) {
    // If col is integral, ensure bounds are integral
    // If the integer variable is fractional, we also strong branch
    if (solver->isInteger(col)) {
      // Is the variable fractional in the solution?
      if (!isVal(origColSolution[col], std::floor(origColSolution[col]),
          params.get(doubleConst::AWAY))
          && !isVal(origColSolution[col], std::ceil(origColSolution[col]),
              params.get(doubleConst::AWAY))) {
        const double origLB = solver->getColLower()[col];
        const double origUB = solver->getColUpper()[col];
        bool downBranchFeasible = true, upBranchFeasible = true;

        // Check down branch
        solver->setColUpper(col, std::floor(origColSolution[col]));
        //solver->solveFromHotStart();
        solveFromHotStart(solver, col, true, origUB, std::floor(origColSolution[col]));
        if (solver->isProvenPrimalInfeasible()) {
          downBranchFeasible = false;
        }
        solver->setColUpper(col, origUB);

        // Return to original state
        solver->solveFromHotStart();

        // Check up branch
        solver->setColLower(col, std::ceil(origColSolution[col]));
        //solver->solveFromHotStart();
        solveFromHotStart(solver, col, false, origLB, std::ceil(origColSolution[col]));
        if (solver->isProvenPrimalInfeasible()) {
          upBranchFeasible = false;
        }
        solver->setColLower(col, origLB);

        // Return to original state
        solver->solveFromHotStart();

        // Check if some side of the split is infeasible
        if (!downBranchFeasible || !upBranchFeasible) {
          if (!downBranchFeasible && !upBranchFeasible) {
            // Infeasible problem
            error_msg(errorstring,
                "Infeasible problem due to integer variable %d (value %e).\n",
                col, origColSolution[col]);
            writeErrorToLog(errorstring, params.logfile);
            exit(1);
          }
          numSBFixed++;
          is_clean = false;
          if (!downBranchFeasible) {
            boundSolver->setColLower(col, std::ceil(origColSolution[col]));
          } else {
            boundSolver->setColUpper(col, std::floor(origColSolution[col]));
          }
        }
      } /* end strong branching */

      // Check integrality of LB
      if (!isVal(solver->getColLower()[col],
          std::floor(solver->getColLower()[col]))) {
        numBoundsChanged++;
        is_clean = false;
        boundSolver->setColLower(col, std::ceil(solver->getColLower()[col]));
      }

      // Check integrality of UB
      if (!isVal(solver->getColUpper()[col],
          std::ceil(solver->getColUpper()[col]))) {
        numBoundsChanged++;
        is_clean = false;
        boundSolver->setColUpper(col, std::floor(solver->getColUpper()[col]));
      }
    } /* check if integer */

    // Check if LB can be tightened
    boundSolver->setObjCoeff(col, 1.);
    boundSolver->resolve();
    if (checkSolverOptimality(boundSolver, false)) {
      if (greaterThanVal(boundSolver->getColSolution()[col],
          solver->getColLower()[col])) {
        numBoundsChanged++;
        is_clean = false;
        boundSolver->setColLower(col, boundSolver->getColSolution()[col]);
      }
    }

    // Check if UB can be tightened
    boundSolver->setObjCoeff(col, -1.);
    boundSolver->resolve();
    if (checkSolverOptimality(boundSolver, false)) {
      if (lessThanVal(boundSolver->getColSolution()[col],
          solver->getColUpper()[col])) {
        numBoundsChanged++;
        is_clean = false;
        boundSolver->setColUpper(col, boundSolver->getColSolution()[col]);
      }
    }

    // Reset
    boundSolver->setObjCoeff(col, 0.);
  } /* end iterating over columns */
  solver->unmarkHotStart();
  solver->disableFactorization();

  // Change upper and lower bounds wherever needed
  for (int col = 0; col < solver->getNumCols(); col++) {
    solver->setColLower(col, boundSolver->getColLower()[col]);
    solver->setColUpper(col, boundSolver->getColUpper()[col]);
  } /* end iterating over columns to update solver */

  if (boundSolver) {
    delete boundSolver;
  }

  return is_clean;
} /* cleanProblem */

void printPreprocessingHeader(const VPCParametersNamespace::VPCParameters& params, const char SEP) {
  FILE* logfile = params.logfile;
  if (logfile == NULL)
    return;

  // First line of the header details the categories of information displayed
  std::string tmpstring = "";
  fprintf(logfile, "%c", SEP); // instance name
  fprintf(logfile, "%c", SEP); // bb strategy
  fprintf(logfile, "%s", "BOUND INFO");
  tmpstring.assign(countBoundInfoEntries, SEP);
  fprintf(logfile, "%s", tmpstring.c_str());
  fprintf(logfile, "%s", "BB INFO");
  tmpstring.assign(countSummaryBBInfoEntries, SEP);
  fprintf(logfile, "%s", tmpstring.c_str());
  fprintf(logfile, "%s", "ORIG PROB INFO");
  tmpstring.assign(countOrigProbEntries, SEP);
  fprintf(logfile, "%s", tmpstring.c_str());
  fprintf(logfile, "%s", "CLEANED PROB INFO");
  tmpstring.assign(countCleanedProbEntries, SEP);
  fprintf(logfile, "%s", tmpstring.c_str());
  fprintf(logfile, "%s", "FULL BB INFO");
  tmpstring.assign(countFullBBInfoEntries, SEP);
  fprintf(logfile, "%s", tmpstring.c_str());
  fprintf(logfile, "%s", "END");
  fprintf(logfile, "\n");

  fprintf(logfile, "%s%c", "INSTANCE", SEP);
  fprintf(logfile, "%s%c", "STRATEGY", SEP);
  { // BOUND INFO
    int count = 0;
    fprintf(logfile, "%s%c", "ORIG LP OBJ", SEP); count++;  // 1
    fprintf(logfile, "%s%c", "CLEANED LP OBJ", SEP); count++;  // 1
    assert(count == countBoundInfoEntries);
  }
  { // BB INFO
    int count = 0;
    std::vector<std::string> nameVec = {"NODES", "TIME"};
    for (auto name : nameVec) {
      fprintf(logfile, "%s%c", ("ORIG FIRST GUR " + name).c_str(), SEP); count++;
      fprintf(logfile, "%s%c", ("CLEANED FIRST GUR " + name).c_str(), SEP); count++;
      fprintf(logfile, "%s%c", ("ORIG BEST GUR " + name).c_str(), SEP); count++;
      fprintf(logfile, "%s%c", ("CLEANED BEST GUR " + name).c_str(), SEP); count++;
    }
    assert(count == countSummaryBBInfoEntries);
  } // BB INFO
  { // ORIG PROB
    int count = 0;
    fprintf(logfile, "%s%c", "ORIG ROWS", SEP); count++;
    fprintf(logfile, "%s%c", "ORIG COLS", SEP); count++;
    fprintf(logfile, "%s%c", "ORIG NUM FRAC", SEP); count++;
    fprintf(logfile, "%s%c", "ORIG MIN FRACTIONALITY", SEP); count++;
    fprintf(logfile, "%s%c", "ORIG MAX FRACTIONALITY", SEP); count++;
    fprintf(logfile, "%s%c", "ORIG EQ ROWS", SEP); count++;
    fprintf(logfile, "%s%c", "ORIG BOUND ROWS", SEP); count++;
    fprintf(logfile, "%s%c", "ORIG ASSIGN ROWS", SEP); count++;
    fprintf(logfile, "%s%c", "ORIG FIXED COLS", SEP); count++;
    fprintf(logfile, "%s%c", "ORIG GEN INT", SEP); count++;
    fprintf(logfile, "%s%c", "ORIG BINARY", SEP); count++;
    fprintf(logfile, "%s%c", "ORIG CONTINUOUS", SEP); count++;
    fprintf(logfile, "%s%c", "ORIG A-DENSITY", SEP); count++;
    assert(count == countOrigProbEntries);
  } // ORIG PROB
  { // CLEANED PROB
    int count = 0;
    fprintf(logfile, "%s%c", "CLEANED ROWS", SEP); count++;
    fprintf(logfile, "%s%c", "CLEANED COLS", SEP); count++;
    fprintf(logfile, "%s%c", "CLEANED NUM FRAC", SEP); count++;
    fprintf(logfile, "%s%c", "CLEANED MIN FRACTIONALITY", SEP); count++;
    fprintf(logfile, "%s%c", "CLEANED MAX FRACTIONALITY", SEP); count++;
    fprintf(logfile, "%s%c", "CLEANED EQ ROWS", SEP); count++;
    fprintf(logfile, "%s%c", "CLEANED BOUND ROWS", SEP); count++;
    fprintf(logfile, "%s%c", "CLEANED ASSIGN ROWS", SEP); count++;
    fprintf(logfile, "%s%c", "CLEANED FIXED COLS", SEP); count++;
    fprintf(logfile, "%s%c", "CLEANED GEN INT", SEP); count++;
    fprintf(logfile, "%s%c", "CLEANED BINARY", SEP); count++;
    fprintf(logfile, "%s%c", "CLEANED CONTINUOUS", SEP); count++;
    fprintf(logfile, "%s%c", "CLEANED A-DENSITY", SEP); count++;
    fprintf(logfile, "%s%c", "CLEANED BOUNDS CHANGED", SEP); count++;
    fprintf(logfile, "%s%c", "CLEANED NUM SB FIXED", SEP); count++;
    assert(count == countCleanedProbEntries);
  } // CLEANED PROB
  { // FULL BB INFO
    int count = 0;
    for (std::string name : BB_INFO_CONTENTS) {
      fprintf(logfile, "%s%c", ("ORIG FIRST GUR " + name).c_str(), SEP); count++;
      fprintf(logfile, "%s%c", ("CLEANED FIRST GUR " + name).c_str(), SEP); count++;
      fprintf(logfile, "%s%c", ("ORIG BEST GUR " + name).c_str(), SEP); count++;
      fprintf(logfile, "%s%c", ("CLEANED BEST GUR " + name).c_str(), SEP); count++;
      fprintf(logfile, "%s%c", ("ORIG AVG GUR " + name).c_str(), SEP); count++;
      fprintf(logfile, "%s%c", ("CLEANED AVG GUR " + name).c_str(), SEP); count++;
    }
    for (std::string name : BB_INFO_CONTENTS) {
      fprintf(logfile, "%s%c", ("ORIG ALL GUR " + name).c_str(), SEP); count++;
    }
    for (std::string name : BB_INFO_CONTENTS) {
      fprintf(logfile, "%s%c", ("CLEANED ALL GUR " + name).c_str(), SEP); count++;
    }
    assert(count == countFullBBInfoEntries);
  } // FULL BB INFO

  fprintf(logfile, "%s%c", "end_time_string", SEP);
  fprintf(logfile, "%s%c", "time elapsed", SEP);
  fprintf(logfile, "\n");
  fflush(logfile);
} /** printPreprocessingHeader */
