/**
 * @file preprocess.cpp
 * @author A. M. Kazachkov
 * @date 2019-03-16
 */
#include "preprocess.hpp"

//#include <filesystem> // for std::copy, C++17 introduces std::copy_file

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
const int countVersionInfoEntries = 4;
const int countExtraInfoEntries = 4;

/**
 * @details Do preprocessing (by CPLEX or Gurobi, as chosen), 
 * as well as possibly using cleanProblem() if \p CLEANING_MODE > 1,
 * and print information at end
 */
void performCleaning(
    const VPCParametersNamespace::VPCParameters& params,
    const OsiSolverInterface* const solver, 
    /// cleaned instance will be saved to filename_stub + "_presolved"/"_cleaned"
    /// (depending on if \p CLEANING_MODE = 1 or 2)
    const std::string& filename_stub,
    /// IP optimal value
    const double ip_obj,
    /// 0: no cleaning, 1: clean using preprocessing of CPLEX/Gurobi, 2: also clean with cleanProblem()
    const int CLEANING_MODE,
    /// when printing to log, what is the separation character to use (default is comma-separated)
    const char SEP) {
  FILE* logfile = params.logfile;
  if (logfile == NULL) {
    error_msg(errorstring,
        "Cleaning operations currently will not be performed when logfile == NULL.\n");
    writeErrorToLog(errorstring, stdout);
    exit(1);
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
  printf("\n## Performing branch-and-bound on original instance. ##\n");
  runBBTests(params, &orig_info, NULL, NULL,
      params.get(stringParam::FILENAME),
      solver, ip_obj, NULL, NULL);

  //********* Now we do the cleaning **********
  SolverInterface* cleanedSolver = new SolverInterface;
  setLPSolverParameters(cleanedSolver, params.get(VERBOSITY), params.get(TIMELIMIT));
  std::string presolved_name_stub =
      (CLEANING_MODE <= 1) ? filename_stub + "_presolved" : "";
  const std::string cleaned_name_stub =
      (CLEANING_MODE <= 1) ?
          presolved_name_stub : filename_stub + "_cleaned";

  // First get presolved opt using commercial solver of choice
  double presolvedLPOpt;
  int numBoundsChanged = 0;
  if (use_bb_option(strategy, BB_Strategy_Options::gurobi)) {
#ifdef USE_GUROBI
    printf("\n## Presolve model with Gurobi. ##\n");
    presolveModelWithGurobi(params, strategy, params.get(stringParam::FILENAME).c_str(),
        presolvedLPOpt, presolved_name_stub, ip_obj); // returns mps.gz
#else
    error_msg(errorstring, "Requesting to use Gurobi, but macro USE_GUROBI not set.\n");
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
#endif // use_gurobi
  } // gurobi
  if (use_bb_option(strategy, BB_Strategy_Options::cplex)) {
#ifdef USE_CPLEX
    printf("\n## Presolve model with CPLEX. ##\n");
    presolveModelWithCplexCallable(params, strategy, params.get(stringParam::FILENAME).c_str(),
        presolvedLPOpt, presolved_name_stub, ip_obj); // returns mps.gz
#else
    error_msg(errorstring, "Requesting to use CPLEX, but macro USE_CPLEX not set.\n");
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
#endif // use_cplex
  } // cplex

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
    error_msg(errorstring, "After initial solve, presolved model is not optimal.\n");
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }
  presolvedLPOpt = cleanedSolver->getObjValue();

  // Now (perhaps) clean using our own methods
  int numSBFixed = 0;
  if (CLEANING_MODE > 1) {
    remove(presolved_name_stub.c_str()); // remove the temporary file
    bool is_clean = false; // 2 (or higher) = do our own cleaning
    const int MAX_ITER = 5;
    int iter = 0;
    while (!is_clean && iter < MAX_ITER) {
      printf("\n## Clean model with custom method (remove one-sided split disjunctions and tighten bounds, iteratively). ##\n");
      is_clean = cleanProblem(params, cleanedSolver, numBoundsChanged, numSBFixed);
      iter++;
    }
    cleanedSolver->writeMps(cleaned_name_stub.c_str(), "mps", cleanedSolver->getObjSense());

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
    printf("\n## Performing branch-and-bound on cleaned instance. ##\n");
    VPCParametersNamespace::VPCParameters cleaned_params = params;
    cleaned_params.set(stringParam::FILENAME, cleaned_name_stub + ".mps.gz");
    std::string solfile = params.get(stringParam::SOLFILE);
    if (!solfile.empty()) {
      std::string dir, instname, in_ext;
      parseFilename(dir, instname, in_ext, solfile, params.logfile);
      if ((in_ext.compare(".sol") == 0) || (in_ext.compare(".sol.gz") == 0)
          || (in_ext.compare(".mst") == 0) || (in_ext.compare(".mst.gz") == 0)) {
        solfile = dir + "/"  + instname + ((CLEANING_MODE <= 1) ? "_presolved" : "_cleaned") + in_ext;
        cleaned_params.set(stringParam::SOLFILE, solfile);
      }
    }
    runBBTests(cleaned_params, &cleaned_info, NULL, NULL,
        cleaned_name_stub + ".mps",
        cleanedSolver, ip_obj, NULL, NULL);
  } // check if any cleaning was performed
  else if (cleanedNumNonZero == 0) {
    warning_msg(warnstring, "Presolve left an empty problem.\n");
    presolvedLPOpt = ip_obj;
//    initializeBBInfo(cleanedBBInfo, ip_obj);
  } // check if cleaning yielded an empty problem
  else {
    printf("Presolve did not change instance.\n");
    cleaned_info = orig_info;

    // Copy over the mip start / solution information (TODO)
    if (use_temp_option(params.get(intParam::TEMP), TempOptions::SAVE_IP_OPT)) {
      std::string solver = "";
#ifdef USE_GUROBI
      if (use_bb_option(strategy, BB_Strategy_Options::gurobi)) {
        solver = "_gurobi";
      }
#endif
#ifdef USE_CPLEX
      if (use_bb_option(strategy, BB_Strategy_Options::cplex)) {
        solver = "_cplex";
      }
#endif
      const std::string solfile = filename_stub + solver + ".mst.gz";
      const std::string solfile_cleaned = cleaned_name_stub + solver + ".mst.gz";
      if (fexists(solfile.c_str())) {
//        std::copy_file(solfile.c_str(), solfile_cleaned.c_str()); // needs C++17
      }
    }
  } // strong branching did nothing so do not repeat the experiments

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
 * @details Do strong branching to check if any variable bounds can be fixed
 *
 * @return true if no cleaning was done
 */
bool cleanProblem(
    const VPCParametersNamespace::VPCParameters& params, 
    /// [in,out] instance to be "cleaned"; outputted solver may have some bounds that are changed
    OsiSolverInterface* solver,
    /// [out] number of variable bounds that are changed (including those that are fixed)
    int& numBoundsChanged, 
    /// [out] number of variable bounds that are fixed after strong branching
    int& numSBFixed) {
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
      } // end strong branching

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
    } // check if integer

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
  } // end iterating over columns
  solver->unmarkHotStart();
  solver->disableFactorization();

  // Change upper and lower bounds wherever needed
  for (int col = 0; col < solver->getNumCols(); col++) {
    solver->setColLower(col, boundSolver->getColLower()[col]);
    solver->setColUpper(col, boundSolver->getColUpper()[col]);
  } // end iterating over columns to update solver

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
  fprintf(logfile, "%s", "VERSION INFO");
  tmpstring.assign(countVersionInfoEntries, SEP);
  fprintf(logfile, "%s", tmpstring.c_str());
  fprintf(logfile, "%s", "WRAP UP INFO");
  tmpstring.assign(countExtraInfoEntries, SEP);
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
  { // VERSION INFO
    fprintf(logfile, "%s%c", "vpc_version", SEP);
    fprintf(logfile, "%s%c", "cbc_version", SEP);
    fprintf(logfile, "%s%c", "clp_version", SEP);
    fprintf(logfile, "%s%c", "gurobi_version", SEP);
    fprintf(logfile, "%s%c", "cplex_version", SEP);
  } // VERSION INFO
  { // WRAP UP INFO
    fprintf(logfile, "%s%c", "ExitReason", SEP);
    fprintf(logfile, "%s%c", "end_time_string", SEP);
    fprintf(logfile, "%s%c", "time elapsed", SEP);
    fprintf(logfile, "%s%c", "instname", SEP);
  } // WRAP UP INFO
  fprintf(logfile, "\n");
  fflush(logfile);
} /* printPreprocessingHeader */
