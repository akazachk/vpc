/**
 * @file StrongBranchingDisjunction.cpp
 * @author Sean Kelley
 * @date 2023-04-20
 */
#include "StrongBranchingDisjunction.hpp"

// standard library
#include <algorithm> // min
#include <cmath> // round, power
#include <math.h> // abs, log2

// Project files
#include "CglVPC.hpp" // get timer information too
#include "OsiChooseStrongCustom.hpp"
#include "SolverHelper.hpp"
#include "utility.hpp"
#include "VPCEventHandler.hpp"
#include "VPCParameters.hpp"

#ifdef TRACE
#include "vpc_debug.hpp"
#endif

#ifdef USE_CLP
#include <OsiClpSolverInterface.hpp>
#endif

/****************** PUBLIC  **********************/
/** Param constructor */
StrongBranchingDisjunction::StrongBranchingDisjunction(const VPCParametersNamespace::VPCParameters& param) {
  initialize(NULL, &param);
} /* param constructor */

/** Copy and param constructor */
StrongBranchingDisjunction::StrongBranchingDisjunction(const StrongBranchingDisjunction& source, const VPCParametersNamespace::VPCParameters& param) {
  initialize(&source, &param);
} /* copy & param constructor */

/** Default constructor */
StrongBranchingDisjunction::StrongBranchingDisjunction() {
  initialize(NULL, NULL);
} /* default constructor */

/** Copy constructor */
StrongBranchingDisjunction::StrongBranchingDisjunction(const StrongBranchingDisjunction& source) {
  initialize(&source, NULL);
} /* copy constructor */

/** Destructor */
StrongBranchingDisjunction::~StrongBranchingDisjunction() {
} /* destructor */

/** Assignment operator */
StrongBranchingDisjunction& StrongBranchingDisjunction::operator=(const StrongBranchingDisjunction& source) {
  if (this != &source) {
    initialize(&source, NULL);
  }
  return *this;
} /* assignment operator */

/** Clone */
StrongBranchingDisjunction* StrongBranchingDisjunction::clone() const {
  return new StrongBranchingDisjunction(*this);
} /* clone */

/**
 * Set up the disjunction class as new
 * (except the timer pointer, and do not reset params;
 * in addition, here we do not reset var, in case user has set it)
 */
void StrongBranchingDisjunction::setupAsNew() {
  VPCDisjunction::setupAsNew();
} /* setupAsNew */

/**
 * @brief Prepare a new disjunction
 *
 * This will throw away all the information from the old disjunction, except it will not reset the timer
 */
DisjExitReason StrongBranchingDisjunction::prepareDisjunction(const OsiSolverInterface* const si) {
  // Reset things in case we are reusing the class for some reason
  setupAsNew();

  // start the timer
  if (timer)
    timer->start_timer(CglVPC::VPCTimeStatsName[static_cast<int>(CglVPC::VPCTimeStats::DISJ_GEN_TIME)]);

  // get number of variables to select
  int numTerms = params.get(VPCParametersNamespace::intParam::DISJ_TERMS);
  if (!isVal(std::log2(numTerms), std::round(std::log2(numTerms)), 1e-6)) {
    error_msg(errorstring,
        "Number of terms in disjunction must be a power of 2.\n");
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }
  int numVarsToSelect = std::round(std::log2(numTerms));
  DisjExitReason retVal = DisjExitReason::UNKNOWN;

  // get info from strong branching
  std::vector<int> fracIntVars;
  std::vector<double> downBranchObjs;
  std::vector<double> upBranchObjs;
  std::vector<double> scores;
  SolverInterface* solver;
  OsiChooseStrongCustom metric;
  try {
    solver = dynamic_cast<SolverInterface*>(si->clone());
  } catch (std::exception& e) {
    error_msg(errorstring,
        "Unable to clone solver into desired SolverInterface.\n");
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }
#ifdef USE_CLP
  try {
    setupClpForStrongBranching(dynamic_cast<OsiClpSolverInterface*>(solver));
  } catch (std::exception& e) {
    // It's okay, we can continue
  }
#endif

  // make sure solver is optimal
  if (!solver->isProvenOptimal()) {
    solver->initialSolve();
    if (!solver->isProvenOptimal()) {
      error_msg(errorstring, "Solver must be optimal to create disjunction.\n");
      writeErrorToLog(errorstring, params.logfile);
      exit(1);
    }
  }

  solver->enableFactorization();
  solver->markHotStart();
  const double * origColSolution = solver->getColSolution();
  int numberColumns = solver->getNumCols();
  double origObjVal = solver->getObjValue();

  for (int col = 0; col < numberColumns; col++) {
    // Is this an integer variable with fractional value?
    if (solver->isInteger(col) &&
        !isVal(origColSolution[col], std::floor(origColSolution[col]),
               params.get(VPCParametersNamespace::doubleConst::AWAY)) &&
        !isVal(origColSolution[col], std::ceil(origColSolution[col]),
               params.get(VPCParametersNamespace::doubleConst::AWAY))) {

      const double origLB = solver->getColLower()[col];
      const double origUB = solver->getColUpper()[col];
      bool downBranchFeasible = true, upBranchFeasible = true;
      fracIntVars.push_back(col);

      // Check down branch
      solver->setColUpper(col, std::floor(origColSolution[col]));
      solveFromHotStart(solver, col, true, origUB, std::floor(origColSolution[col]));
      if (solver->isProvenPrimalInfeasible()) {
        downBranchFeasible = false;
        downBranchObjs.push_back(1e50);
      } else {
        downBranchObjs.push_back(solver->getObjValue());
      }

      // Return to original state
      solver->setColUpper(col, origUB);
      solver->solveFromHotStart();

      // Check up branch
      solver->setColLower(col, std::ceil(origColSolution[col]));
      solveFromHotStart(solver, col, false, origLB, std::ceil(origColSolution[col]));
      if (solver->isProvenPrimalInfeasible()) {
        upBranchFeasible = false;
        upBranchObjs.push_back(1e50);
      } else {
        upBranchObjs.push_back(solver->getObjValue());
      }

      // Return to original state
      solver->setColLower(col, origLB);
      solver->solveFromHotStart();

      // calculate score - minimum of down and up branches
      scores.push_back(metric.calculateValue(2, col, downBranchObjs.back(),
                                             upBranchObjs.back(), 0));

      // Make sure one side of the split is feasible
      if (!downBranchFeasible && !upBranchFeasible) {
        // Infeasible problem
        error_msg(errorstring, "Infeasible problem due to integer variable %d.\n",
                  col, origColSolution[col]);
        writeErrorToLog(errorstring, params.logfile);
        exit(1);
      }
    } // check if integer
  } // end iterating over columns

  if (fracIntVars.size() == 0) {
    if (solver)
      delete solver;
    return DisjExitReason::NO_DISJUNCTION_EXIT;
  }

  // Sort by scores
  std::vector<int> sortedIndices(fracIntVars.size());
  for (int i = 0; i < fracIntVars.size(); i++) {
    sortedIndices[i] = i;
  }
  std::sort(sortedIndices.begin(), sortedIndices.end(),
            [&](const int i, const int j)
            { return scores[i] > scores[j]; });

  // Select the top numVarsToSelect variables
  int numFracIntVars = fracIntVars.size();
  std::vector<int> selectedIndices;
  for (int i = 0; i < std::min(numFracIntVars, numVarsToSelect); i++) {
    selectedIndices.push_back(fracIntVars[sortedIndices[i]]);
  }

  // get the branching information for the selected variables
  std::vector<std::vector<int> > waysByVar;
  std::vector<std::vector<double> > valuesByVar;
  for (int i : selectedIndices) {
    std::vector<int> ways;
    std::vector<double> vals;
    // down branch is 1 because it's upper bounded
    // (i.e. a rhs bound, like to the right on the number line)
    ways.push_back(1);
    vals.push_back(std::floor(origColSolution[i]));
    // analogously, up branch is 0
    ways.push_back(0);
    vals.push_back(std::ceil(origColSolution[i]));
    waysByVar.push_back(ways);
    valuesByVar.push_back(vals);
  }

  // take cartesian product of ways and values by variable to get them listed by disjunctive term
  // recall, the goal here is to just make a partial BB tree from the cartesian product of best branching decisions
  std::vector<int> tempOutput;
  std::vector<double> tempOutput2;
  std::vector<std::vector<int> > waysByTerm;
  std::vector<std::vector<double> > valuesByTerm;
  cartesianProduct(waysByTerm, tempOutput, waysByVar.begin(), waysByVar.end());
  cartesianProduct(valuesByTerm, tempOutput2, valuesByVar.begin(), valuesByVar.end());

  // add each term to the disjunction
  int actualNumTerms = std::pow(2, selectedIndices.size());
  for (int i = 0; i < actualNumTerms; i++) {
    addTerm(selectedIndices, waysByTerm[i], valuesByTerm[i], solver);
  }

  // turn off hot start and factorization
  solver->unmarkHotStart();
  solver->disableFactorization();
  if (solver)
    delete solver;

  // check to make sure we have a feasible term
  for (int i = 0; i < actualNumTerms; i++) {
    if (terms[i].feasible){
      return DisjExitReason::SUCCESS_EXIT;
    }
  }
  return DisjExitReason::NO_DISJUNCTION_EXIT;;
} /* prepareDisjunction */

/****************** PROTECTED **********************/
void StrongBranchingDisjunction::initialize(const StrongBranchingDisjunction* const source,
    const VPCParametersNamespace::VPCParameters* const params) {
  VPCDisjunction::initialize(source, params);
  if (source) {
    // copy the source
  } else {
    // initialize to default values
  }
} /* initialize */

void StrongBranchingDisjunction::setCgsName(const int var, const double val) {
  return;
} /* setCgsName */

void StrongBranchingDisjunction::addTerm(const std::vector<int>& branching_variables,
    const std::vector<int>& branching_ways, const std::vector<double>& branching_values,
    OsiSolverInterface* solver) {

  std::vector<double> oldBounds(branching_variables.size());
  std::string disjTermName = "";
  int valxk;

  // set the bounds for this term on the solver and resolve
  for (int i = 0; i < branching_variables.size(); i++) {
    disjTermName += "x" + std::to_string(branching_variables[i]);
    if (branching_ways[i]) {
      // down branch
      oldBounds[i] = solver->getColUpper()[branching_variables[i]];
      solver->setColUpper(branching_variables[i], branching_values[i]);
      disjTermName += " <= ";
    } else {
      // up branch
      oldBounds[i] = solver->getColLower()[branching_variables[i]];
      solver->setColLower(branching_variables[i], branching_values[i]);
      disjTermName += " >= ";
    }
    valxk = static_cast<int>(branching_values[i]);
    disjTermName += std::to_string(valxk);
    disjTermName += i < branching_variables.size() - 1 ? " ^ " : "";
  }
  solver->solveFromHotStart();

  // create the disjunctive term object and add it to the disjunction
  DisjunctiveTerm term;
  term.basis = dynamic_cast<CoinWarmStartBasis*>(solver->getWarmStart());
  term.obj = solver->getObjValue();
  term.changed_var = branching_variables;
  term.changed_bound = branching_ways;
  term.changed_value = branching_values;
  term.feasible = solver->isProvenPrimalInfeasible() ? false : true;
  this->terms.push_back(term);
  this->num_terms++;
  Disjunction::setCgsName(this->name, disjTermName);

  // reset the solver to its original bounds
  for (int i = 0; i < branching_variables.size(); i++) {
    if (branching_ways[i]) {
      // down branch
      solver->setColUpper(branching_variables[i], oldBounds[i]);
    } else {
      // up branch
      solver->setColLower(branching_variables[i], oldBounds[i]);
    }
  }
  solver->solveFromHotStart();

} /* addTerm */

// cartesian product of vector of vectors
template <typename T>
void StrongBranchingDisjunction::cartesianProduct(
    /// [out] empty vector of vectors that will return as cartesian product
    std::vector< std::vector<T> >& finalResult,
    /// [out] empty vector that will have work done in it
    std::vector<T>& currentResult,
    /// [in] begin iterator for input vector of vectors to cartesian product
    typename std::vector< std::vector<T> >::const_iterator currentInput,
    /// [in] end iterator for input vector of vectors to cartesian product
    typename std::vector< std::vector<T> >::const_iterator lastInput) {

  if(currentInput == lastInput) {
    // terminal condition of the recursion. We no longer have
    // any input vectors to manipulate. Add the current result (currentResult)
    // to the total set of results (finalResult).
    finalResult.push_back(currentResult);
    return;
  }

  // need an easy name for my vector-of-T
  const std::vector<T>& temp = *currentInput;
  for(typename std::vector<T>::const_iterator it = temp.begin(); it != temp.end(); it++) {
    currentResult.push_back(*it);  // add currentInput
    cartesianProduct(finalResult, currentResult, currentInput+1, lastInput);
    currentResult.pop_back(); // clean currentInput off for next round
  }
}
