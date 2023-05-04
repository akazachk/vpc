/**
 * @file CompleteDisjunction.cpp
 * @author Sean Kelley
 * @date 2023-04-20
 */
#include "CompleteDisjunction.hpp"

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
/** Default constructor */
CompleteDisjunction::CompleteDisjunction() {
  initialize(NULL);
} /* default constructor */

/** Copy constructor */
CompleteDisjunction::CompleteDisjunction(const CompleteDisjunction& source) {
  initialize(&source);
} /* copy constructor */

/** Destructor */
CompleteDisjunction::~CompleteDisjunction() {
} /* destructor */

/** Assignment operator */
CompleteDisjunction& CompleteDisjunction::operator=(const CompleteDisjunction& source) {
  if (this != &source) {
    initialize(&source);
  }
  return *this;
} /* assignment operator */

/** Clone */
CompleteDisjunction* CompleteDisjunction::clone() const {
  return new CompleteDisjunction(*this);
} /* clone */

/**
 * Set up the disjunction class as new
 */
void CompleteDisjunction::setupAsNew() {
  Disjunction::setupAsNew();
  partialDisj = NULL;
} /* setupAsNew */

/**
 * @brief Get disjunction
 *
 * This will throw away all the information from the old disjunction (if one exists)
 */
DisjExitReason CompleteDisjunction::prepareDisjunction(const OsiSolverInterface* const si,
                                                       const Disjunction* const partialDisj) {
  setupAsNew();
  this->partialDisj = partialDisj->clone();
  prepareDisjunction(si);
}

/****************** PROTECTED **********************/
void CompleteDisjunction::initialize(const CompleteDisjunction* const source) {
  Disjunction::initialize(source);
  if (source) {
    partialDisj = source->partialDisj;
  } else {
    partialDisj = NULL;
  }
} /* initialize */

/**
 * @brief Prepare a new disjunction
 *
 * This will throw away all the information from the old disjunction, except it will not reset the timer
 */
DisjExitReason CompleteDisjunction::prepareDisjunction(const OsiSolverInterface* const si) {
  return DisjExitReason::NO_DISJUNCTION_EXIT;
} /* prepareDisjunction */

void CompleteDisjunction::addTerm(const std::vector<int>& branching_variables,
                                  const std::vector<int>& branching_ways,
                                  const std::vector<double>& branching_values,
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
void CompleteDisjunction::cartesianProduct(
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