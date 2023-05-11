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
#include <vector> // vector

// Project files
#include "CglVPC.hpp" // get timer information too
#include "SolverInterface.hpp" // SolverInterface
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
  VPCDisjunction::setupAsNew();
  partialDisj = NULL;
} /* setupAsNew */

/**
 * @brief Get disjunction
 *
 * This will throw away all the information from the old disjunction (if one exists)
 */
DisjExitReason CompleteDisjunction::prepareDisjunction(const OsiSolverInterface* const si,
                                                       const VPCDisjunction* const partialDisj) {
  setupAsNew();
  VPCDisjunction::initialize(NULL, &partialDisj->params);
  this->partialDisj = partialDisj->clone();
  return prepareDisjunction(si);
}

/****************** PROTECTED **********************/
void CompleteDisjunction::initialize(const CompleteDisjunction* const source) {
  if (source) {
    VPCDisjunction::initialize(source, &source->params);
    partialDisj = source->partialDisj;
  } else {
    VPCDisjunction::initialize(NULL, NULL);
    partialDisj = NULL;
  }
} /* initialize */

/**
 * @brief Prepare a new disjunction
 *
 * This will throw away all the information from the old disjunction, except it will not reset the timer
 */
DisjExitReason CompleteDisjunction::prepareDisjunction(const OsiSolverInterface* const si) {

  // get the solver
  std::shared_ptr<SolverInterface> solver = getSolver(si);

  // expand the tree for the common fixed terms
  addCommonFixedTerms(solver.get());

  // expand the tree for the fixes in each disjunction
  addDisjunctionFixedTerms(solver.get());

  // remove terms that contain other terms (i.e. are not leaf nodes)
  removeRedundantTerms(si);

  // add objective information
  best_obj = partialDisj->best_obj;
  worst_obj = partialDisj->worst_obj;
  integer_obj = partialDisj->integer_obj;
  integer_sol = partialDisj->integer_sol;

  // set the name
  for (int i = 0; i < num_terms; i++) {
    Disjunction::setCgsName(this->name, this->terms[i].name);
  }

  // clean up the solver before it's removal
  solver->unmarkHotStart();
  solver->disableFactorization();

  return DisjExitReason::SUCCESS_EXIT;
} /* prepareDisjunction */

/// @brief add the terms arising from the tree that created the fixed variables
/// common for all disjunctive terms
void CompleteDisjunction::addCommonFixedTerms(OsiSolverInterface* solver){
  std::vector<std::vector<int> > variablesByTerm;
  std::vector<std::vector<int> > boundsByTerm;
  std::vector<std::vector<double> > valuesByTerm;

  // add the terms for the common fixed variables
  getTerms(variablesByTerm, boundsByTerm, valuesByTerm, partialDisj->common_changed_var,
           partialDisj->common_changed_bound, partialDisj->common_changed_value);
  // the first term will be further broken up by terms in the partial disjunction
  for (int i = 1; i < boundsByTerm.size(); i++) {
    addTerm(variablesByTerm[i], boundsByTerm[i], valuesByTerm[i], solver);
  }
}

/// @brief add the terms arising from the tree that created fixed variables for each term
void CompleteDisjunction::addDisjunctionFixedTerms(OsiSolverInterface* solver){

  std::vector<std::vector<int> > variablesByTerm;
  std::vector<std::vector<int> > boundsByTerm;
  std::vector<std::vector<double> > valuesByTerm;

  for (int i = 0; i < partialDisj->num_terms; i++) {

    // get the info for the terms
    getTerms(variablesByTerm, boundsByTerm, valuesByTerm, partialDisj->terms[i].changed_var,
             partialDisj->terms[i].changed_bound, partialDisj->terms[i].changed_value,
             &partialDisj->common_changed_var, &partialDisj->common_changed_bound,
             &partialDisj->common_changed_value);

    // create each term
    for (int j = 0; j < boundsByTerm.size(); j++) {
      addTerm(variablesByTerm[j], boundsByTerm[j], valuesByTerm[j], solver);
    }

    // reset the info for the next term
    variablesByTerm.clear();
    boundsByTerm.clear();
    valuesByTerm.clear();
  }
}

/// @brief Complete the tree encoded in the disjunction
void CompleteDisjunction::getTerms(
    /// [out] the variables to fix for each term
    std::vector<std::vector<int> >& variablesByTerm,
    /// [out] the bounds to fix for each term
    std::vector<std::vector<int> >& boundsByTerm,
    /// [out] the values to fix respective bounds for each term
    std::vector<std::vector<double> >& valuesByTerm,
    /// [in] which variable had a bound changed
    const std::vector<int>& changedVariables,
    /// [in] which bound was changed (0 - lower; 1 - upper)
    const std::vector<int>& changedBounds,
    /// [in] which value was the bound changed to
    const std::vector<double>& changedValues,
    /// [in] which variables are fixed for these terms
    const std::vector<int>* fixedVariables,
    /// [in] which bounds are fixed for these terms
    const std::vector<int>* fixedBounds,
    /// [in] which values are the fixed bounds for these terms set to
    const std::vector<double>* fixedValues){

  // validate inputs
  if ((fixedBounds == NULL) != (fixedValues == NULL) ||
      (fixedValues == NULL) != (fixedVariables == NULL)) {
    error_msg(errorstring, "fixedValues, fixedBounds, and fixedVariables must"
                           "either be all null or all nonnull.\n");
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }
  if (changedVariables.size() != changedBounds.size() ||
      changedBounds.size() != changedValues.size()) {
    error_msg(errorstring, "changedVariables, changedBounds, and changedValues"
                           "must be the same size.\n");
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }
  if ((fixedBounds != NULL) && (fixedValues != NULL) && (fixedBounds != NULL) &&
      ((fixedBounds->size() != fixedValues->size()) || (fixedValues->size() != fixedVariables->size()))) {
    error_msg(errorstring, "fixedBounds, fixedValues, and fixedVariables must"
                           "be the same size.\n");
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }

  // if no inputs, just return
  if ((changedVariables.size() == 0) && (changedValues.size() == 0) &&
      (changedBounds.size() == 0)) {
    return;
  }

  // containers for creating disjunctive terms
  std::vector<int> workingVariables = changedVariables;
  std::vector<int> workingBounds = changedBounds;
  std::vector<double> workingValues = changedValues;
  int newBound;
  double newVal;

  // add the existing term first
  variablesByTerm.push_back(workingVariables);
  boundsByTerm.push_back(workingBounds);
  valuesByTerm.push_back(workingValues);

  // then negate and remove branching decisions one at a time
  do {
    // swap the last variable's branching direction
    newBound = workingBounds.back() == 0 ? 1 : 0;
    // if newBound is upper bound, get the next value below the current bound and vice versa
    newVal = workingValues.back() + (newBound == 1 ? -1 : 1);
    workingBounds.pop_back();
    workingValues.pop_back();
    workingBounds.push_back(newBound);
    workingValues.push_back(newVal);

    // add the term
    variablesByTerm.push_back(workingVariables);
    boundsByTerm.push_back(workingBounds);
    valuesByTerm.push_back(workingValues);

    // pop the last variable, moving a level up in the tree
    workingVariables.pop_back();
    workingBounds.pop_back();
    workingValues.pop_back();

  } while (workingVariables.size() > 0);

  // prepend the fixed bounds and values to each term
  if ((fixedBounds != NULL) && (fixedValues != NULL) && (fixedVariables != NULL)) {
    for (int i = 0; i < boundsByTerm.size(); i++) {
      variablesByTerm[i].insert(variablesByTerm[i].begin(), fixedVariables->begin(), fixedVariables->end());
      boundsByTerm[i].insert(boundsByTerm[i].begin(), fixedBounds->begin(), fixedBounds->end());
      valuesByTerm[i].insert(valuesByTerm[i].begin(), fixedValues->begin(), fixedValues->end());
    }
  }
} /* getTerms */

/// @brief Add disjunctive term - split this into creating the term and adding it
/// to the disjunction. Then inbetween check if the term contains another term
void CompleteDisjunction::addTerm(const std::vector<int>& branching_variables,
                                  const std::vector<int>& branching_bounds,
                                  const std::vector<double>& branching_values,
                                  OsiSolverInterface* solver) {

  std::vector<double> oldBounds(branching_variables.size());
  std::string disjTermName = "";
  int valxk;

  // set the bounds for this term on the solver and resolve
  for (int i = 0; i < branching_variables.size(); i++) {
    disjTermName += "x" + std::to_string(branching_variables[i]);
    if (branching_bounds[i]) {
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
    disjTermName += i < branching_variables.size() - 1 ? "; " : "";
  }
  solver->solveFromHotStart();

  // create the disjunctive term object and add it to the disjunction
  DisjunctiveTerm term;
  term.basis = dynamic_cast<CoinWarmStartBasis*>(solver->getWarmStart());
  term.obj = solver->getObjValue();
  term.changed_var = branching_variables;
  term.changed_bound = branching_bounds;
  term.changed_value = branching_values;
  term.feasible = solver->isProvenPrimalInfeasible() ? false : true;
  term.name = disjTermName;
  this->terms.push_back(term);
  this->num_terms++;

  // reset the solver to its original bounds
  for (int i = 0; i < branching_variables.size(); i++) {
    if (branching_bounds[i]) {
      // down branch
      solver->setColUpper(branching_variables[i], oldBounds[i]);
    } else {
      // up branch
      solver->setColLower(branching_variables[i], oldBounds[i]);
    }
  }
  solver->solveFromHotStart();
} /* addTerm */

/// @brief finds and removes redundant terms from the disjunction
/// (i.e. terms with variable bounds that contain another term's variable bounds)
void CompleteDisjunction::removeRedundantTerms(const OsiSolverInterface* const si){

  std::vector<std::vector<int> > containsTerms(num_terms);
  std::vector<std::vector<int> > equalsTerms(num_terms);

  // get each term's solver
  std::vector<OsiSolverInterface*> solver(num_terms);
  for (int i = 0; i < num_terms; i++) {
    getSolverForTerm(solver[i], i, si, true, .001, params.logfile, false, false, false);
  }

  // check for containment
  bool i_contains_j, j_contains_i;
  for (int i = 0; i < num_terms; i++) {
    for (int j = i + 1; j < num_terms; j++) {
      i_contains_j = variableBoundsContained(solver[i], solver[j]);
      j_contains_i = variableBoundsContained(solver[j], solver[i]);
      if (i_contains_j && j_contains_i) {
        equalsTerms[j].push_back(i); // use the larger index so we can remove in reverse order
      } else if (i_contains_j) {
        containsTerms[i].push_back(j); // term i contains term j
      } else if (j_contains_i) {
        containsTerms[j].push_back(i); // term j contains term i
      }
    }
  }

  // remove all redundant terms
  for (int i = num_terms - 1; i >= 0; i--) {
    if (containsTerms[i].size() > 0 || equalsTerms[i].size() > 0) {
      this->terms.erase(this->terms.begin() + i);
      this->num_terms--;
    }
  }
}
