// Name:     SplitDisjunction.cpp
// Author:   A. M. Kazachkov
// Date:     2018-02-22
//-----------------------------------------------------------------------------
#include "SplitDisjunction.hpp"

// COIN-OR

// Project files
#include "BBHelper.hpp"
#include "CglVPC.hpp" // get timer information too
#include "SolverHelper.hpp"
#include "utility.hpp"
#include "VPCEventHandler.hpp"

#ifdef TRACE
#include "debug.hpp"
#endif

#ifdef USE_CLP
#include <OsiClpSolverInterface.hpp>
#endif

/****************** PUBLIC  **********************/
/** Handle parameters */
SplitDisjunction::SplitDisjunction(const VPCParameters& param) {
  initialize(NULL, &param);
} /* param constructor */
void SplitDisjunction::setParams(const VPCParameters& param) {
  this->params = param;
} /* setParams */

/** Default constructor */
SplitDisjunction::SplitDisjunction() {
  initialize();
} /* default constructor */

/** Copy constructor */
SplitDisjunction::SplitDisjunction(const SplitDisjunction& source) : Disjunction(source) {
  initialize(&source);
} /* copy constructor */

/** Destructor */
SplitDisjunction::~SplitDisjunction() {
} /* destructor */

/** Assignment operator */
SplitDisjunction& SplitDisjunction::operator=(const SplitDisjunction& source) {
  if (this != &source) {
    initialize(&source);
  }
  return *this;
} /* assignment operator */

/** Clone */
SplitDisjunction* SplitDisjunction::clone() const {
  return new SplitDisjunction(*this);
} /* clone */

/** Set up the disjunction class as new (except the timer pointer, and do not reset params) */
void SplitDisjunction::setupAsNew() {
  Disjunction::setupAsNew();
} /* setupAsNew */

/**
 * @brief Prepare a new disjunction
 *
 * This will throw away all the information from the old disjunction, except it will not reset the timer
 */
ExitReason SplitDisjunction::prepareDisjunction(OsiSolverInterface* const si) {
  // Reset things in case we are reusing the class for some reason
  setupAsNew();
  if (!timer) {
    error_msg(errorstring, "Timer is not set.\n");
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }

  if (timer)
    timer->start_timer(CglVPC::VPCTimeStatsName[static_cast<int>(CglVPC::VPCTimeStats::DISJ_GEN_TIME)]);

  OsiVectorInt fracCore = si->getFractionalIndices(params.get(doubleConst::AWAY));

  // Set up solver for hot start
  SolverInterface* solver;
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
  solver->enableFactorization();
  solver->markHotStart();
  for (int var : fracCore) {
    if (num_terms >= params.get(DISJ_TERMS))
      break;

    const double val = solver->getColSolution()[var];
    const double floorxk = std::floor(val);
    const double ceilxk = std::ceil(val);

    if (!si->isInteger(var)) {
      error_msg(errorstring, "Chosen variable %d is not an integer variable.\n", var);
      writeErrorToLog(errorstring, params.logfile);
      exit(1);
    }
    if (isVal(val, floorxk, params.get(doubleConst::AWAY))
          || isVal(val, ceilxk, params.get(doubleConst::AWAY))) {
      error_msg(errorstring, "Chosen variable %d is not fractional (value: %1.6e).\n", var, val);
      writeErrorToLog(errorstring, params.logfile);
      exit(1);
    }

    const double origLB = solver->getColLower()[var];
    const double origUB = solver->getColUpper()[var];
    bool downBranchFeasible = true, upBranchFeasible = true;

    // Check down branch
    solver->setColUpper(var, floorxk);
    solveFromHotStart(solver, var, true, origUB, floorxk);
    if (solver->isProvenOptimal()) {
      addTerm(var, 1, floorxk, solver);
    } else if (solver->isProvenPrimalInfeasible()) {
      downBranchFeasible = false;
    } else {
      // Something strange happened
      error_msg(errorstring,
          "Down branch is neither optimal nor primal infeasible on variable %d (value %e).\n",
          var, val);
      writeErrorToLog(errorstring, params.logfile);
      exit(1);
    }
    solver->setColUpper(var, origUB);

    // Return to previous state
    solver->solveFromHotStart();

    // Check up branch
    solver->setColLower(var, ceilxk);
    solveFromHotStart(solver, var, false, origLB, ceilxk);
    if (solver->isProvenOptimal()) {
      addTerm(var, 0, ceilxk, solver);
    } else if (solver->isProvenPrimalInfeasible()) {
      upBranchFeasible = false;
    } else {
      // Something strange happened
      error_msg(errorstring,
          "Up branch is neither optimal nor primal infeasible on variable %d (value %e).\n",
          var, val);
      writeErrorToLog(errorstring, params.logfile);
      exit(1);
    }
    solver->setColLower(var, origLB);

    // Return to original state
    solver->solveFromHotStart();

    // Check if some side of the split is infeasible
    if (!downBranchFeasible || !upBranchFeasible) {
      if (!downBranchFeasible && !upBranchFeasible) {
        // Infeasible problem
        error_msg(errorstring,
            "Infeasible problem due to integer variable %d (value %e).\n",
            var, val);
        writeErrorToLog(errorstring, params.logfile);
        exit(1);
      }
      // If one side infeasible, can delete last term added and fix variables instead 
      this->terms.resize(this->num_terms - 1);
      this->num_terms--;
      this->common_changed_var.push_back(var);
      if (!downBranchFeasible) { // set lb to ceilxk
        this->common_changed_bound.push_back(0);
        this->common_changed_value.push_back(ceilxk);
      }
      if (!upBranchFeasible) { // set ub to floorxk
        this->common_changed_bound.push_back(1);
        this->common_changed_value.push_back(floorxk);
      }
    }
  } // loop through fractional core
  solver->unmarkHotStart();
  solver->disableFactorization();
  if (solver)
    delete solver;

  if (timer)
    timer->end_timer(CglVPC::VPCTimeStatsName[static_cast<int>(CglVPC::VPCTimeStats::DISJ_GEN_TIME)]);
  return ExitReason::SUCCESS_EXIT;
} /* prepareDisjunction */

/****************** PROTECTED **********************/
void SplitDisjunction::initialize(const SplitDisjunction* const source,
    const VPCParameters* const params) {
  if (params != NULL) {
    setParams(*params);
  }
  if (source) {
    Disjunction::initialize(source);
    if (!params) {
      setParams(source->params);
    }
  } else {
    setupAsNew();
  }
} /* initialize */

void SplitDisjunction::addTerm(const int branching_variable, 
    const int branching_way, const double branching_value, 
    const OsiSolverInterface* const solver) {
  DisjunctiveTerm term;
  term.basis = dynamic_cast<CoinWarmStartBasis*>(solver->getWarmStart());
  term.obj = solver->getObjValue();
  term.changed_var.resize(1, branching_variable); 
  term.changed_bound.resize(1, branching_way); 
  term.changed_value.resize(1, branching_value); 
  this->terms.push_back(term);
  this->num_terms++;
} /* addTerm */
