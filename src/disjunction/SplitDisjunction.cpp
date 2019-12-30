/**
 * @file SplitDisjunction.cpp
 * @author A. M. Kazachkov
 * @date 2018-02-22
 */
#include "SplitDisjunction.hpp"

// COIN-OR

// Project files
#include "CglVPC.hpp" // get timer information too
#include "SolverHelper.hpp"
#include "utility.hpp"
#include "VPCEventHandler.hpp"

#ifdef TRACE
#include "vpc_debug.hpp"
#endif

#ifdef USE_CLP
#include <OsiClpSolverInterface.hpp>
#endif

/****************** PUBLIC  **********************/
/** Param constructor */
SplitDisjunction::SplitDisjunction(const VPCParametersNamespace::VPCParameters& param) {
  initialize(NULL, &param);
} /* param constructor */

/** Copy and param constructor */
SplitDisjunction::SplitDisjunction(const SplitDisjunction& source, const VPCParametersNamespace::VPCParameters& param) {
  initialize(&source, &param);
} /* copy & param constructor */

/** Default constructor */
SplitDisjunction::SplitDisjunction() {
  initialize(NULL, NULL);
} /* default constructor */

/** Copy constructor */
SplitDisjunction::SplitDisjunction(const SplitDisjunction& source) {
  initialize(&source, NULL);
} /* copy constructor */

/** Destructor */
SplitDisjunction::~SplitDisjunction() {
} /* destructor */

/** Assignment operator */
SplitDisjunction& SplitDisjunction::operator=(const SplitDisjunction& source) {
  if (this != &source) {
    initialize(&source, NULL);
  }
  return *this;
} /* assignment operator */

/** Clone */
SplitDisjunction* SplitDisjunction::clone() const {
  return new SplitDisjunction(*this);
} /* clone */

/**
 * Set up the disjunction class as new
 * (except the timer pointer, and do not reset params;
 * in addition, here we do not reset var, in case user has set it)
 */
void SplitDisjunction::setupAsNew() {
  VPCDisjunction::setupAsNew();
} /* setupAsNew */

/**
 * @brief Prepare a new disjunction
 *
 * This will throw away all the information from the old disjunction, except it will not reset the timer
 */
DisjExitReason SplitDisjunction::prepareDisjunction(const OsiSolverInterface* const si) {
  // Reset things in case we are reusing the class for some reason
  setupAsNew();
//  if (!timer) {
//    error_msg(errorstring, "Timer is not set.\n");
//    writeErrorToLog(errorstring, params.logfile);
//    exit(1);
//  }

  if (timer)
    timer->start_timer(CglVPC::VPCTimeStatsName[static_cast<int>(CglVPC::VPCTimeStats::DISJ_GEN_TIME)]);

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

  // If var >= 0, that means user has picked a variable already
  // Otherwise, var < 0, and we need to pick one
  DisjExitReason retVal = DisjExitReason::UNKNOWN;
  if (var >= 0) {
    if (checkVar(solver, var)) {
      retVal = DisjExitReason::SUCCESS_EXIT;
    } else {
      retVal = DisjExitReason::NO_DISJUNCTION_EXIT;
    }
  }
  else { // loop through variables and find a split disjunction we can use
    std::vector<int> fracCore;
    std::vector<double> fractionality;
    for (int col = 0; col < solver->getNumCols(); col++) {
      if (!solver->isInteger(col))
        continue;
      const double val = solver->getColSolution()[col];
      const double floorxk = std::floor(val);
      const double ceilxk = std::ceil(val);
      const double frac = CoinMin(val - floorxk, ceilxk - val);
      if (!isVal(frac, 0., params.get(VPCParametersNamespace::doubleConst::AWAY))) {
        fracCore.push_back(col);
        fractionality.push_back(frac);
      }
    }
    
    if (fracCore.size() == 0) {
      if (solver)
        delete solver;
      return DisjExitReason::NO_DISJUNCTION_EXIT;
    }

    // Sort by fractionality
    std::vector<unsigned> sortIndex(fracCore.size());
    for (unsigned i = 0; i < fracCore.size(); i++) {
      sortIndex[i] = i;
    }
    std::sort(sortIndex.begin(), sortIndex.end(),
        [&](const int i, const int j)
        { return fractionality[i] > fractionality[j]; }); 
        //index_cmp_dsc<const std::vector<double>&>(fractionality)); // descending

    for (unsigned i : sortIndex) {
      int col = fracCore[i];
      if (checkVar(solver, col)) {
#ifdef TRACE
        printf("Selected var %d (value %1.3f) for split variable to branch on.\n", col, solver->getColSolution()[col]);
#endif
        retVal = DisjExitReason::SUCCESS_EXIT;
        break; // stop when we find a disjunction; might want to sort by strong branching value instead, but that can be too expensive
      }
    } // loop through fractional core
  } // loop through variables and find a split disjunction we can use
  solver->unmarkHotStart();
  solver->disableFactorization();
  if (solver)
    delete solver;

  if (timer)
    timer->end_timer(CglVPC::VPCTimeStatsName[static_cast<int>(CglVPC::VPCTimeStats::DISJ_GEN_TIME)]);
  if (retVal == DisjExitReason::SUCCESS_EXIT) {
    setCgsName(var, si->getColSolution()[var]);
  }
  return retVal;
} /* prepareDisjunction */

/****************** PROTECTED **********************/
void SplitDisjunction::initialize(const SplitDisjunction* const source,
    const VPCParametersNamespace::VPCParameters* const params) {
  VPCDisjunction::initialize(source, params);
  if (source) {
    this->var = source->var;
  } else {
    this->var = -1;
  }
} /* initialize */

/**
 * @brief Checks whether we can add a term from the given variable
 *
 * It is assumed that the solver is already set up for hot start
 */
bool SplitDisjunction::checkVar(OsiSolverInterface* si, int col) {
  SolverInterface* solver = dynamic_cast<SolverInterface*>(si);
  const double val = solver->getColSolution()[col];
  const double floorxk = std::floor(val);
  const double ceilxk = std::ceil(val);

  if (!si->isInteger(col)) {
    error_msg(errorstring, "Chosen variable %d is not an integer variable.\n", col);
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }
  if (isVal(val, floorxk, params.get(VPCParametersNamespace::doubleConst::AWAY))
        || isVal(val, ceilxk, params.get(VPCParametersNamespace::doubleConst::AWAY))) {
    error_msg(errorstring, "Chosen variable %d is not fractional (value: %1.6e).\n", col, val);
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }

  const double origLB = solver->getColLower()[col];
  const double origUB = solver->getColUpper()[col];
  bool downBranchFeasible = true, upBranchFeasible = true;

  // Check down branch
  solver->setColUpper(col, floorxk);
  solveFromHotStart(solver, col, true, origUB, floorxk);
  if (solver->isProvenOptimal()) {
    addTerm(col, 1, floorxk, solver);
  } else if (solver->isProvenPrimalInfeasible()) {
    downBranchFeasible = false;
  } else {
    // Something strange happened
    error_msg(errorstring,
        "Down branch is neither optimal nor primal infeasible on variable %d (value %e).\n",
        col, val);
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }
  solver->setColUpper(col, origUB);

  // Return to previous state
  solver->solveFromHotStart();

  // Check up branch
  solver->setColLower(col, ceilxk);
  solveFromHotStart(solver, col, false, origLB, ceilxk);
  if (solver->isProvenOptimal()) {
    addTerm(col, 0, ceilxk, solver);
  } else if (solver->isProvenPrimalInfeasible()) {
    upBranchFeasible = false;
  } else {
    // Something strange happened
    error_msg(errorstring,
        "Up branch is neither optimal nor primal infeasible on variable %d (value %e).\n",
        col, val);
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
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
          col, val);
      writeErrorToLog(errorstring, params.logfile);
      exit(1);
    }
    // If one side infeasible, can delete last term added and fix variables instead
    this->terms.resize(this->num_terms - 1);
    this->num_terms--;
    this->common_changed_var.push_back(col);
    if (!downBranchFeasible) { // set lb to ceilxk
      this->common_changed_bound.push_back(0);
      this->common_changed_value.push_back(ceilxk);
    }
    if (!upBranchFeasible) { // set ub to floorxk
      this->common_changed_bound.push_back(1);
      this->common_changed_value.push_back(floorxk);
    }
    return false;
  } else {
    this->var = col;
    return true;
  }
} /* checkVar */

void SplitDisjunction::setCgsName(const int var, const double val) {
  if (var < 0)
    return;
  const int floorxk = static_cast<int>(std::floor(val));
  const int ceilxk = static_cast<int>(std::ceil(val));
  std::string disjTermName = "";
  disjTermName += "x" + std::to_string(var);
  disjTermName += " <= ";
  disjTermName += std::to_string(floorxk);
  Disjunction::setCgsName(this->name, disjTermName);

  disjTermName = "";
//  disjTermName += ") V (";
  disjTermName += "x" + std::to_string(var);
  disjTermName += " >= ";
  disjTermName += std::to_string(ceilxk);
//  disjTermName += ")";
  Disjunction::setCgsName(this->name, disjTermName);
} /* setCgsName */

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
