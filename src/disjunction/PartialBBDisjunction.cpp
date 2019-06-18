// Name:     PartialBBDisjunction.cpp
// Author:   A. M. Kazachkov
// Date:     2018-02-22
//-----------------------------------------------------------------------------
#include "PartialBBDisjunction.hpp"

// COIN-OR
#ifdef USE_CBC
#include <CbcModel.hpp>

// Cbc
#include <CbcSolver.hpp>
#include <CbcTree.hpp>

// Variable selection
#include <CbcBranchDefaultDecision.hpp>
#include <CbcBranchStrongDecision.hpp>
#include <OsiChooseStrongCustom.hpp>
#include <CbcBranchDynamic.hpp>

// Node selection
#include <CbcCompareDefault.hpp>
#include <CbcCompareBFS.hpp>
#include <CbcCompareDepth.hpp>
#include <CbcCompareEstimate.hpp>
#include <CbcCompareObjective.hpp>

// General strategy
#include <CbcStrategy.hpp>
#endif // USE_CBC

// Project files
#include "CglVPC.hpp" // get timer information too
#include "SolverHelper.hpp"
#include "utility.hpp"
#include "VPCEventHandler.hpp" // currently requires Cbc

#ifdef TRACE
#include "debug.hpp"
#endif

/****************** PUBLIC  **********************/
/** Handle parameters */
PartialBBDisjunction::PartialBBDisjunction(const VPCParameters& param) {
  initialize(NULL, &param);
} /* param constructor */
void PartialBBDisjunction::setParams(const VPCParameters& param) {
  this->params = param;
} /* setParams */

/** Default constructor */
PartialBBDisjunction::PartialBBDisjunction() {
  initialize();
} /* default constructor */

/** Copy constructor */
PartialBBDisjunction::PartialBBDisjunction(const PartialBBDisjunction& source) : Disjunction(source) {
  initialize(&source);
} /* copy constructor */

/** Destructor */
PartialBBDisjunction::~PartialBBDisjunction() {
} /* destructor */

/** Assignment operator */
PartialBBDisjunction& PartialBBDisjunction::operator=(const PartialBBDisjunction& source) {
  if (this != &source) {
    initialize(&source);
  }
  return *this;
} /* assignment operator */

/** Clone */
PartialBBDisjunction* PartialBBDisjunction::clone() const {
  return new PartialBBDisjunction(*this);
} /* clone */

/** Set up the disjunction class as new (except the timer pointer, and do not reset params) */
void PartialBBDisjunction::setupAsNew() {
  Disjunction::setupAsNew();
  this->data.num_nodes_on_tree = 0;
  this->data.num_partial_bb_nodes = 0;
  this->data.num_pruned_nodes = 0;
  this->data.min_node_depth = std::numeric_limits<int>::max();
  this->data.max_node_depth = 0;
  this->data.num_fixed_vars = 0;
  //    this->data.stats.resize(0);
  //    this->data.pruned_stats.resize(0);
  //    this->data.node_id.resize(0);
} /* setupAsNew */

/**
 * @brief Prepare a new disjunction
 *
 * This will throw away all the information from the old disjunction, except it will not reset the timer
 */
ExitReason PartialBBDisjunction::prepareDisjunction(const OsiSolverInterface* const si) {
  // Reset things in case we are reusing the class for some reason
  setupAsNew();
//  if (!timer) {
//    error_msg(errorstring, "Timer is not set.\n");
//    writeErrorToLog(errorstring, params.logfile);
//    exit(1);
//  }
  if (this->params.get(intParam::DISJ_TERMS) < 2) {
    return ExitReason::NO_DISJUNCTION_EXIT;
  }

  // Set up solver
  SolverInterface* BBSolver;
  try {
    BBSolver = dynamic_cast<SolverInterface*>(si->clone());
  } catch (std::exception& e) {
    error_msg(errorstring,
        "Unable to clone solver into desired SolverInterface.\n");
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }
//  setupClpForCbc(BBSolver);

  // Setup LP for use in partial tree generation
  setLPSolverParameters(BBSolver);
  BBSolver->setHintParam(OsiDoPresolveInInitial, false);
  BBSolver->setHintParam(OsiDoPresolveInResolve, false);
  BBSolver->setIntParam(OsiMaxNumIterationHotStart, std::numeric_limits<int>::max());
#ifdef USE_CLP
  try {
    dynamic_cast<OsiClpSolverInterface*>(BBSolver)->setSpecialOptions(16); // use standard strong branching rather than clp's
    // Do not switch from dual to primal, or something to this effect;
    // This allows infeasible branches to be fixed during strong branching
    dynamic_cast<OsiClpSolverInterface*>(BBSolver)->getModelPtr()->setMoreSpecialOptions(
        dynamic_cast<OsiClpSolverInterface*>(BBSolver)->getModelPtr()->moreSpecialOptions() + 256);
  } catch (std::exception& e) {
    std::cerr << "Unable to cast solver as OsiClpSolverInterface." << std::endl;
    exit(1);
  }
#endif

  printf("## Generating partial branch-and-bound tree. ##\n");
  CbcModel* cbc_model = new CbcModel; // for obtaining the disjunction
  cbc_model->swapSolver(BBSolver);
  cbc_model->setModelOwnsSolver(true); // solver will be deleted with cbc object
  setIPSolverParameters(cbc_model, params.get(VERBOSITY));

  if (timer)
    timer->start_timer(CglVPC::VPCTimeStatsName[static_cast<int>(CglVPC::VPCTimeStats::DISJ_GEN_TIME)]);
  int num_strong = params.get(intParam::PARTIAL_BB_NUM_STRONG);
  if (num_strong == -1) {
    num_strong = si->getNumCols();
  }
  if (num_strong == -2) {
    num_strong = static_cast<int>(std::ceil(std::sqrt(si->getNumCols())));
  }
  const int num_before_trusted = std::numeric_limits<int>::max(); // 10;
  generatePartialBBTree(this, cbc_model, si, params.get(intParam::DISJ_TERMS),
      num_strong, num_before_trusted);
  printf("Finished generating partial branch-and-bound tree");
  if (timer) {
    timer->end_timer(CglVPC::VPCTimeStatsName[static_cast<int>(CglVPC::VPCTimeStats::DISJ_GEN_TIME)]);
    printf(" (%.3f seconds)", timer->get_time(CglVPC::VPCTimeStatsName[static_cast<int>(CglVPC::VPCTimeStats::DISJ_GEN_TIME)]));
  }
  printf("\n");

  VPCEventHandler* eventHandler = NULL;
  try {
    eventHandler = dynamic_cast<VPCEventHandler*>(cbc_model->getEventHandler());
  } catch (std::exception& e) {
  }
  if (!eventHandler) {
    error_msg(errstr, "Could not get event handler.\n");
    writeErrorToLog(errstr, params.logfile);
    exit(1);
  }
  
#ifdef TRACE
  const int TEMP_VAL = params.get(intParam::TEMP);
  if (std::abs(TEMP_VAL) >= static_cast<int>(TempOptions::GEN_TIKZ_STRING_WITH_VPCS)
      && std::abs(TEMP_VAL) <= static_cast<int>(TempOptions::GEN_TIKZ_STRING_AND_EXIT)) {
    generateTikzTreeString(eventHandler, params, params.get(intParam::PARTIAL_BB_STRATEGY), si->getObjValue(), true);
    if (std::abs(TEMP_VAL) == static_cast<int>(TempOptions::GEN_TIKZ_STRING_AND_RETURN)) {
      // Free
      if (cbc_model) { delete cbc_model; }
      return ExitReason::SUCCESS_EXIT;
    }
    if (std::abs(TEMP_VAL) == static_cast<int>(TempOptions::GEN_TIKZ_STRING_AND_EXIT)) {
      exit(1); // this is during debug and does not free memory
    }
  }
#endif

  // If branch-and-bound finished (was not stopped by a user event), check why and exit
  if (cbc_model->status() == 0 || cbc_model->status() == 1
      || cbc_model->status() == 2 || this->num_terms <= 1) {
    if (cbc_model->isProvenOptimal()) {
      // Set strong branching lb and ub to both be this value
      this->best_obj = cbc_model->getObjValue();
      this->worst_obj = cbc_model->getObjValue();
      const double* sol = cbc_model->getColSolution();
      this->integer_sol.assign(sol, sol + si->getNumCols());
    } else {
      warning_msg(warnstr,
          "Giving up on getting cuts from the partial branch-and-bound tree (bad status or too few terms). Model status is %d.\n",
          cbc_model->status());
    }

    return ExitReason::PARTIAL_BB_OPTIMAL_SOLUTION_FOUND_EXIT;
  } // exit out early if cbc_model status is 0 or insufficiently many disjunctive terms

  // Make sure that the right number of terms has been saved
  if ((num_terms != eventHandler->getNumLeafNodes())
      || (num_terms != static_cast<int>(terms.size() + eventHandler->isIntegerSolutionFound()))) {
    error_msg(errstr,
        "Number of terms does not match: num terms = %d, num leaf nodes = %d, num bases = %d, found_integer_sol = %d\n",
        num_terms, eventHandler->getNumLeafNodes(), static_cast<int>(terms.size()),
        eventHandler->isIntegerSolutionFound());
    writeErrorToLog(errstr, params.logfile);
    exit(1);
  }

#ifdef TRACE
  std::vector<NodeStatistics> stats = eventHandler->getStatsVector();
  printNodeStatistics(stats, false);
  if (eventHandler->getPrunedStatsVector().size() > 0) {
    printf("\n");
    printNodeStatistics(eventHandler->getPrunedStatsVector(), false);
  }
#endif

  if (cbc_model) { delete cbc_model; }
  return ExitReason::SUCCESS_EXIT;
} /* prepareDisjunction */

/****************** PROTECTED **********************/
void PartialBBDisjunction::initialize(const PartialBBDisjunction* const source,
    const VPCParameters* const params) {
  if (params != NULL) {
    setParams(*params);
  }
  if (source) {
    Disjunction::initialize(source);
    if (!params) {
      setParams(source->params);
    }
    this->data = source->data;
  } else {
    setupAsNew();
  }
} /* initialize */

#ifdef USE_CBC
/**
 * Set parameters for Cbc used for VPCs, as well as the custom branching decision
 */
void setCbcParametersForPartialBB(const VPCParameters& params,
    CbcModel* const cbc_model, CbcEventHandler* eventHandler,
    const int numStrong, const int numBeforeTrusted, const double max_time) {
  setIPSolverParameters(cbc_model);
  cbc_model->solver()->setIntParam(OsiMaxNumIterationHotStart, 100);

  // What is the partial strategy?
  const int strategy = std::abs(params.get(intParam::PARTIAL_BB_STRATEGY));
  const int sign = (params.get(intParam::PARTIAL_BB_STRATEGY) < 0) ? -1 : 1;
  const int compare_strategy = strategy % 10; // ones digit
  const int branch_strategy = (strategy % 100 - compare_strategy) / 10; // tens digit
  const int choose_strategy = sign * (strategy % 1000 - compare_strategy - 10 * branch_strategy) / 100; // hundreds digit

  // Branching decision (tens digit)
  // Given a branching variable, which direction to choose?
  // (The choice of branching variable is through OsiChooseVariable)
  // 0: default, 1: dynamic, 2: strong, 3: none
  CbcBranchDecision* branch;
  if (branch_strategy == 1) {
    branch = new CbcBranchDynamicDecision();
  } else if (branch_strategy == 2) {
    branch = new CbcBranchStrongDecision();
  } else if (branch_strategy == 3) {
    branch = NULL;
  } else {
    branch = new CbcBranchDefaultDecision();
  }

  // Set comparison for nodes (ones digit)
  // Given a tree, which node to pick next?
  // 0: default: 1: bfs, 2: depth, 3: estimate, 4: objective, 5: objective_reverse
  CbcCompareBase* compare;
  if (compare_strategy == 1) {
    compare = new CbcCompareBFS();
  } else if (compare_strategy == 2) {
    compare = new CbcCompareDepth();
  } else if (compare_strategy == 3) {
    compare = new CbcCompareEstimate();
  } else if (compare_strategy == 4) {
    compare = new CbcCompareObjective();
  } else {
    compare = new CbcCompareDefault();
  }

  cbc_model->setTypePresolve(0);
  cbc_model->setMaximumSeconds(max_time);
  cbc_model->setMaximumCutPassesAtRoot(0);
  cbc_model->setMaximumCutPasses(0);
  cbc_model->setWhenCuts(0);
  if (numStrong >= 0) {
    // Maximum number of strong branching candidates to consider each time
    cbc_model->setNumberStrong(numStrong);
  }
  if (numBeforeTrusted >= 0) {
    // # before switching to pseudocosts, I think; 0 disables dynamic strong branching, doesn't work well
    cbc_model->setNumberBeforeTrust(numBeforeTrusted);
  }

  if (branch) {
    OsiChooseStrongCustom choose;
    if (numStrong >= 0) {
      choose.setNumberStrong(numStrong);
    }
    if (numBeforeTrusted >= 0) {
      choose.setNumberBeforeTrusted(numBeforeTrusted);
    }
    choose.setMethod(choose_strategy);
    branch->setChooseMethod(choose);
    // From CbcModel::convertToDynamic, we see that the branching decision may be ignored without a choose method
    if ((branch->whichMethod()&1) == 0 && !branch->chooseMethod()) {
      OsiChooseStrong choose;
      choose.setNumberStrong(5);
      choose.setNumberBeforeTrusted(0);
      branch->setChooseMethod(choose);
    }
    cbc_model->setBranchingMethod(*branch);
  } /* check that branch is not NULL */

  if (eventHandler) {
    cbc_model->passInEventHandler(eventHandler);
  }

  cbc_model->setNodeComparison(compare);

  if (branch) {
    delete branch;
  }
  if (compare) {
    delete compare;
  }
} /* setCbcParametersForPartialBB */

/************************************************************/
/**
 * Generate a partial branch-and-bound tree with at most max_leaf_nodes leaf nodes
 */
void generatePartialBBTree(PartialBBDisjunction* const owner, CbcModel* cbc_model,
    const OsiSolverInterface* const solver, const int max_leaf_nodes,
    const int num_strong, const int num_before_trusted) {
  //  const double partial_timelimit = 100 * max_leaf_nodes * solver->getNumCols()
  //      * GlobalVariables::timeStats.get_time(GlobalConstants::INIT_SOLVE_TIME);
  const double partial_timelimit = owner->params.get(PARTIAL_BB_TIMELIMIT); // will be checked manually by the eventHandler

  // Set up options
  VPCEventHandler* eventHandler = new VPCEventHandler(owner, max_leaf_nodes, partial_timelimit);
  eventHandler->setOriginalSolver(solver);

  // This sets branching decision, event handling, etc.
  setCbcParametersForPartialBB(owner->params, cbc_model, eventHandler, num_strong,
      num_before_trusted, std::numeric_limits<double>::max());

#ifdef TRACE
  cbc_model->branchAndBound(3);
#else
  cbc_model->branchAndBound(0);
#endif

  // Free
  // When eventHandler is passed, Cbc currently clones it and does not delete the original
  if (eventHandler) { // in case this behavior gets changed in the future
    delete eventHandler;
  }
} /* generatePartialBBTree */
#endif
