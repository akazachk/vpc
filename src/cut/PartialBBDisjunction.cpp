// Name:     PartialBBDisjunction.cpp
// Author:   A. M. Kazachkov
// Date:     2018-02-22
//-----------------------------------------------------------------------------
#include "PartialBBDisjunction.hpp"

// COIN-OR
#include <CbcModel.hpp>

// Project files
#include "BBHelper.hpp"
#include "CglVPC.hpp" // get timer information too
#include "SolverHelper.hpp"
#include "utility.hpp"
#include "VPCEventHandler.hpp"

#ifdef TRACE
#include "debug.hpp"
#endif

/****************** PUBLIC  **********************/
/** Default constructor */
PartialBBDisjunction::PartialBBDisjunction() {
  initialize();
} /* default constructor */

/** Param constructor */
PartialBBDisjunction::PartialBBDisjunction(const VPCParameters& param) {
  initialize(NULL, &param);
} /* param constructor */

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

/** setParams */
void PartialBBDisjunction::setParams(const VPCParameters& param) {
  this->params = param;
} /* setParams */

ExitReason PartialBBDisjunction::prepareDisjunction(OsiSolverInterface* const si) {
  if (!timer) {
    error_msg(errorstring, "Timer is not set.\n");
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }
//  ExitReason retval = ExitReason::UNKNOWN;
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
#ifdef USE_CLP
  setupClpForCbc(BBSolver);
#endif
//
#ifdef TRACE
  printf("\n## Generating partial branch-and-bound tree. ##\n");
#endif
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
  generatePartialBBTree(params, cbc_model, si, params.get(intParam::DISJ_TERMS),
      num_strong, num_before_trusted);
  if (timer)
    timer->end_timer(CglVPC::VPCTimeStatsName[static_cast<int>(CglVPC::VPCTimeStats::DISJ_GEN_TIME)]);

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

  this->data.num_partial_bb_nodes = cbc_model->getNodeCount(); // save number of nodes looked at
  this->num_terms = eventHandler->getNumLeafNodes()
      + eventHandler->isIntegerSolutionFound();
  this->data.num_pruned_nodes = eventHandler->getPrunedStatsVector().size()
      - eventHandler->isIntegerSolutionFound();

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

//  // Save everything we need from this disjunction
  this->data.stats = eventHandler->getStatsVector();
  this->data.num_nodes_on_tree = eventHandler->getNumNodesOnTree();
  this->data.num_fixed_vars = cbc_model->strongInfo()[1]; // number fixed during b&b
  this->data.node_id.resize(this->data.num_nodes_on_tree);
  this->bases.resize(this->data.num_nodes_on_tree);
  for (int tmp_ind = 0; tmp_ind < this->data.num_nodes_on_tree; tmp_ind++) {
//    this->disjunction.node_id[tmp_ind] = eventHandler->getNodeIndex(tmp_ind);
    this->bases[tmp_ind] =
        (eventHandler->getBasisForNode(tmp_ind))->clone();
  }
  if (eventHandler->isIntegerSolutionFound()) {
    const double* sol = eventHandler->getIntegerFeasibleSolution();
    this->integer_sol.assign(sol, sol + si->getNumCols());
  }

#ifdef TRACE
  std::vector<NodeStatistics> stats = eventHandler->getStatsVector();
  printNodeStatistics(stats, false);
  if (eventHandler->getPrunedStatsVector().size() > 0) {
    printf("\n");
    printNodeStatistics(eventHandler->getPrunedStatsVector(), false);
  }

  if (std::abs(params.get(intParam::TEMP)) >= 10 && std::abs(params.get(intParam::TEMP)) <= 15) {
    generateTikzTreeString(eventHandler, params,
        params.get(intParam::PARTIAL_BB_STRATEGY),
        si->getObjValue(), true);
    if (std::abs(params.get(intParam::TEMP)) == 14) {
      // Free
      if (cbc_model) { delete cbc_model; }
      return ExitReason::SUCCESS_EXIT;
    }
    if (std::abs(params.get(intParam::TEMP)) == 15) {
      exit(1); // this is during debug and does not free memory
    }
  }
#endif

  if (cbc_model) { delete cbc_model; }
  return ExitReason::SUCCESS_EXIT;
} /* prepareDisjunction */

/** Clone */
PartialBBDisjunction* PartialBBDisjunction::clone() const {
  return new PartialBBDisjunction(*this);
} /* clone */

void PartialBBDisjunction::initialize(const PartialBBDisjunction* const source, const VPCParameters* const params) {
  Disjunction::initialize(source);
  if (params != NULL) {
    setParams(*params);
  }
  if (source) {
    if (!params) {
      setParams(source->params);
    }
    this->timer = source->timer;
    this->data = source->data;
  } else {
    this->timer = NULL;
    this->data.num_nodes_on_tree = 0;
    this->data.num_partial_bb_nodes = 0;
    this->data.num_pruned_nodes = 0;
    this->data.min_node_depth = std::numeric_limits<int>::max();
    this->data.max_node_depth = 0;
    this->data.num_fixed_vars = 0;
    this->data.stats.resize(0);
    this->data.pruned_stats.resize(0);
    this->data.node_id.resize(0);
  }
} /* initialize */
