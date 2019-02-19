//============================================================================
// Name        : VPCEventHandler.cpp
// Author      : A. M. Kazachkov
// Version     : 2018-Dec-25
// Description : Custom event handler for Cbc
//============================================================================

// For saving the information on the tree
#include "VPCEventHandler.hpp"
#include "utility.hpp"

#include <CbcModel.hpp>
#include <CbcTree.hpp>
#include <OsiAuxInfo.hpp> // solver characteristics OsiBab
#include <CbcFeasibilityBase.hpp>
#include <CbcSimpleInteger.hpp>
#include <CbcBranchDynamic.hpp> // CbcDynamicPseudoCostBranchingObject

/**
 * Set statistics when a child will be counted
 * The ids and node numbers matching correctly with Cbc's statistics strongly depends on where we call this from
 */
void setNodeStatistics(NodeStatistics& stats, const CbcNode* const node,
    const CbcModel* const model, const std::vector<NodeStatistics>& stats_vec,
    const std::vector<double>& originalLB, const std::vector<double>& originalUB,
    const bool will_create_child, const bool solution_is_found) {
  CbcNodeInfo* nodeInfo = node->nodeInfo();
  if (!nodeInfo) { // say infeasible
    return;
  }
  CbcNodeInfo* parentNodeInfo = nodeInfo->parent();
  if (parentNodeInfo) {
    stats.parent_id = parentNodeInfo->nodeNumber();
  }

  // Cbc is not consistent with way; "down" could be a 0 or a -1
  const OsiBranchingObject* branch = node->branchingObject();
  stats.branch_index = branch->branchIndex();
  const OsiTwoWayBranchingObject* osi_branch = dynamic_cast<const OsiTwoWayBranchingObject*>(branch);
  const CbcBranchingObject* cbc_branch = dynamic_cast<const CbcBranchingObject*>(branch);

  double prevLB, prevUB;
  if (osi_branch) {
    stats.variable = osi_branch->columnNumber();
    stats.way = (stats.branch_index == 0) ? osi_branch->firstBranch() : !osi_branch->firstBranch();

    const OsiSimpleInteger* obj =
        dynamic_cast<const OsiSimpleInteger*>(branch->originalObject());
    prevLB = obj->originalLowerBound();
    prevUB = obj->originalUpperBound();
  } /* osi_branch */
  else {
    stats.way = node->way();
    stats.variable = cbc_branch->variable();

    const CbcSimpleInteger* obj =
        dynamic_cast<const CbcSimpleInteger*>(branch->originalObject());
    if (obj) {
      prevLB = obj->originalLowerBound();
      prevUB = obj->originalUpperBound();
    } else {
      const CbcDynamicPseudoCostBranchingObject* dyn_branch =
          dynamic_cast<const CbcDynamicPseudoCostBranchingObject*>(branch);
      if (dyn_branch && dyn_branch->object()) {
        prevLB = dyn_branch->object()->originalLowerBound();
        prevUB = dyn_branch->object()->originalUpperBound();
      } else {
        // We get here, for example, if model has branchingMethod without a chooseMethod and numberBeforeTrusted = 0,
        // or there is no branchingMethod at all
        // I am not sure the below is right, but let us just pull this information from the solver for now
        prevLB = model->solver()->getColLower()[stats.variable];
        prevUB = model->solver()->getColUpper()[stats.variable];
      }
    }
  } /* cbc_branch */

  if (stats.branch_index == 0) {
    stats.id = nodeInfo->nodeNumber();
    stats.orig_id = stats.id;
    stats.number = model->getNodeCount();
    changedBounds(stats, model->solver(), originalLB, originalUB);
  } else {
    // We are coming from a place where this has not been incremented yet after child will be created
    stats.id = model->getNodeCount2() + will_create_child;
    //stats.number = stats_vec[node->nodeNumber()].number;
    int prev_stat_id = nodeInfo->nodeNumber();
    if (prev_stat_id < 0) {
      error_msg(errorstring, "Branch index is 1 but node has no previous number.\n");
      //writeErrorToLog(errorstring, params->logfile);
      exit(1);
    }
    /*
    if (prev_stat_id < 0) { // I don't think we ever get here...
      const NodeStatistics origStat = stats_vec[stats_vec.size() - 1];
      prev_stat_id = origStat.parent_id;
    }
    */
    stats.orig_id = prev_stat_id;
    stats.number = stats_vec[prev_stat_id].number;
  }

  stats.value = branch->value();
  stats.lb = (stats.way > 0) ? std::ceil(branch->value()) : prevLB;
  stats.ub = (stats.way > 0) ? prevUB : std::floor(branch->value());
  //stats.lb = model->solver()->getColLower()[stats.variable];
  //stats.ub = model->solver()->getColUpper()[stats.variable];

  stats.depth = node->depth();
  stats.obj = node->objectiveValue();
  if (model->getObjValue() < 1e30) {
    stats.integer_obj = model->getObjValue();
  }
  stats.found_integer_solution = solution_is_found;

  if (!stats_vec.empty()) {
    assert(stats.id == stats_vec[stats_vec.size()-1].id + 1);
  }
} /* setNodeStatistics */

/**
 * Set statistics when a child will not be counted (but the node number does go up)
 * Since this node is pruned, we only set id, parent_id, number
 * For the integer-feasible nodes, we record a bit more
 */
void setPrunedNodeStatistics(NodeStatistics& stats,
    const CbcNode* const parent_node, const CbcModel* const model,
    //const std::vector<NodeStatistics>& stats_vec,
    const bool prune_by_integrality, const std::vector<double>& originalLB,
    const std::vector<double>& originalUB, const double obj, const bool solution_is_found) {
  CbcNodeInfo* parentNodeInfo = parent_node->nodeInfo();
  if (!parentNodeInfo) { // say infeasible
    return;
  }

  stats.id = model->getNodeCount2();
  stats.orig_id = stats.id;
  stats.number = model->getNodeCount();
  stats.parent_id = parentNodeInfo->nodeNumber();
  stats.depth = parent_node->depth() + 1;
  if (prune_by_integrality) {
    stats.integer_obj = obj; // when a sol is found during sb, model has not yet been updated to reflect the value
  } else if (model->getObjValue() < 1e30) {
    stats.integer_obj = model->getObjValue();
  }
  stats.found_integer_solution = solution_is_found || prune_by_integrality;

  stats.obj = obj;
  if (prune_by_integrality) {
    changedBounds(stats, model->solver(), originalLB, originalUB);
  }
} /* setPrunedNodeStatistics */

void printNodeStatistics(const std::vector<NodeStatistics>& stats) {
  printNodeStatistics(stats, true, stdout);
} /* printNodeStatistics (useful to have for debug purposes) */

void printNodeStatistics(const std::vector<NodeStatistics>& stats,
    const bool print_bounds, FILE* myfile) {
  if (stats.size() == 0 || myfile == NULL) {
    return;
  }

  const int WIDTH = 10;
  const std::string obj = stringValue(stats[0].obj, "%-1.3f");
  const int OBJ_WIDTH = obj.length() >= WIDTH - 1 ? obj.length() + 2 : WIDTH;
  const std::string int_obj = stringValue(stats[0].integer_obj, "%-1.3f");
  const int INT_OBJ_WIDTH = int_obj.length() >= WIDTH - 1 ? int_obj.length() + 2 : WIDTH;
  const std::string val = stringValue(stats[0].value, "%-1.3f");
  const int VAL_WIDTH = val.length() >= WIDTH - 1 ? val.length() + 2 : WIDTH;
  const std::string lb = stringValue(stats[0].lb, "%.0f");
  const int LB_WIDTH = lb.length() + 2;
  const std::string ub = stringValue(stats[0].ub, "%.0f");
  const int UB_WIDTH = ub.length() + 2;

  fprintf(myfile, "%-*.*s%-*.*s%-*.*s%-*.*s%-*.*s%-*.*s%-*.*s%-*.*s%-*.*s%-*.*s%-*.*s%-*.*s%-*.*s%-*.*s",
          WIDTH, WIDTH, "id",
          WIDTH, WIDTH, "parent",
          WIDTH, WIDTH, "depth",
          WIDTH, WIDTH, "variable",
          WIDTH, WIDTH, "branch",
          WIDTH, WIDTH, "way",
          OBJ_WIDTH, OBJ_WIDTH, "obj",
          VAL_WIDTH, VAL_WIDTH, "value",
          LB_WIDTH, LB_WIDTH, "lb",
          UB_WIDTH, UB_WIDTH, "ub",
          WIDTH, WIDTH, "orig_id",
          WIDTH, WIDTH, "number",
          INT_OBJ_WIDTH, INT_OBJ_WIDTH, "int_obj",
          WIDTH, WIDTH, "int_found");
  if (print_bounds) {
    fprintf(myfile, "%-*.*s", WIDTH, WIDTH, "bounds");
  }
  fprintf(myfile, "\n");
  for (int i = 0; i < (int) stats.size(); i++) {
    printNodeStatistics(stats[i], print_bounds, myfile);
  }
} /* printNodeStatistics (vector) */

void printNodeStatistics(const NodeStatistics& stats, const bool print_bounds,
    FILE* myfile) {
  if (myfile == NULL) {
    return;
  }
  const int WIDTH = 10;
  const std::string obj = stringValue(stats.obj, "%-1.3f");
  const int OBJ_WIDTH = obj.length() >= WIDTH - 1 ? obj.length() + 2 : WIDTH;
  const std::string int_obj = stringValue(stats.integer_obj, "%-1.3f");
  const int INT_OBJ_WIDTH = int_obj.length() >= WIDTH - 1 ? int_obj.length() + 2 : WIDTH;
  const std::string val = stringValue(stats.value, "%-1.3f");
  const int VAL_WIDTH = val.length() >= WIDTH - 1 ? val.length() + 2 : WIDTH;
  const std::string lb = stringValue(stats.lb, "%.0f");
  const int LB_WIDTH = lb.length() + 2;
  const std::string ub = stringValue(stats.ub, "%.0f");
  const int UB_WIDTH = ub.length() + 2;

  fprintf(myfile, "% -*d% -*d% -*d% -*d% -*d%-*.*s%-*.*s%-*.*s%-*.*s%-*.*s% -*d% -*d%-*.*s% -*d",
          WIDTH, stats.id,
          WIDTH, stats.parent_id,
          WIDTH, stats.depth,
          WIDTH, stats.variable,
          WIDTH, stats.branch_index,
          WIDTH, WIDTH, ((stats.way <= 0) ? "down" : "up"),
          OBJ_WIDTH, OBJ_WIDTH, obj.c_str(),
          VAL_WIDTH, VAL_WIDTH, val.c_str(),
          LB_WIDTH, LB_WIDTH, lb.c_str(),
          UB_WIDTH, UB_WIDTH, ub.c_str(),
          WIDTH, stats.orig_id,
          WIDTH, stats.number,
          INT_OBJ_WIDTH, INT_OBJ_WIDTH, int_obj.c_str(),
          WIDTH, (int) stats.found_integer_solution);
  if (print_bounds) {
    fprintf(myfile, "%-*s", WIDTH, stats.changed_bounds.c_str());
  }
  fprintf(myfile, "\n");
} /* printNodeStatistics (single) */

/**
 * We are assuming here that none of the original columns have been deleted or anything,
 * and that we are only tracking the bound changes for the original columns
 */
void changedBounds(NodeStatistics& stats, const OsiSolverInterface* const solver,
    const std::vector<double>& originalLB, const std::vector<double>& originalUB) {
  stats.changed_bounds = "";
  const double* lower = solver->getColLower();
  const double* upper = solver->getColUpper();
  const int num_cols = originalLB.size();
  for (int i = 0; i < num_cols; i++) {
    const std::string variable = stringValue(i);
    if (!isVal(lower[i], originalLB[i])) {
      const std::string value = stringValue(lower[i], "%.0f");
      stats.changed_bounds += "(" + variable + ", l, " + value + "), ";
      stats.changed_var.push_back(i);
      stats.changed_bound.push_back(0);
      stats.changed_value.push_back(lower[i]);
    }
    if (!isVal(upper[i], originalUB[i])) {
      const std::string value = stringValue(upper[i], "%.0f");
      stats.changed_bounds += "(" + variable + ", u, " + value + "), ";
      stats.changed_var.push_back(i);
      stats.changed_bound.push_back(1);
      stats.changed_value.push_back(upper[i]);
    }
  }
} /* changedBounds (solver version) */

/************************************************************/
/* Implement the methods for VPCEventHandler which we use to exit early and save information */
VPCEventHandler::VPCEventHandler () : CbcEventHandler() {
  copyOurStuff(NULL);
} /* default constructor */

VPCEventHandler::VPCEventHandler (const int maxNumLeafNodes, const double maxTime, const VPCParameters* params) : CbcEventHandler() {
  copyOurStuff(NULL);
  this->params = params;
  this->maxNumLeafNodes_ = maxNumLeafNodes;
  this->maxTime_ = maxTime;
} /* VPCEventHandler-specific constructor */

VPCEventHandler::VPCEventHandler (const VPCEventHandler & rhs) : CbcEventHandler(rhs) {
  copyOurStuff(&rhs);
} /* copy constructor */

VPCEventHandler::VPCEventHandler(CbcModel * model) : CbcEventHandler(model) {
  try {
    copyOurStuff(dynamic_cast<VPCEventHandler*>(model->getEventHandler()));
  } catch (std::exception& e) {
    copyOurStuff(NULL);
  }
} /* constructor with pointer to model */

VPCEventHandler::~VPCEventHandler () {
  if (originalSolver_) {
    delete originalSolver_;
  }
  for (int i = 0; i < (int) bases_.size(); i++) {
    if (bases_[i]) {
      delete bases_[i];
    }
  }
  bases_.resize(0);
  //integerFeasibleSolutions_.resize(0);
} /* destructor */

VPCEventHandler& VPCEventHandler::operator=(const VPCEventHandler& rhs) {
  if (this != &rhs) {
    CbcEventHandler::operator=(rhs);
  }
  copyOurStuff(&rhs);
  return *this;
} /* assignment */

CbcEventHandler* VPCEventHandler::clone() const {
  return new VPCEventHandler(*this);
} /* clone */

void VPCEventHandler::copyOurStuff(const VPCEventHandler* const rhs) {
  if (rhs) {
    this->params = rhs->params;
    this->maxNumLeafNodes_ = rhs->maxNumLeafNodes_;
    this->maxTime_ = rhs->maxTime_;
    this->numLeafNodes_ = rhs->numLeafNodes_;
    this->numNodesOnTree_ = rhs->numNodesOnTree_;
    this->numNodes_ = rhs->numNodes_;
//    this->originalBasis_ = dynamic_cast<CoinWarmStartBasis*>(rhs->originalBasis_->clone());
    this->originalSolver_ = dynamic_cast<SolverInterface*>(rhs->originalSolver_->clone());
    this->originalLB_ = rhs->originalLB_;
    this->originalUB_ = rhs->originalUB_;
    this->bases_ = rhs->bases_;
    this->stats_ = rhs->stats_;
    this->pruned_stats_ = rhs->pruned_stats_;
    this->finalNodeIndices_ = rhs->finalNodeIndices_;
    this->savedSolution_ = rhs->savedSolution_;
    /*
    this->lb = rhs->lb;
    this->ub = rhs->ub;
    */
  } else {
    this->maxNumLeafNodes_ = 0;
    this->maxTime_ = 0;
    this->numLeafNodes_ = 0;
    this->numNodesOnTree_ = 0;
    this->numNodes_ = 0;
//    this->originalBasis_ = NULL;
    this->originalSolver_ = NULL;
  }
} /* copyOurStuff */

/*
 * Returns CbcAction based on one of the following CbcEvents
    enum CbcEvent {
      *! Processing of the current node is complete.
      node = 200,
      *! A tree status interval has arrived.
      treeStatus,
      *! A solution has been found.
      solution,
      *! A heuristic solution has been found.
      heuristicSolution,
      *! A solution will be found unless user takes action (first check).
      beforeSolution1,
      *! A solution will be found unless user takes action (thorough check).
      beforeSolution2,
      *! After failed heuristic.
      afterHeuristic,
      *! On entry to small branch and bound.
      smallBranchAndBound,
      *! After a pass of heuristic.
      heuristicPass,
      *! When converting constraints to cuts.
      convertToCuts,
      *! End of search.
      endSearch
    } ;
 */
CbcEventHandler::CbcAction
VPCEventHandler::event(CbcEvent whichEvent) {
//  printf("\n####################### WHERE FROM: %d #######\n", (int) whichEvent);
  if (reachedEnd_ || numLeafNodes_ >= maxNumLeafNodes_) {
    return stop;
  }

  if (whichEvent == treeStatus) {
    // Update statistics _only_ for the first node
    if (model_->getNodeCount2() == 1) {
      CbcNode* node = model_->currentNode(); // points to node just created, I think
      NodeStatistics currNodeStats;
      setNodeStatistics(currNodeStats, node, model_, stats_, originalLB_, originalUB_, true, foundSolution_);
      stats_.push_back(currNodeStats);
      numNodes_ = 1;

      // If any bounds were changed, let us update the originalLB and originalUB
      const int num_changed_bounds = stats_[0].changed_var.size();
      for (int i = 0; i < num_changed_bounds; i++) {
        if (stats_[0].changed_bound[i] <= 0) {
          originalLB_[stats_[0].changed_var[i]] = stats_[0].changed_value[i];
        } else {
          originalUB_[stats_[0].changed_var[i]] = stats_[0].changed_value[i];
        }
      }
    }

    // Save the nodes currently on the tree so we can figure out the parent easily
    numNodesOnTree_ = model_->tree()->size();
    numLeafNodes_ = isIntegerSolutionFound(); // the best integer feasible solution should be a leaf
    for (int i = 0; i < numNodesOnTree_; i++) {
      numLeafNodes_ += model_->tree()->nodePointer(i)->nodeInfo()->numberBranchesLeft();
    }
#ifdef TRACE
    printf("\n## Number nodes on tree: %d. Number leaf nodes: %d. Depth: %d. ##\n", numNodesOnTree_, numLeafNodes_, model_->currentDepth());
#endif

    currentNodes_.clear();
    pruneNode_ = 0;
    foundSolution_ = false;

    // Check if we are done
    if (numLeafNodes_ >= maxNumLeafNodes_ || model_->getCurrentSeconds() >= maxTime_) {
#ifdef TRACE
      printf(
          "\n## Reached limit (leaf nodes = %d with limit = %d; time = %.2f with limit = %.2f). Exiting. ##\n",
          numLeafNodes_, maxNumLeafNodes_,
          model_->getCurrentSeconds(), maxTime_);
#endif

      // Save information
      saveInformation();
      reachedEnd_ = true;
      return stop;
    } else {
      // Else save current tree to be used during the node event
      currentNodes_.resize(numNodesOnTree_);
      for (int i = 0; i < numNodesOnTree_; i++) {
        currentNodes_[i] = model_->tree()->nodePointer(i);
      }
      if (!model_->branchingMethod() || !model_->branchingMethod()->chooseMethod()) {
        // So far, when branchingMethod exists and has a chooseMethod, there have been no problems
        // Without those, we may run into a situation in which solver is infeasible during node event
        // but actually the node exists and is added to the tree
        // e.g., with bm23 -47
        if (model_->getNodeCount2() != stats_[stats_.size() - 1].id + 1) {
          // The misplaced node will be in the pruned_stats_ vector
          NodeStatistics stats = pruned_stats_[pruned_stats_.size() - 1];
          pruned_stats_.resize(pruned_stats_.size() - 1); // delete it from pruned_stats_
          const int missed_node_number = stats.id;
          int i; // find where this misplaced node is on the tree
          for (i = 0; i < numNodesOnTree_; i++) {
            const int curr_node_number = currentNodes_[i]->nodeInfo()->nodeNumber();
            if (curr_node_number == missed_node_number) {
              break;
            }
          }
          assert (i < numNodesOnTree_); // make sure we found it

          // Now adjust the stats (need to set variable, branch_index, way, obj, value, lb, ub, bounds)
          stats.obj = currentNodes_[i]->objectiveValue();
          const OsiBranchingObject* branch = currentNodes_[i]->branchingObject();
          stats.branch_index = branch->branchIndex(); // should be 0, but be safe
          stats.value = branch->value();

          // I think this situation is only going to happen on the Cbc side of things
          const CbcBranchingObject* cbc_branch = dynamic_cast<const CbcBranchingObject*>(branch);
          assert(cbc_branch != NULL);
          stats.way = currentNodes_[i]->way(); // if it were the Osi side, way would not be set this way
          stats.variable = cbc_branch->variable();

          const CbcSimpleInteger* obj =
            dynamic_cast<const CbcSimpleInteger*>(branch->originalObject());
          double prevLB, prevUB;
          if (obj) {
            prevLB = obj->originalLowerBound();
            prevUB = obj->originalUpperBound();
          } else {
            const CbcDynamicPseudoCostBranchingObject* dyn_branch =
              dynamic_cast<const CbcDynamicPseudoCostBranchingObject*>(branch);
            if (dyn_branch && dyn_branch->object()) {
              prevLB = dyn_branch->object()->originalLowerBound();
              prevUB = dyn_branch->object()->originalUpperBound();
            } else {
              error_msg(errorstring, "Cannot access original object.\n");
              writeErrorToLog(errorstring, params->logfile);
              exit(1);
            }
          }
          stats.lb = (stats.way > 0) ? std::ceil(branch->value()) : prevLB;
          stats.ub = (stats.way > 0) ? prevUB : std::floor(branch->value());

          // For the bounds, I am not sure the solver is reflecting the node soluton any more... but let's hope
          changedBounds(stats, model_->solver(), originalLB_, originalUB_);

          // Adjust stats vector
          NodeStatistics last_stats = stats_[stats_.size() - 1];
          if (last_stats.id == stats.id) {
            // This means that the parent of the missed node had a branch left
            last_stats.id = stats.id + 1;
            stats_[stats_.size() - 1] = stats;
            stats_.push_back(last_stats);
          } else {
            stats_.push_back(stats);
          }
        }
      }
    }
  } /* (whichEvent == treeStatus) */
  else if (whichEvent == node) {
    // Update statistics after a node has been branched on
    // Deal with child
    child_ = model_->currentNode();
    // 2017-07-09: Changed from isProvenOptimal because in some cases,
    // we hit an iteration limit during hot starting in OsiChooseVariable.cpp:399
    // Might be good to check dual infeasibility too, but these subproblems
    // should not be unbounded...
    // 2017-07-23: NB: solver_feasible = false may be if child is integer-feasible
    const bool solver_feasible = !(model_->solver()->isProvenPrimalInfeasible());
#ifdef TRACE
    if (model_->solver()->isProvenDualInfeasible()) {
      error_msg(errorstring, "Did not think that dual infeasibility in the model solver could happen. Check this.\n");
      writeErrorToLog(errorstring, params->logfile);
      exit(1);
    }
#endif
    // When a solution is found, newNode is not null, but it has no branching object
    bool solution_found =
        (solver_feasible && (pruneNode_ == 0) &&
        (child_ && child_->numberUnsatisfied() == 0));
    // The node will be pruned if it exists but pruneNode_ > 0 or a solution is found
    const bool prune_node = (pruneNode_ > 0 || solution_found);
    // The node may also simply not exist... (when solver is infeasible)
    // I think that everything after !prune_node should be unnecessary... but let us leave it in
    bool will_create_child = solver_feasible && !prune_node && child_ &&
        child_->active() && child_->nodeInfo() && child_->nodeInfo()->parent() &&
        child_->branchingObject();

    if (will_create_child) {
      // I think our assumption right now is that a new child is created with branch index set to zero
      if (child_->branchingObject()->branchIndex() > 0) {
        error_msg(errorstring,
            "Our assumption that a child is created with branch index set to zero is wrong. "
            "Node number 1: %d. Node number 2: %d. Parent index: %d.\n",
            model_->getNodeCount(), model_->getNodeCount2(),
            (child_->nodeInfo()->parent() ? child_->nodeInfo()->parent()->nodeNumber() : -1));
        writeErrorToLog(errorstring, params->logfile);
        exit(1);
      }
      NodeStatistics currNodeStats;
      setNodeStatistics(currNodeStats, child_, model_, stats_, originalLB_, originalUB_, true, foundSolution_);
      stats_.push_back(currNodeStats);
      parentInfo_ = child_->nodeInfo()->parent();
      if (currNodeStats.number >= numNodes_) {
        numNodes_ = currNodeStats.number + 1;
      }
    } /* node is added to the tree */
    else {
      // Child will not be created, and parent is not available
      // Though not created, the node number does go up
      // (whether or not the stats id increases)

      // Find parent (node removed from tree)
      int tmpNumNodesOnTree = model_->tree()->size();
      parentInfo_ = NULL;
      int i, j;
      for (i = 0; i < numNodesOnTree_; i++) {
        for (j = 0; j < tmpNumNodesOnTree; j++) {
          if (model_->tree()->nodePointer(j) == currentNodes_[i]) {
            break;
          }
        }
        if (j == tmpNumNodesOnTree) {
          // We found the parent
          parentInfo_ = currentNodes_[i]->nodeInfo();
          break;
        }
      } /* find the parent */

#ifdef TRACE
      if (!parentInfo_) {
        error_msg(errorstring, "Could not find parent node.\n");
        writeErrorToLog(errorstring, params->logfile);
        exit(1);
      }
#endif

      // I am not sure we are capturing the situation when newNode is active but branch_ is null
      // This should correspond to the node being integer feasible
      // Perhaps this happens when we do not use strong branching...
      if (solver_feasible && !solution_found && (pruneNode_ == 0) && child_ && child_->active() && child_->branchingObject()) {
        error_msg(errorstring, "We should not get here, I think. Is it an integer-feasible solution? If so, it should be saved.\n");
        writeErrorToLog(errorstring, params->logfile);
        exit(1);
      }

      // When !solver_feasible, this may be because it is actually integer-feasible
      // which happens at CbcModel.cpp:14411, during resolve
      // (called from resolve, called from solveWithCuts)
      // (see email to cbc mailing list on 2017-07-28)
      // 2017-08-11: However, it could also be due to truly being infeasible
      // When infeasible, the garbage primal solution might end up being actually integer feasible...
      // E.g., bm23 -2 _after_ cuts
      // 2017-08-23: Might happen that it is truly infeasible, but the integer part is not only integer-feasible,
      // but also actually feasible, as with neos-880324_cleaned -32, so need to check all cols, not just integer objects
      if (pruneNode_ != 3 && !solution_found && !solver_feasible) {
        solution_found = true;
        for (int col = 0; col < model_->getNumCols(); col++) {
          const double val = model_->solver()->getColSolution()[col];
          const double lb = model_->solver()->getColLower()[col];
          const double ub = model_->solver()->getColUpper()[col];
          if (lessThanVal(val, lb) || greaterThanVal(val, ub)) {
            solution_found = false;
            break;
          }
          if (originalSolver_->isInteger(col)) {
            const double floor_xk = std::floor(val);
            const double infeasibility = (val - floor_xk < 0.5) ? (val - floor_xk) : (floor_xk + 1 - val);
            if (!isZero(infeasibility, params->get(doubleParam::EPS))) {
              solution_found = false;
              break;
            }
          }
        }
      } /* check if a solution found after all, when solver is "infeasible" */

      // Ensure the solution is really integer-feasible and not cut off
      double objectiveValue = model_->solver()->getObjValue(); // sometimes this does not have the right value stored
      if (pruneNode_ == 3 || (solution_found &&
          !greaterThanVal(objectiveValue, model_->getCutoff() + model_->getCutoffIncrement()))) {
        const double* sol_ptr = (pruneNode_ == 3) ? savedSolution_.data() : model_->solver()->getColSolution();
        std::vector<double> sol;
        sol.reserve(model_->getNumCols() + model_->getNumRows());
        sol.insert(sol.end(), sol_ptr, sol_ptr + model_->getNumCols());
        int error = -1;
        double violation = 0.;
        // Check bounds on variables and integer constraints satisfied
        // Also, deal with some round-off errors with integers
        for (int col_ind = 0; col_ind < model_->getNumCols(); col_ind++) {
          violation = originalSolver_->getColLower()[col_ind] - sol[col_ind];
          if (greaterThanVal(violation, 0.)) {
            error = col_ind;
            break;
          } else if (greaterThanVal(violation, 0., 0.)) {
            sol[col_ind] += violation;
          }
          violation = sol[col_ind] - originalSolver_->getColUpper()[col_ind];
          if (greaterThanVal(violation, 0.)) {
            error = col_ind;
            break;
          } else if (greaterThanVal(violation, 0., 0.)) {
            sol[col_ind] -= violation;
          }
          if (originalSolver_->isInteger(col_ind)) {
            const double val = sol[col_ind];
            const double floor_xk = std::floor(val);
            violation = (val - floor_xk < 0.5) ? val - floor_xk : val - (floor_xk + 1);
            if (!isZero(violation)) {
              error = col_ind;
              violation = std::abs(violation);
              break;
            } else if (!isZero(violation, 0.)) {
              sol[col_ind] -= violation;
            }
          }
        } /* check col bounds */

        // Check row bounds are satisfied
        // And save slack solution
        for (int row_ind = 0; row_ind < model_->getNumRows() && !error; row_ind++) {
          // We recompute activity due to possibly adjusting the variables above
          const double ref_activity = model_->solver()->getRowActivity()[row_ind];
          const CoinShallowPackedVector& vec = model_->getMatrixByRow()->getVector(row_ind);
          const int size_vec = vec.getNumElements();
          const int* vec_ind = vec.getIndices();
          const double* vec_el = vec.getElements();
          const double activity = dotProduct(size_vec, vec_ind, vec_el, sol.data());
          if (!isVal(activity, ref_activity)) {
            error_msg(errorstring, "Error calculating activity. Activity: %f. Reference activity: %f.\n", activity, ref_activity);
            writeErrorToLog(errorstring, params->logfile);
            exit(1);
          }
          double rhs = model_->solver()->getRightHandSide()[row_ind];
          violation = std::abs(rhs - activity);
          sol.push_back(violation);
          const char sense = model_->solver()->getRowSense()[row_ind];
          if (sense == 'G') {
            if (lessThanVal(activity, rhs, params->get(doubleParam::EPS))) {
              error = row_ind + model_->getNumCols();
              break;
            }
          } else if (sense == 'L') {
            if (greaterThanVal(activity, rhs, params->get(doubleParam::EPS))) {
              error = row_ind + model_->getNumCols();
              break;
            }
          } else if (sense == 'E') {
            if (!isVal(activity, rhs, params->get(doubleParam::EPS))) {
              error = row_ind + model_->getNumCols();
              break;
            }
          } else { // ranged; not 'N'
            // Need to check upper and lower bounds
            // rhs returns rowupper() in this case
            if (greaterThanVal(activity, rhs, params->get(doubleParam::EPS))) {
              error = row_ind + model_->getNumCols();
              break;
            }
            rhs = model_->solver()->getColLower()[row_ind];
            if (lessThanVal(activity, rhs, params->get(doubleParam::EPS))) {
              error = row_ind + model_->getNumCols();
              break;
            }
          }
        } /* check row bounds */

        if (error >= 0) {
          error_msg(errorstring, "An integer-feasible solution was supposedly found, but it is infeasible somehow. Variable: %d. Violation: %f.\n", error, violation);
          writeErrorToLog(errorstring, params->logfile);
          exit(1);
        }

        // Save the integer feasible solution if it cannot pruned by bound
        // We're being a little extra safe here by allowing an extra epsilon
        double objOffset = 0.;
        model_->solver()->getDblParam(OsiDblParam::OsiObjOffset, objOffset);
        objectiveValue = dotProduct(sol.data(), model_->getObjCoefficients(), model_->getNumCols()) - objOffset;
        if (pruneNode_ == 3 && !isVal(objectiveValue, model_->getObjValue())) {
          error_msg(errorstring, "We are assuming that if pruning by integrality, the objective we calculated should be the same as the objective of the solver.\n");
          writeErrorToLog(errorstring, params->logfile);
          exit(1);
        }
        if (pruneNode_ == 3 || (!greaterThanVal(objectiveValue, model_->getCutoff() + model_->getCutoffIncrement()))) {
          savedSolution_ = sol;
        }
      } /* is an integer-feasible solution found? */

      // Also update the infeasible nodes stats
      NodeStatistics prunedNodeStats;
      setPrunedNodeStatistics(prunedNodeStats, parentInfo_->owner(), model_,
      /*stats_,*/solution_found || pruneNode_ == 3, originalLB_, originalUB_,
          objectiveValue, foundSolution_);
      pruned_stats_.push_back(prunedNodeStats);
      if (prunedNodeStats.number >= numNodes_) {
        numNodes_ = prunedNodeStats.number + 1;
      }
    } /* node is pruned */

    // Update parent if it has branches left
    if (parentInfo_->numberBranchesLeft() > 0) {
      NodeStatistics currNodeStats;
      setNodeStatistics(currNodeStats, parentInfo_->owner(), model_, stats_,
          originalLB_, originalUB_, will_create_child);
      stats_.push_back(currNodeStats);
    }
  }  /* (whichEvent == node) */
  else if (whichEvent == beforeSolution2) {
    // This event may cause newNode to be deleted, so we remember not to access deleted memory
    // This action is happening in CbcModel.cpp::chooseBranch (lines 14845+ in Cbc-2.9.9)
    // If the return code is -2, the node will be deleted, so that is what we are trying to catch

    // If the obj value of the new solution (stored at bestSolution()) is < cutoff,
    // the cutoff will be updated and pruning will happen
    // (cutoff = best integer objective value - 0.999 or so;
    // used to cut any nodes with a worse objective value than this cutoff,
    // or integer solutions with worse value)
    foundSolution_ = true;
    child_ = model_->currentNode();

    double objectiveValue = model_->savedSolutionObjective(0); // eventHandler replaces bestObjective_ temporarily, so we are accessing the solution's objective here
    double cutoff = model_->getCutoff(); // this is the same as the previous bestObjective_
    double cutoffIncrement= model_->getCutoffIncrement();
    double newCutoff = objectiveValue - cutoffIncrement;
    /*
    // Should really be more careful with rounding and such
    // (see CbcModel.cpp:13130: adjust to ``allow for rounding errors'')
    if (cutoffIncrement == 1e-5) {
      saveObjectiveValue = CoinMax(saveObjectiveValue, bestObjective_ - 0.0000001 * std::abs(bestObjective_));
      cutoff = CoinMin(bestObjective_, saveObjectiveValue) - 1.0e-5;
      if (std::abs(cutoff + 1.0e-5 - std::floor(cutoff + 0.5)) < 1.0e-8)
        cutoff -= 2.0e-5;
    }
    */

    if (!model_->solverCharacteristics()->mipFeasible()
        || (model_->problemFeasibility()->feasible(model_, 0) < 0)) {
      pruneNode_ = 1;
    } /* mip infeasible */
    else if (child_->objectiveValue() >= model_->getCutoff()) {
      // CbcModel.cpp:15283
      pruneNode_ = 2;
    } /* pruned by objective */
    else if (objectiveValue < cutoff &&
        child_->objectiveValue() >= newCutoff) {
      // Check whether the new solution gives a strictly better bound on the objective
      // and whether this new value implies that the new node will be cut off
      // It is when CbcModel.cpp:13115 is true, so that cutoff is updated, then -2 returned at CbcModel.cpp:15283; e.g., bm23 -70
      // 2017-08-11: Can reach here from CbcModel.cpp:16844, e.g., with bm23 -8 after cuts, from CbcNode::chooseOsiBranch finding a feasible solution during strong branching, so solver never gets the good solution (it is only in choose)
      pruneNode_ = 3;
      savedSolution_.assign(model_->bestSolution(), model_->bestSolution() + model_->getNumCols());
    } /* mip feasible */
  } /* (whichEvent == beforeSolution2) */
  else if (whichEvent == endSearch) {
    // If we get here,
    // this means that a solution has been found prior to the maximum number of leaf nodes,
    // or perhaps a (time, node, solution) limit has been reached?
    // Desired information has not been saved
#ifdef TRACE
    printf(
        "\n## Reached end of search prematurely with model status %d. Exiting. ##\n",
        model_->status());
#endif
    if (model_->status() == 1) {
      // Should not really get here any more since we are setting time limit to twice the real time limit
      // But if we do, I guess we have to abandon the search
      numNodesOnTree_ = 0;
      numLeafNodes_ = 0;
    } /* max nodes, max sols, or max time */
    else if (model_->status() == 2) {
      // When we are here it is too late to salvage the search, I think
      numNodesOnTree_ = 0;
      numLeafNodes_ = 0;
    } /* difficulties so run was abandoned */
    else {
      assert(model_->status() == 0);
      numNodesOnTree_ = 0;
      numLeafNodes_ = 0;
    }
    reachedEnd_ = true;
  } /* (whichEvent == endSearch) */

//  printNodeStatistics(stats_);
  return noAction;
} /* VPCEventHandler::event  */

/**
 * Save all relevant information before it gets deleted by BB finishing
 */
void VPCEventHandler::saveInformation() {
  /*
  // 2017-08-11: Only one integer-feasible solution needs to be kept (the best one), so we drop this
  // Prune integer feasible solutions that are above cutoff
  // Will use a lambda function because I want to try it (with the erase-remove idiom)
  const double cutoff = model_->getCutoff() + model_->getCutoffIncrement();
  const double* obj = model_->getObjCoefficients();
  const int num_cols = model_->getNumCols();
  double objOffset = 0.;
  model_->solver()->getDblParam(OsiDblParam::OsiObjOffset, objOffset);
  integerFeasibleSolutions_.erase(
      std::remove_if(integerFeasibleSolutions_.begin(),
          integerFeasibleSolutions_.end(),
          [&](const std::vector<double>& x)
          {
            return (greaterThanVal(dotProduct(x.data(), obj, num_cols) - offset, cutoff));
          }), integerFeasibleSolutions_.end());
  */
  // For each node on the tree, save the warm start basis and branching directions
  bases_.resize(numNodesOnTree_);
  CoinWarmStartBasis* tmp = dynamic_cast<CoinWarmStartBasis*>(originalSolver_->getWarmStart());
  for (int i = 0; i < numNodesOnTree_; i++) {
    bases_[i] = dynamic_cast<CoinWarmStartBasis*>(tmp->clone());
    CbcNode* node = model_->tree()->nodePointer(i);
    CbcNodeInfo* nodeInfo = node->nodeInfo();
    nodeInfo->buildRowBasis(*bases_[i]);
    finalNodeIndices_.push_back(stats_[nodeInfo->nodeNumber()].id);

    /*
    // Might be partial or full
    try {
    } catch (std::exception& e) {

    }
    */

    /*
    // Figure out the variables fixed to get to this node
    // Work backwards via parent
    CbcNodeInfo* parent = nodeInfo->parent();
    printf("Node number: %d.\n", node->nodeNumber());
    printf("\tBranching on: %d.\n", node->branchingObject()->columnNumber());
    while (parent) {
    const CbcNode* parentNode = parent->owner();
    printf("\tFrom node %d branching on %d.\n", parentNode->nodeNumber(), parentNode->branchingObject()->columnNumber());
    parent = parentNode->nodeInfo()->parent();
    }
    printf("\n");
    */
  }
  if (tmp) {
    delete tmp;
  }

  /*
     const int numObjects = model_->numberObjects();

  // Resize arrays
  lb.resize(numObjects);
  ub.resize(numObjects);

  // Get new bounds
  //const SolverInterface* solver = dynamic_cast<SolverInterface*>(model->
  for (int i = 0; i < model_->numberObjects(); i++) {
  try {
  // Usually we will be in the (new?) Osi side
  const OsiSimpleInteger* obj = dynamic_cast<const OsiSimpleInteger*>(model_->object(i));
  lb[i] = obj->originalLowerBound();
  ub[i] = obj->originalUpperBound();
  } catch (std::exception& e) {
// Maybe we are in the Cbc side of the branching process
const CbcSimpleInteger* obj = dynamic_cast<const CbcSimpleInteger*>(model_->object(i));
lb[i] = obj->originalLowerBound();
ub[i] = obj->originalUpperBound();
}
}
*/
} /* saveInformaton */

//CoinWarmStartBasis* VPCEventHandler::setOriginalBasis(const CoinWarmStart* const copyBasis, const bool return_old) {
//  CoinWarmStartBasis* tmp;
//  if (originalBasis_ && return_old) {
//    tmp = dynamic_cast<CoinWarmStartBasis*>(originalBasis_->clone());
//  }
//  originalBasis_ = dynamic_cast<CoinWarmStartBasis*>(copyBasis->clone());
//  if (return_old) {
//    return tmp;
//  } else {
//    return NULL;
//  }
//} /* setOriginalBasis */

SolverInterface* VPCEventHandler::setOriginalSolver(
    const OsiSolverInterface* const copySolver, const bool return_old) {
  SolverInterface* tmp = NULL;
  if (return_old && originalSolver_) {
    tmp = dynamic_cast<SolverInterface*>(originalSolver_->clone());
  }
  try {
    originalSolver_ = dynamic_cast<SolverInterface*>(copySolver->clone());
  } catch (std::exception& e) {
    error_msg(errorstring, "Unable to clone solver as SolverInterface.\n");
    writeErrorToLog(errorstring, params->logfile);
    exit(1);
  }
  this->setOriginalLB(originalSolver_->getNumCols(), originalSolver_->getColLower());
  this->setOriginalUB(originalSolver_->getNumCols(), originalSolver_->getColUpper());
  if (return_old) {
    return tmp;
  } else {
    return NULL;
  }
} /* setOriginalSolver */

void VPCEventHandler::setOriginalLB(const int num_cols, const double* const vec) {
  originalLB_.assign(vec, vec + num_cols);
} /* setOriginalLB */

void VPCEventHandler::setOriginalUB(const int num_cols, const double* const vec) {
  originalUB_.assign(vec, vec + num_cols);
} /* setOriginalUB */
