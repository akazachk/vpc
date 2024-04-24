/**
 * @file VPCEventHandler.cpp
 * @brief Custom event handler for Cbc
 * @author A. M. Kazachkov
 * @date 2019-02-24
 */
// For saving the information on the tree
#include "VPCEventHandler.hpp"
#include "PartialBBDisjunction.hpp"
#include "SolverHelper.hpp"
#include "utility.hpp"
#include <algorithm> // for std::find
#include <memory> // unique_ptr, make_unique (C++14)

#include <OsiAuxInfo.hpp> // solver characteristics OsiBab

#ifdef USE_CBC
#include <CbcModel.hpp>
#include <CbcTree.hpp>
#include <CbcFeasibilityBase.hpp>
#include <CbcSimpleInteger.hpp>
#include <CbcBranchDynamic.hpp> // CbcDynamicPseudoCostBranchingObject
#include <OsiChooseVariable.hpp>
#endif // USE_CBC

#include "CglVPC.hpp"

using namespace VPCParametersNamespace;

/**
 * @details Set statistics when a child will be counted.
 * The ids and node numbers matching correctly with Cbc's statistics strongly depends on where we call this from.
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
  if (!branch) {
    error_msg(errorstring, "Expected branch to not be NULL in setNodeStatistics.\n");
    //writeErrorToLog(errorstring, owner->params.logfile);
    exit(1);
  }
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
      //writeErrorToLog(errorstring, owner->params.logfile);
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

  stats.status = "normal";
} /* setNodeStatistics */

/**
 * @details Set statistics when a child will not be counted (but the node number does go up).
 * 
 * Since this node is pruned, we only set id, orig_id, parent_id, number, depth, obj.
 * For the integer-feasible nodes, we record the integer_obj.
 * We store the changed bounds but this may be inaccurate.
 */
void setPrunedNodeStatistics(NodeStatistics& stats,
    const CbcNode* const parent_node, const CbcModel* const model,
    //const std::vector<NodeStatistics>& stats_vec,
    const bool prune_by_integrality, const std::vector<double>& originalLB,
    const std::vector<double>& originalUB, const double obj,
    const bool solution_is_found,
    const VPCEventHandler::PruneNodeOption& pruneReason) {
  CbcNodeInfo* parentNodeInfo = parent_node->nodeInfo();
  if (!parentNodeInfo) { // say infeasible
    return;
  }

  stats.id = model->getNodeCount2();
  stats.orig_id = stats.id;
  stats.number = model->getNodeCount();
  stats.parent_id = parentNodeInfo->nodeNumber();
  stats.depth = parent_node->depth() + 1;
  stats.branch_index = 0;
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

  stats.status = "unknown";
  if (pruneReason == VPCEventHandler::PruneNodeOption::PRUNE_BY_INFEASIBILITY) {
    stats.status = "infeas";
  }
  else if (pruneReason == VPCEventHandler::PruneNodeOption::PRUNE_BY_BOUND) {
    stats.status = "bound";
  }
  else if (pruneReason == VPCEventHandler::PruneNodeOption::PRUNE_BY_INTEGRALITY) {
    stats.status = "integ";
  }
} /* setPrunedNodeStatistics */

void printNodeStatistics(const std::vector<NodeStatistics>& stats) {
  printNodeStatistics(stats, true, stdout);
} /* printNodeStatistics (useful to have for debug purposes) */

/// @details Prints header for table of node statistics,
/// then calls printNodeStatistics(const NodeStatistics&, const bool, FILE*).
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

  fprintf(myfile, "%-*.*s%-*.*s%-*.*s%-*.*s%-*.*s%-*.*s%-*.*s%-*.*s%-*.*s%-*.*s%-*.*s%-*.*s%-*.*s%-*.*s%-*.*s",
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
          WIDTH, WIDTH, "int_found",
          WIDTH, WIDTH, "status");
  if (print_bounds) {
    fprintf(myfile, "%-*.*s", WIDTH, WIDTH, "bounds");
  }
  fprintf(myfile, "\n");
  for (int i = 0; i < (int) stats.size(); i++) {
    printNodeStatistics(stats[i], print_bounds, myfile);
  }
} /* printNodeStatistics (vector) */

/// @details Must have width of each statistic printed line up with the header (assumed to be) printed earlier
void printNodeStatistics(
    /// [in] node that we print statistics for
    const NodeStatistics& stats, 
    /// [in] should the string NodeStatistics::changed_bounds be printed?
    const bool print_bounds,
    /// [in] output file
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

  fprintf(myfile, "% -*d% -*d% -*d% -*d% -*d%-*.*s%-*.*s%-*.*s%-*.*s%-*.*s% -*d% -*d%-*.*s% -*d%-*.*s",
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
          WIDTH, (int) stats.found_integer_solution,
          WIDTH, WIDTH, stats.status.c_str());
  if (print_bounds) {
    fprintf(myfile, "%-*s", WIDTH, stats.changed_bounds.c_str());
  }
  fprintf(myfile, "\n");
} /* printNodeStatistics (single) */

/**
 * @details We are assuming here that none of the original columns have been deleted or anything,
 * and that we are only tracking the bound changes for the original columns.
 * The result may be incorrect for pruned nodes when, for example, an integer solution is found during strong branching.
 */
void changedBounds(
    /// Current node statistics
    NodeStatistics& stats,
    /// Current solver state in branch-and-bound
    const OsiSolverInterface* const solver,
    /// Original lower bounds on variables
    const std::vector<double>& originalLB,
    /// Original upper bounds on variables
    const std::vector<double>& originalUB) {
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
// <!---------- Implement the methods for VPCEventHandler which we use to exit early and save information ---------->
VPCEventHandler::VPCEventHandler () : CbcEventHandler() {
  initialize(NULL);
} /* default constructor */

VPCEventHandler::VPCEventHandler(PartialBBDisjunction* const disj,
    const int maxNumLeafNodes, const double maxTime) :
    CbcEventHandler() {
  initialize(NULL);
  this->owner = disj;
  this->maxNumLeafNodes_ = maxNumLeafNodes;
  this->maxTime_ = maxTime;
} /* VPCEventHandler-specific constructor */

VPCEventHandler::VPCEventHandler (const VPCEventHandler & rhs) : CbcEventHandler(rhs) {
  initialize(&rhs);
} /* copy constructor */

VPCEventHandler::VPCEventHandler(CbcModel * model) : CbcEventHandler(model) {
  try {
    initialize(dynamic_cast<VPCEventHandler*>(model->getEventHandler()));
  } catch (std::exception& e) {
    initialize(NULL);
  }
} /* constructor with pointer to model */

VPCEventHandler::~VPCEventHandler () {
  if (originalSolver_) {
    delete originalSolver_;
    originalSolver_ = NULL;
  }
//  if (cuts_) {
//    delete cuts_;
//    cuts_ = NULL;
//  }
//  for (int i = 0; i < (int) bases_.size(); i++) {
//    if (bases_[i]) {
//      delete bases_[i];
//    }
//  }
//  bases_.resize(0);
  //integerFeasibleSolutions_.resize(0);
} /* destructor */

VPCEventHandler& VPCEventHandler::operator=(const VPCEventHandler& rhs) {
  if (this != &rhs) {
    CbcEventHandler::operator=(rhs);
  }
  initialize(&rhs);
  return *this;
} /* assignment */

CbcEventHandler* VPCEventHandler::clone() const {
  return new VPCEventHandler(*this);
} /* clone */

/**
 * @details We need to track:
 *   1. the dual bound at each node
 *   2. branching variable and direction at each node
 *   3. which nodes are pruned (by bound, infeasibility, or integrality)
 *   4. the best integer-feasible solution found and its objective value
 *       -- only the best, because all others are effectively pruned by bound then
 *
 * We will stop when either the limit on leaf nodes (disjunctive terms) or time is reached.
 *
 * Questions / Challenges / Comments:
 *  * Is it necessary to track beforeSolution2? Maybe solution is enough...
 *  * If at any point (newNode->objectiveValue() >= getCutoff()), then that node will be pruned
 *    We have access to newNode in only a few places...
 *  * anyAction = chooseBranch

  Returns CbcAction based on one of the following CbcEvents
    enum CbcEvent {
      *! Processing of the current node is complete. amk: Called within CbcModel::doOneNode.
      node = 200,
      *! A tree status interval has arrived.
      *! amk: Called within CbcModel::branchAndBound.
      *! amk: At this point, we know the next step will be picking the next leaf node to branch on using tree_->bestNode(cutoff);
      *! amk: We may not reach a treeStatus interval if stoppingCriterionReached() is true after doOneNode is called (below the treeStatus event).
      *! amk: If stoppingCriterionReached(), the tree is destroyed in tree_->cleanTree().
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
      *! Having generated cuts, allows user to think. // amk: added in Cbc 2.10
      generatedCuts,
      *! End of search.
      endSearch
    } ;
 */
CbcEventHandler::CbcAction
VPCEventHandler::event(CbcEvent whichEvent) {
//  printf("\n####################### WHERE FROM: %d #######\n", (int) whichEvent);
  if (reachedEnd_ || numLeafNodes_ >= maxNumLeafNodes_) {
    return CbcEventHandler::CbcAction::stop;
  }

  if (whichEvent == CbcEventHandler::treeStatus) {
    // Here we will update statistics for the first node and save the current leaf nodes as currentNodes_
    // Exceptionally, we also catch a case that we incorrectly labeled a node as pruned

    // Update statistics _only_ for the first node
    if (model_->getNodeCount2() == 1) {
      CbcNode* node = model_->tree()->nodePointer(0);
      NodeStatistics currNodeStats;
      // set node stats for root
      setNodeStatistics(currNodeStats, node, model_, stats_, originalLB_,
                        originalUB_, true, foundSolution_);
      stats_.push_back(currNodeStats);
      // increment number of pruned nodes by those implicitly pruned at the root from strong branching
      for (int i = 0; i < currNodeStats.changed_var.size(); i++) {
        if (model_->solver()->isInteger(currNodeStats.changed_var[i])){
          numPrunedNodes_++;
        }
      }
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
    printf("\n## Number nodes on tree: %d. Number leaf nodes: %d. Depth: %d. Node count: %d. ##\n",
        numNodesOnTree_, numLeafNodes_, model_->currentDepth(), model_->getNodeCount2());
#endif

    currentNodes_.clear();
    pruneNode_ = PruneNodeOption::NO;
    foundSolution_ = false;

    // Check if we are done
    if (model_->getCurrentSeconds() >= maxTime_ ||
        (numLeafNodes_ >= maxNumLeafNodes_ && !owner->params.get(PARTIAL_BB_KEEP_PRUNED_NODES)) ||
        (numLeafNodes_ + numPrunedNodes_ >= maxNumLeafNodes_ && owner->params.get(PARTIAL_BB_KEEP_PRUNED_NODES))) {
      // Save information
      const int status = owner->params.get(PARTIAL_BB_KEEP_PRUNED_NODES) ?
          saveInformationWithPrunes() : saveInformation();

      if (status == 0) {
#ifdef TRACE
      printf(
          "\n## Reached limit (leaf nodes = %d with limit = %d; time = %.2f with limit = %.2f). Exiting. ##\n",
          numLeafNodes_, maxNumLeafNodes_,
          model_->getCurrentSeconds(), maxTime_);
#endif
        reachedEnd_ = true;
        return CbcEventHandler::CbcAction::stop;
      }
    }
    
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
            writeErrorToLog(errorstring, owner->params.logfile);
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
      } // if (model_->getNodeCount2() != stats_[stats_.size() - 1].id + 1)
    } // if (!model_->branchingMethod() || !model_->branchingMethod()->chooseMethod())
  } // (whichEvent == treeStatus)
  else if (whichEvent == CbcEventHandler::node) {
    // Update statistics after a node has been branched on
    // Deal with child
    // 2020-07-15: child_ is set to newNode in CbcModel, which *may* be NULL, when anyAction == -2 after chooseBranch is called within CbcModel::doOneNode
    child_ = model_->currentNode(); // currently broken in Cbc-2.10, because currentNode_ is set to NULL within chooseNode
    // 2017-07-09: Changed from isProvenOptimal because in some cases,
    // we hit an iteration limit during hot starting in OsiChooseVariable.cpp:399
    // Might be good to check dual infeasibility too, but these subproblems
    // should not be unbounded...
    // 2017-07-23: NB: solver_feasible = false may be if child is integer-feasible
    // 2020-08-25: and maybe when child is pruned by bound?
    const bool solver_feasible = !(model_->solver()->isProvenPrimalInfeasible());
#ifdef TRACE
    if (model_->solver()->isProvenDualInfeasible()) {
      error_msg(errorstring, "Did not think that dual infeasibility in the model solver could happen. Check this.\n");
      writeErrorToLog(errorstring, owner->params.logfile);
      exit(1);
    }
#endif
    // When a solution is found, newNode is not null, but it has no branching object
    bool solution_found =
        (solver_feasible && (pruneNode_ == PruneNodeOption::NO) &&
        (child_ && child_->numberUnsatisfied() == 0));
    // The node will be pruned if it exists but pruneNode_ > 0 or a solution is found
    const bool prune_node = (pruneNode_ != PruneNodeOption::NO || solution_found);
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
        writeErrorToLog(errorstring, owner->params.logfile);
        exit(1);
      }
      NodeStatistics currNodeStats;
      // set node stats for child
      setNodeStatistics(currNodeStats, child_, model_, stats_, originalLB_,
                        originalUB_, true, foundSolution_);
      stats_.push_back(currNodeStats);

      // increment number of implicitly pruned nodes from strong branching
      const int p_id = currNodeStats.parent_id;
      const int p_orig_id = stats_[p_id].orig_id;
      std::vector<int> p_var = stats_[p_orig_id].changed_var;
      std::vector<int> p_bound = stats_[p_orig_id].changed_bound;
      std::vector<double> p_value = stats_[p_orig_id].changed_value;
      p_var.push_back(stats_[p_id].variable);
      p_bound.push_back(stats_[p_id].way <= 0);
      p_value.push_back(stats_[p_id].way <= 0 ? stats_[p_id].ub : stats_[p_id].lb);
      std::vector<int> idx_diff = findIndicesOfDifference(
          currNodeStats.changed_var, p_var, currNodeStats.changed_bound, p_bound,
          currNodeStats.changed_value, p_value);
      for (int i = 0; i < idx_diff.size(); i++) {
        // ignore tightened bounds on continuous variables
        if (model_->solver()->isInteger(currNodeStats.changed_var[idx_diff[i]])){
          numPrunedNodes_++;
        }
      }

      parentInfo_ = child_->nodeInfo()->parent(); // 2020-07-21: same as model_->parentNode()->nodeInfo() in trunk code
      if (currNodeStats.number >= numNodes_) {
        numNodes_ = currNodeStats.number + 1;
      }
    } // node is added to the tree
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
#ifdef CBC_VERSION_210PLUS
          assert(parentInfo_ == model_->parentNode()->nodeInfo());
#endif
          break;
        }
      } // find the parent

#ifdef TRACE
      if (!parentInfo_) {
        error_msg(errorstring, "Could not find parent node.\n");
        writeErrorToLog(errorstring, owner->params.logfile);
        exit(1);
      }
#endif

      // I am not sure we are capturing the situation when newNode is active but branch_ is null
      // This should correspond to the node being integer feasible
      // Perhaps this happens when we do not use strong branching...
      if (solver_feasible && !solution_found && (pruneNode_ == PruneNodeOption::NO) && child_ && child_->active() && child_->branchingObject()) {
        error_msg(errorstring, "We should not get here, I think. Is it an integer-feasible solution? If so, it should be saved.\n");
        writeErrorToLog(errorstring, owner->params.logfile);
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
      if (pruneNode_ != PruneNodeOption::PRUNE_BY_INTEGRALITY && !solution_found && !solver_feasible) {
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
            if (!isZero(infeasibility, owner->params.get(doubleParam::EPS))) {
              solution_found = false;
              break;
            }
          }
        }
      } // check if a solution found after all, when solver is "infeasible"

      // Ensure the solution is really integer-feasible and not cut off
      double objectiveValue = model_->solver()->getObjValue(); // sometimes this does not have the right value stored
      if (pruneNode_ == PruneNodeOption::PRUNE_BY_INTEGRALITY || (solution_found &&
          !greaterThanVal(objectiveValue, model_->getCutoff() + model_->getCutoffIncrement()))) {
        const double* sol_ptr = (pruneNode_ == PruneNodeOption::PRUNE_BY_INTEGRALITY) ? savedSolution_.data() : model_->solver()->getColSolution();
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
        } // check col bounds

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
            writeErrorToLog(errorstring, owner->params.logfile);
            exit(1);
          }
          double rhs = model_->solver()->getRightHandSide()[row_ind];
          violation = std::abs(rhs - activity);
          sol.push_back(violation);
          const char sense = model_->solver()->getRowSense()[row_ind];
          if (sense == 'G') {
            if (lessThanVal(activity, rhs, owner->params.get(doubleParam::EPS))) {
              error = row_ind + model_->getNumCols();
              break;
            }
          } else if (sense == 'L') {
            if (greaterThanVal(activity, rhs, owner->params.get(doubleParam::EPS))) {
              error = row_ind + model_->getNumCols();
              break;
            }
          } else if (sense == 'E') {
            if (!isVal(activity, rhs, owner->params.get(doubleParam::EPS))) {
              error = row_ind + model_->getNumCols();
              break;
            }
          } else { // ranged; not 'N'
            // Need to check upper and lower bounds
            // rhs returns rowupper() in this case
            if (greaterThanVal(activity, rhs, owner->params.get(doubleParam::EPS))) {
              error = row_ind + model_->getNumCols();
              break;
            }
            rhs = model_->solver()->getColLower()[row_ind];
            if (lessThanVal(activity, rhs, owner->params.get(doubleParam::EPS))) {
              error = row_ind + model_->getNumCols();
              break;
            }
          }
        } // check row bounds

        if (error >= 0) {
          error_msg(errorstring, "An integer-feasible solution was supposedly found, but it is infeasible somehow. Variable: %d. Violation: %f.\n", error, violation);
          writeErrorToLog(errorstring, owner->params.logfile);
          exit(1);
        }

        // Save the integer feasible solution if it cannot pruned by bound
        // We're being a little extra safe here by allowing an extra epsilon
        double objOffset = 0.;
        model_->solver()->getDblParam(OsiDblParam::OsiObjOffset, objOffset);
        objectiveValue = dotProduct(sol.data(), model_->getObjCoefficients(), model_->getNumCols()) - objOffset;
        if (pruneNode_ == PruneNodeOption::PRUNE_BY_INTEGRALITY && !isVal(objectiveValue, model_->getObjValue())) {
          error_msg(errorstring, "We are assuming that if pruning by integrality, the objective we calculated should be the same as the objective of the solver.\n");
          writeErrorToLog(errorstring, owner->params.logfile);
          exit(1);
        }
        // Code below is wrong, I think (2020-07-21)
        //if (pruneNode_ == PruneNodeOption::PRUNE_BY_INTEGRALITY || (!greaterThanVal(objectiveValue, model_->getCutoff() + model_->getCutoffIncrement()))) {
        //  savedSolution_ = sol;
        //  this->obj_ = objectiveValue;
        //}
      } /* is an integer-feasible solution found? */

      // Also update the infeasible nodes stats
      NodeStatistics prunedNodeStats;
#ifdef CBC_VERSION_210PLUS
      const CbcNode* parent_node = model_->parentNode();
#else
      const CbcNode* parent_node = parentInfo_->owner();
#endif
      if (!parent_node) {
        // 2020-07-21: We will get here, for example, with rlp2_presolved -d 4
        // The owner is the first node in this case, but I am not sure how to get to it
        error_msg(errorstring, "Owner CbcNode of parentInfo_ not found!\n");
        writeErrorToLog(errorstring, owner->params.logfile);
        exit(1);
      }

      // If no solution was found, then we are pruning by bound or by infeasibility
      // The issue is that isProvenPrimalInfeasible is not reliable
      // TODO decide by bound or infeasibility

      setPrunedNodeStatistics(prunedNodeStats, parent_node, model_,
          solution_found || pruneNode_ == PruneNodeOption::PRUNE_BY_INTEGRALITY,
          originalLB_, originalUB_, objectiveValue, foundSolution_, pruneNode_);
      if (pruneNode_ != PruneNodeOption::NO && model_->branchingMethod() && model_->branchingMethod()->chooseMethod()) {
        // I think this means something happened during strong branching, usually an integer-feasible solution found
        prunedNodeStats.variable = model_->branchingMethod()->chooseMethod()->bestObjectIndex();
        prunedNodeStats.way = model_->branchingMethod()->chooseMethod()->bestWhichWay();
      }
      pruned_stats_.push_back(prunedNodeStats);
      numPrunedNodes_++;
      if (prunedNodeStats.number >= numNodes_) {
        numNodes_ = prunedNodeStats.number + 1;
      }
    } // node is pruned

    // Update parent if it has branches left
    if (parentInfo_->numberBranchesLeft() > 0) {
      NodeStatistics currNodeStats;
#ifdef CBC_VERSION_210PLUS
      const CbcNode* parent_node = model_->parentNode();
#else
      const CbcNode* parent_node = parentInfo_->owner();
#endif
      // set node stats for sibling node
      setNodeStatistics(currNodeStats, parent_node, model_, stats_,
          originalLB_, originalUB_, will_create_child);
      stats_.push_back(currNodeStats);
    }
  } // (whichEvent == node)
  else if (whichEvent == CbcEventHandler::beforeSolution2) {
    // This event is called from CbcModel::setBestSolution
    //
    // This event may cause newNode to be deleted, so we should remember not to access deleted memory
    // This action (of deleting) is happening in CbcModel.cpp::chooseBranch
    // in the latest Cbc, anyAction is set at ~line 15357:
    //     anyAction = newNode->chooseOsiBranch(this, oldNode, &usefulInfo, branchingState)
    // after which currentNode_ is set to NULL (hence we need the custom parentNode_ code enabled by SAVE_NODE_INFO)
    // and anyAction is set to -2 if feasible != true
    // If the return code is -2, the node will be deleted, so that is what we are trying to catch

    // If the obj value of the new solution (stored at bestSolution()) is < cutoff,
    // the cutoff will be updated and pruning will happen
    // (cutoff = best integer objective value - 0.999 or so;
    // used to cut any nodes with a worse objective value than this cutoff,
    // or integer solutions with worse value)
    foundSolution_ = true;
    child_ = model_->currentNode(); // this can be NULL if we are in a strong branching loop

    if (!child_) {
      // child should exist...
      error_msg(errorstring,
          "child_ reference should not be NULL within beforeSolution2. model_->getNodeCount2() = %d.\n",
          model_->getNodeCount2());
      writeErrorToLog(errorstring, owner->params.logfile);
      exit(1);
    }

    double objectiveValue = model_->savedSolutionObjective(0); // eventHandler replaces bestObjective_ temporarily, so we are accessing the solution's objective here
    double cutoff = model_->getCutoff(); // this is the same as the previous bestObjective_
    double cutoffIncrement = model_->getCutoffIncrement();
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
      pruneNode_ = PruneNodeOption::PRUNE_BY_INFEASIBILITY;
    } // mip infeasible
    else if (child_->objectiveValue() >= model_->getCutoff()) {
      // CbcModel.cpp:15283
      // if (newNode->objectiveValue() >= getCutoff()) { anyAction = -2; }
      pruneNode_ = PruneNodeOption::PRUNE_BY_BOUND;
    } // pruned by objective
    else if (objectiveValue < cutoff && child_->objectiveValue() >= newCutoff) {
      // Check whether the new solution gives a strictly better bound on the objective
      // and whether this new value implies that the new node will be cut off
      // It is when CbcModel.cpp:13115 is true, so that cutoff is updated, then -2 returned at CbcModel.cpp:15283; e.g., bm23 -70
      // 2017-08-11: Can reach here from CbcModel.cpp:16844, e.g., with bm23 -8 after cuts, from CbcNode::chooseOsiBranch finding a feasible solution during strong branching, so solver never gets the good solution (it is only in choose)
      pruneNode_ = PruneNodeOption::PRUNE_BY_INTEGRALITY;
      savedSolution_.assign(model_->bestSolution(), model_->bestSolution() + model_->getNumCols());
      this->obj_ = objectiveValue;
    } // mip feasible
    else if (!child_->branchingObject() && model_->branchingMethod() && model_->branchingMethod()->chooseMethod() && model_->branchingMethod()->chooseMethod()->goodSolution()) {
      // 2020-07-15 #c630ed3:
      // Catch when an integer solution is found through strong branching (goodSolution() exists),
      // causing pruning by integrality (CbcNode::branch_ is not created inside of chooseOsiBranch())
      // but also a split is found in which both sides have dual bound worse than the integer-feasible solution
      // (so we have solved the problem)
      // This occurs in rlp2_presolved -4 for example, with a tree that is
      //       node 0
      //       obj val 14
      //       var 13
      //       /    |
      // node 2    node 1 (will be branched on)
      // obj 19    obj 17
      //           var 9
      //          /    |
      //        both nodes pruned by integrality with value 19 (node 1 is NOT integral)
      pruneNode_ = PruneNodeOption::PRUNE_BY_INTEGRALITY;
      const double* sol = model_->branchingMethod()->chooseMethod()->goodSolution();
      const double objval = model_->branchingMethod()->chooseMethod()->goodObjectiveValue();
      if (objval < cutoff) {
        savedSolution_.assign(sol, sol + model_->getNumCols());
        this->obj_ = objval;
      }
      /*{
        // 2020-08-23 amk:
        //  If no more nodes will be explored, then we will never encounter another treeStatus event
        //  Keep this node in pruned_stats_
        //  but be sure that this is not going to be added twice
        //  TODO: not sure how to capture the fixed variables here
        NodeStatistics prunedNodeStats;
        const CbcNode* parent_node = NULL;
#ifdef CBC_VERSION_210PLUS
        parent_node = model_->parentNode();
#endif
        if (!parent_node) {
          error_msg(errorstring, "Owner CbcNode of child_ not found!\n");
          writeErrorToLog(errorstring, owner->params.logfile);
          exit(1);
        }
        setPrunedNodeStatistics(prunedNodeStats, parent_node, model_,
            true, originalLB_, originalUB_, objectiveValue, foundSolution_);
        prunedNodeStats.variable = model_->branchingMethod()->chooseMethod()->bestObjectIndex();
        prunedNodeStats.way = model_->branchingMethod()->chooseMethod()->bestWhichWay();
        prunedNodeStats.branch_index = 0;
        pruned_stats_.push_back(prunedNodeStats);
      }*/
    } // no branching object; feasible solution found
  } // (whichEvent == beforeSolution2)
#ifdef CBC_VERSION_210PLUS
  else if (whichEvent == CbcEventHandler::generatedCuts) {
    // Maybe we want to generate cuts without exiting branching
    const bool GEN_CUTS_WHILE_BRANCHING = false;
    if (GEN_CUTS_WHILE_BRANCHING) {
      // Temporarily replace the disjunction associated with this event handler
      PartialBBDisjunction* old = this->owner;
      //std::shared_ptr<PartialBBDisjunction> disj = std::make_shared<PartialBBDisjunction>(old);
      PartialBBDisjunction* disj = new PartialBBDisjunction;
      this->owner = disj;
      saveInformation(); // save necessary info from current leaf nodes into disj

      // Generate cuts from current branch-and-bound leaf nodes
      OsiCuts cuts;
      CglVPC gen(old->params);
      gen.setDisjunction(disj, false);
      gen.generateCuts(*model_->solver(), cuts);
      this->numCuts_ += cuts.sizeCuts();
      //OsiCuts* theseCuts = dynamic_cast<OsiCuts*>(model_->getApplicationData());
      void* appData = model_->getApplicationData();
      OsiCuts* theseCuts = static_cast<OsiCuts*>(appData);
      theseCuts->insert(cuts);
      this->cuts_->insert(cuts);

      // Reset everything to as before
      if (disj) {
        delete disj;
      }
      this->owner = old;
    }
  } // (whichEvent == generatedCuts)
#endif
  else if (whichEvent == CbcEventHandler::endSearch) {
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
  } // (whichEvent == endSearch)

//  printNodeStatistics(stats_);
  return CbcEventHandler::CbcAction::noAction;
} /* VPCEventHandler::event */

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
    writeErrorToLog(errorstring, owner->params.logfile);
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

/****************** PROTECTED **********************/

void VPCEventHandler::initialize(const VPCEventHandler* const source) {
  if (source) {
    this->owner = source->owner;
    this->maxNumLeafNodes_ = source->maxNumLeafNodes_;
    this->numCuts_ = source->numCuts_;
    this->maxTime_ = source->maxTime_;
    this->numNodesOnTree_ = source->numNodesOnTree_;
    this->numLeafNodes_ = source->numLeafNodes_;
    this->numNodes_ = source->numNodes_;
    this->obj_ = source->obj_;
    this->originalSolver_ = dynamic_cast<SolverInterface*>(source->originalSolver_->clone());
    this->originalLB_ = source->originalLB_;
    this->originalUB_ = source->originalUB_;
    this->stats_ = source->stats_;
    this->pruned_stats_ = source->pruned_stats_;
    this->finalNodeIndices_ = source->finalNodeIndices_;
    this->savedSolution_ = source->savedSolution_;
    this->cuts_ = std::make_unique<OsiCuts>(*source->cuts_);
    // Temporary members
    this->currentNodes_ = source->currentNodes_;
    this->parentInfo_ = source->parentInfo_;
    this->child_ = source->child_;
    this->pruneNode_ = source->pruneNode_;
    this->reachedEnd_ = source->reachedEnd_;
    this->foundSolution_ = source->foundSolution_;
    this->checked_nodes_ = source->checked_nodes_;
    this->sorted_nodes_ = source->sorted_nodes_;
    this->numPrunedNodes_ = source->numPrunedNodes_;
  } else {
    this->owner = NULL;
    this->maxNumLeafNodes_ = 0;
    this->numCuts_ = 0;
    this->maxTime_ = 0;
    this->numNodesOnTree_ = 0;
    this->numLeafNodes_ = 0;
    this->numNodes_ = 0;
    this->obj_ = std::numeric_limits<double>::max();
    this->originalSolver_ = NULL;
    this->originalLB_.resize(0);
    this->originalUB_.resize(0);
    this->stats_.resize(0);
    this->pruned_stats_.resize(0);
    this->finalNodeIndices_.resize(0);
    this->savedSolution_.resize(0);
    this->cuts_ = std::make_unique<OsiCuts>();
    // Temporary members
    this->currentNodes_.resize(0);
    this->parentInfo_ = NULL;
    this->child_ = NULL;
    this->pruneNode_ = PruneNodeOption::NO;
    this->reachedEnd_ = false;
    this->foundSolution_ = false;
    this->checked_nodes_ = std::set<int>();
    this->sorted_nodes_ = std::set<int>();
    this->numPrunedNodes_ = 0;
  }
} /* initialize */


/**
 * Update the dual bound of the disjunction if the term worsens it
 *
 * @return void
 */
void VPCEventHandler::updateDualBound(const int orig_node_id,
                                      const DisjunctiveTerm& term){
  // Update the root node information
  // First, identify if this node is on the down or up part of the root split
  // stats_[0].way tells us which is the first child of the root node
  // The other node with 0 as a parent is therefore the second child
  bool is_down = false, is_up = false;
  int curr_id = orig_node_id;
  const int way = stats_[0].way;
  while (!is_down && !is_up) {
    if (stats_[curr_id].parent_id == -1) {
      if (curr_id == 0) { // check if this is the first child of the root
        if (way <= 0)
          is_down = true;
        else
          is_up = true;
      } else { // if not the first child, then it is the second child
        if (way <= 0)
          is_up = true;
        else
          is_down = true;
      }
    } else {
      curr_id = stats_[curr_id].parent_id;
    }
  }
  // Now update the bound if possible
  double& bound = is_down ? this->owner->root.boundD : this->owner->root.boundU;
  if (lessThanVal(term.obj, bound)) {
    bound = term.obj;
  }
}

bool VPCEventHandler::setupDisjunctiveTerm(const int node_id,
    const int branching_variable, const int branching_way,
    const double branching_value, const SolverInterface* const tmpSolverBase,
    const int curr_num_changed_bounds,
    std::vector<std::vector<int> >& commonTermIndices,
    std::vector<std::vector<double> >& commonTermCoeff,
    std::vector<double>& commonTermRHS) {
  bool isFeasible = false;
  const int ind[] = {branching_variable};
  const double coeff[] = {(branching_way <= 0) ? -1. : 1.};
  const double rhs = coeff[0] * branching_value;
  const int orig_node_id = stats_[node_id].orig_id;

  SolverInterface* tmpSolver = dynamic_cast<SolverInterface*>(tmpSolverBase->clone());
  if (branching_way <= 0)
    tmpSolver->setColUpper(branching_variable, branching_value);
  else
    tmpSolver->setColLower(branching_variable, branching_value);
  //tmpSolver->addRow(1, ind, coeff, rhs, tmpSolver->getInfinity());
  tmpSolver->resolve();
  if (checkSolverOptimality(tmpSolver, true)) {
    enableFactorization(tmpSolver, owner->params.get(doubleParam::EPS)); // this may change the solution slightly
    this->owner->num_terms++;
    Disjunction::setCgsName(this->owner->name, 1, ind, coeff, rhs);
    if (curr_num_changed_bounds > 0)
      Disjunction::setCgsName(this->owner->name, curr_num_changed_bounds,
          commonTermIndices, commonTermCoeff, commonTermRHS, true);
    DisjunctiveTerm term;
    term.basis = dynamic_cast<CoinWarmStartBasis*>(tmpSolver->getWarmStart());
    term.obj = tmpSolver->getObjValue();
    term.changed_var = stats_[orig_node_id].changed_var;
    term.changed_var.push_back(branching_variable);
    term.changed_bound = stats_[orig_node_id].changed_bound;
    term.changed_bound.push_back(branching_way == 1 ? 0 : 1);
    term.changed_value = stats_[orig_node_id].changed_value;
    term.changed_value.push_back(branching_value);
    owner->terms.push_back(term);
    isFeasible = true;

    updateDualBound(orig_node_id, term);
  } else {
    this->numLeafNodes_--;
  }
  if (tmpSolver)
    delete tmpSolver;
  return isFeasible;
} /* setupDisjunctiveTerm */

/**
 * Set the warm start information for the parent node
 *
 * @param solver: the solver to set the warm start information for
 * @param basis: the warm start information to set
 * @param tmp_ind: the index of the parent node amongst all parent nodes being iterated over
 * @return void
 */
void VPCEventHandler::setWarmStart(
    SolverInterface* solver, CoinWarmStartBasis* basis, int tmp_ind, int node_id,
    const int curr_num_changed_bounds, const std::vector<std::vector<int> >& commonTermIndices,
    const std::vector<std::vector<double> >& commonTermCoeff, const std::vector<double>& commonTermRHS){

  // Set the parent node warm start
  if (!(solver->setWarmStart(basis))) {
    error_msg(errorstring,
              "Warm start information not accepted for parent node %d.\n", tmp_ind);
    writeErrorToLog(errorstring, owner->params.logfile);
    exit(1);
  }

  // Resolve
#ifdef TRACE
  printf("\n## Solving for parent node %d/%d. ##\n", tmp_ind + 1, numNodesOnTree_);
#endif
  solver->resolve();
  if (!checkSolverOptimality(solver, false)) {
    error_msg(errorstring, "Solver not proven optimal for node %d.\n",
              node_id);
    writeErrorToLog(errorstring, owner->params.logfile);
    exit(1);
  }
  // If called from saveInformationWithPrunes, we can't guarantee we get the same
  // objective for the node that generated the basis as we may have reset continuous
  // variable bounds. Also node 0 shouldn't match original solver LP objective
  // because no bounds have been tightened
  if (tmp_ind < 0 && node_id < 0) {
    return;
  }
  // Sometimes we run into a few issues getting the "right" value
  double ratio = solver->getObjValue() / stats_[node_id].obj;
  if (ratio < 1.) {
    ratio = 1. / ratio;
  }
  if (greaterThanVal(ratio, 1.03)) {
    solver->resolve();
    ratio = solver->getObjValue() / stats_[node_id].obj;
    if (ratio < 1.) {
      ratio = 1. / ratio;
    }
  }

#ifdef TRACE
  double curr_nb_obj_val = solver->getObjValue() - originalSolver_->getObjValue();
    printf("DEBUG: Node: %d .......... Obj val: %.3f .......... NB obj val: %.3f\n", node_id,
        solver->getObjValue(), curr_nb_obj_val);
    printNodeStatistics(stats_[node_id], true);
#endif
  if (!isVal(solver->getObjValue(), stats_[node_id].obj,
             owner->params.get(doubleConst::DIFFEPS))) {
#ifdef TRACE
    std::string commonName;
      Disjunction::setCgsName(commonName, curr_num_changed_bounds, commonTermIndices,
          commonTermCoeff, commonTermRHS, false);
      printf("Bounds changed: %s.\n", commonName.c_str());
#endif
    // Allow it to be up to 3% off without causing an error
    if (greaterThanVal(ratio, 1.03)) {
      error_msg(errorstring,
                "Objective at parent node %d/%d (node id %d) is incorrect. During BB, it was %s, now it is %s.\n",
                tmp_ind + 1, this->numNodesOnTree_, node_id,
                stringValue(stats_[node_id].obj, "%1.3f").c_str(),
                stringValue(solver->getObjValue(), "%1.3f").c_str());
      writeErrorToLog(errorstring, owner->params.logfile);
      exit(1);
    } else {
      warning_msg(warnstring,
                  "Objective at parent node %d/%d (node id %d) is somewhat incorrect. During BB, it was %s, now it is %s.\n",
                  tmp_ind + 1, this->numNodesOnTree_, node_id,
                  stringValue(stats_[node_id].obj, "%1.3f").c_str(),
                  stringValue(solver->getObjValue(), "%1.3f").c_str());
    }
  } // check that the parent node objective matches
}

/** 
 * Reset before saveInformation
 */
void VPCEventHandler::clearInformation() {
  owner->setupAsNew();
  checked_nodes_.clear();
  sorted_nodes_.clear();
} /* clearInformation */

/**
 * Save all relevant information before it gets deleted by BB finishing
 *
 * @return status: 0 if everything is okay, otherwise 1 (e.g., if all but one of the terms is infeasible)
 */
int VPCEventHandler::saveInformation() {
  int status = 0;
  const bool hitTimeLimit = model_->getCurrentSeconds() >= maxTime_;
  const bool hitHardNodeLimit = false;
  //const bool hitHardNodeLimit = model_->getNodeCount2() > maxNumLeafNodes_ * 10;

  clearInformation();
  this->owner->data.num_nodes_on_tree = this->getNumNodesOnTree();
  this->owner->data.num_partial_bb_nodes = model_->getNodeCount(); // save number of nodes looked at
  this->owner->data.num_pruned_nodes = this->getPrunedStatsVector().size() - this->isIntegerSolutionFound();
  this->owner->data.num_fixed_vars = model_->strongInfo()[1]; // number fixed during b&b

  // Save variables with bounds that were changed at the root
  const int num_stats = this->stats_.size();
  if (num_stats > 0) {
    this->owner->common_changed_bound = this->stats_[0].changed_bound;
    this->owner->common_changed_value = this->stats_[0].changed_value;
    this->owner->common_changed_var = this->stats_[0].changed_var;
    this->owner->root_obj = this->stats_[0].obj;
    this->owner->root.var = this->stats_[0].variable;
    this->owner->root.val = this->stats_[0].value;
  }

  // If an integer solution was found, save it
  if (isIntegerSolutionFound()) { // todo: abstract for adding all integer solutions
    // 2017-08-11: Only one integer-feasible solution needs to be kept (the best one),
    // So we drop old code that would prune integer feasible solutions that are above cutoff
    this->owner->num_terms++;
    this->owner->integer_sol = savedSolution_;
    std::string feasname = "feasSol0";
    Disjunction::setCgsName(this->owner->name, feasname);
  }

  this->owner->terms.reserve(2 * numNodesOnTree_);

  // Set up original basis including bounds changed at root
  CoinWarmStartBasis* original_basis = dynamic_cast<CoinWarmStartBasis*>(originalSolver_->getWarmStart());

  // For each node on the tree, use the warm start basis and branching directions to save the disjunctive node
  for (int tmp_ind = 0; tmp_ind < numNodesOnTree_; tmp_ind++) {
    // Exit early if needed
    if (!hitTimeLimit && !hitHardNodeLimit 
        && (this->owner->num_terms + 2 * (numNodesOnTree_ - tmp_ind) <= 0.5 * maxNumLeafNodes_)) {
      status = 1;
      break;
    }

    CoinWarmStartBasis* parent_basis = dynamic_cast<CoinWarmStartBasis*>(original_basis->clone());
    CbcNode* node = model_->tree()->nodePointer(tmp_ind);
    CbcNodeInfo* nodeInfo = node->nodeInfo();
    nodeInfo->buildRowBasis(*parent_basis);
    finalNodeIndices_.push_back(stats_[nodeInfo->nodeNumber()].id);

    const int node_id = stats_[nodeInfo->nodeNumber()].id;
    // the node's id the first time we saw it (this differs from node_id if one
    // branching direction of the node has already been explored)
    const int orig_node_id = stats_[node_id].orig_id;

    SolverInterface* tmpSolverBase =
        dynamic_cast<SolverInterface*>(originalSolver_->clone());
    setLPSolverParameters(tmpSolverBase, owner->params.get(VERBOSITY), owner->params.get(TIMELIMIT));
//    tmpSolverBase->disableFactorization(); // seg fault

    // Change bounds in the solver
    // First the root node bounds
    const int init_num_changed_bounds = this->owner->common_changed_var.size();
    for (int i = 0; i < init_num_changed_bounds; i++) {
      const int col = this->owner->common_changed_var[i];
      const double val = this->owner->common_changed_value[i];
      if (this->owner->common_changed_bound[i] <= 0) {
        tmpSolverBase->setColLower(col, val);
      } else {
        tmpSolverBase->setColUpper(col, val);
      }
    }

    // Collect and apply the parent changed node bounds
    const int parent_num_changed_bounds =
        (orig_node_id == 0) ? 0 : stats_[orig_node_id].changed_var.size();
    std::vector < std::vector<int> > parentTermIndices(parent_num_changed_bounds);
    std::vector < std::vector<double> > parentTermCoeff(parent_num_changed_bounds);
    std::vector<double> parentTermRHS(parent_num_changed_bounds);
    for (int i = 0; i < parent_num_changed_bounds; i++) {
      const int col = stats_[orig_node_id].changed_var[i];
      const double coeff = (stats_[orig_node_id].changed_bound[i] <= 0) ? 1. : -1.;
      const double val = stats_[orig_node_id].changed_value[i];
      parentTermIndices[i].resize(1, col);
      parentTermCoeff[i].resize(1, coeff);
      parentTermRHS[i] = coeff * val;
      if (stats_[orig_node_id].changed_bound[i] <= 0) {
        tmpSolverBase->setColLower(col, val);
      } else {
        tmpSolverBase->setColUpper(col, val);
      }
    }

    // Set the parent node warm start and make sure we're kosher
    // start from the root node LP relaxation optimal basis to speed things up
    setWarmStart(tmpSolverBase, parent_basis, tmp_ind, node_id, parent_num_changed_bounds,
                 parentTermIndices, parentTermCoeff, parentTermRHS);

    // Now we check the branch or branches left for this node
    bool hasFeasibleChild = false;
    const int branching_index = stats_[node_id].branch_index;
    const int branching_variable = stats_[node_id].variable;
    int branching_way = stats_[node_id].way;
    double branching_value = (branching_way <= 0) ? stats_[node_id].ub : stats_[node_id].lb;

#ifdef TRACE
    printf("\n## Solving first child of parent %d/%d for term %d. ##\n",
        tmp_ind + 1, this->numNodesOnTree_, this->owner->num_terms);
#endif
    hasFeasibleChild = setupDisjunctiveTerm(node_id, branching_variable,
        branching_way, branching_value, tmpSolverBase, parent_num_changed_bounds,
        parentTermIndices, parentTermCoeff, parentTermRHS);

    // Check if we exit early
    if (!hitTimeLimit && !hitHardNodeLimit 
        && (this->owner->num_terms + 2 * (numNodesOnTree_ - tmp_ind) - 1 <= 0.5 * maxNumLeafNodes_)) {
      status = 1;
    }

    if (status == 0 && branching_index == 0) { // should we compute a second branch?
#ifdef TRACE
      printf("\n## Solving second child of parent %d/%d for term %d. ##\n",
          tmp_ind + 1, this->numNodesOnTree_, this->owner->num_terms);
#endif
      branching_way = (stats_[node_id].way <= 0) ? 1 : 0;
      branching_value =
          (branching_way <= 0) ? branching_value - 1 : branching_value + 1;
      const bool childIsFeasible = setupDisjunctiveTerm(node_id, branching_variable, branching_way,
              branching_value, tmpSolverBase, parent_num_changed_bounds,
              parentTermIndices, parentTermCoeff, parentTermRHS);
      hasFeasibleChild = hasFeasibleChild || childIsFeasible;
    } // end second branch computation

    if (hasFeasibleChild) {
      if (stats_[node_id].depth + 1 < owner->data.min_node_depth) {
        owner->data.min_node_depth = stats_[node_id].depth + 1;
      }
      if (stats_[node_id].depth + 1 > owner->data.max_node_depth) {
        owner->data.max_node_depth = stats_[node_id].depth + 1;
      }
    }

    if (tmpSolverBase)
      delete tmpSolverBase;
    if (parent_basis)
      delete parent_basis;
  } // loop over num nodes on tree // DONE

  if (original_basis) {
    delete original_basis;
  }

  // If number of real (feasible) terms is too few, we should keep going, unless we have been branching for too long
  if (status == 0 && !hitTimeLimit && !hitHardNodeLimit && this->owner->num_terms <= 0.5 * maxNumLeafNodes_) {
    status = 1;
#ifdef TRACE
      warning_msg(warnstring,
          "Save information exiting with status %d because number of feasible terms is calulated to be %d "
          "(whereas requested was %d).#\n",
          status, this->owner->num_terms, maxNumLeafNodes_);
  #endif
  }
  return status;
} /* saveInformation */

/**
 * @details drop all tightened bounds occurring on continuous variables. We can't
 * create disjunctive terms for these variables since we assume disjunctive
 * constraints are placed on integer variables.
 */
void VPCEventHandler::removeContinuousVariableTightenedBounds(std::vector<NodeStatistics>& stat_vec){

  // iterate over each node to remove tightened continuous bounds from our disjunction
  for (int stat_idx = 0; stat_idx < stat_vec.size(); stat_idx++){

    int branch_var = stat_vec[stat_idx].variable;
    verify(stat_vec[stat_idx].found_integer_solution || branch_var < 0 ||
           originalSolver_->isInteger(branch_var), "Branch variables must be integer if provided.");

    // new vectors to store the tightened bounds
    std::vector<int> var;
    std::vector<int> bound;
    std::vector<double> value;

    // iterate over tightened variables
    for (int var_idx = 0; var_idx < stat_vec[stat_idx].changed_var.size(); var_idx++){

      // only keep tightened variables that are integer
      if (originalSolver_->isInteger(stat_vec[stat_idx].changed_var[var_idx])){
        var.push_back(stat_vec[stat_idx].changed_var[var_idx]);
        bound.push_back(stat_vec[stat_idx].changed_bound[var_idx]);
        value.push_back(stat_vec[stat_idx].changed_value[var_idx]);
      }
    }

    // update the stat_vec with the tightened bounds
    stat_vec[stat_idx].changed_var = var;
    stat_vec[stat_idx].changed_bound = bound;
    stat_vec[stat_idx].changed_value = value;
  }
} /* removeContinuousVariableTightenedBounds */


/**
 * @details reorder the branching decisions leading to each node in the tree in
 * the order that they occurred
 */
void VPCEventHandler::sortBranchingDecisions(const int node_id) {

  verify(0 <= node_id && node_id < stats_.size(), "node_id must be in stats_.");

  // declare variables - need orig_id's for the branching decisions
  int orig_node_id = stats_[node_id].orig_id; // to get tightenings at this node
  // both siblings have same parent so querying with orig_id is fine
  // we sometimes get the wrong node for parent id if we use node_id because parent_id
  // reflects the state of the parent node when this node was created and the parent may
  // have been branched on already leading to a wrong id for our purposes
  int parent_id = stats_[orig_node_id].parent_id; // to get parent branching decision
  int orig_parent_id = parent_id >= 0 ? stats_[parent_id].orig_id : -1; // to get tightenings at the parent

  // if the node has not already been checked
  if (sorted_nodes_.find(orig_node_id) == sorted_nodes_.end()) {

    // if the unchecked node has a parent that isn't a root child, recurse to it
    if (orig_parent_id > 0){
      sortBranchingDecisions(parent_id);
    }

    // start with the branching decisions decided before the parent node branched
    std::vector<int> parent_var;
    std::vector<int> parent_bound;
    std::vector<double> parent_value;
    // ignore the branching decisions prior to branching on children of root since
    // they're stored in common_x
    if (orig_parent_id > 0){
      parent_var = stats_[orig_parent_id].changed_var;
      parent_bound = stats_[orig_parent_id].changed_bound;
      parent_value = stats_[orig_parent_id].changed_value;
    }

    // follow up with parent branching choice
    if (orig_parent_id >= 0){
      int bound = stats_[parent_id].way == 1 ? 0 : 1;
      double val = bound == 0 ? stats_[parent_id].lb : stats_[parent_id].ub;
      parent_var.push_back(stats_[parent_id].variable);
      parent_bound.push_back(bound);
      parent_value.push_back(val);
    }

    // end with bound tightening that occurred at this node prior to branching
    std::vector<int> differing_indices =
        findIndicesOfDifference(stats_[orig_node_id].changed_var, parent_var,
                                stats_[orig_node_id].changed_bound, parent_bound,
                                stats_[orig_node_id].changed_value, parent_value);
    std::vector<int> node_exclusive_var =
        subselectVector(stats_[orig_node_id].changed_var, differing_indices);
    std::vector<int> node_exclusive_bound =
        subselectVector(stats_[orig_node_id].changed_bound, differing_indices);
    std::vector<double> node_exclusive_val =
        subselectVector(stats_[orig_node_id].changed_value, differing_indices);

    // update the stats_ with the sorted branching decisions
    stats_[orig_node_id].changed_var = parent_var;
    stats_[orig_node_id].changed_bound = parent_bound;
    stats_[orig_node_id].changed_value = parent_value;
    stats_[orig_node_id].changed_var.insert(stats_[orig_node_id].changed_var.end(),
                                            node_exclusive_var.begin(), node_exclusive_var.end());
    stats_[orig_node_id].changed_bound.insert(stats_[orig_node_id].changed_bound.end(),
                                              node_exclusive_bound.begin(), node_exclusive_bound.end());
    stats_[orig_node_id].changed_value.insert(stats_[orig_node_id].changed_value.end(),
                                              node_exclusive_val.begin(), node_exclusive_val.end());
  }

  // mark this node as checked
  sorted_nodes_.insert(orig_node_id);
} /* orderBranchingDecisions */

/**
 * @details Creates a disjunctive term in a generalized enough manner to allow
 * for pruned terms
 *
 * @param tmpSolverParent: the solver with the parent node's bounds already applied
 * @param term_var: the indices of variables with tightened bounds in tmpSolverParent
 * @param term_bound: the updated bounds for the variables already tightened in tmpSolverParent
 * @param term_val: the values that have already been applied to the tightened variables' bounds
 * @param branching_variable: the index of the variable that was branched on.
 * If -1, then the term did not have a branching decision made (i.e. was pruned)
 * @param branching_bound: the bound that was changed on the branching variable
 * @param branching_value: the value that the branching variable was set to
 * @param term_type: Whether or not this term represents a node "pruned" by branch
 * and bound, a "complement" to a tightened variable between branching, or a "leaf"
 * @param orig_node_id: the id of this term's parent node the first time the parent was processed
 *
 * @return whether or not the term represents a feasible leaf node
 */
bool VPCEventHandler::setupDisjunctiveTerm(
    const SolverInterface* const tmpSolverParent, const std::vector<int>& term_var,
    const std::vector<int>& term_bound, const std::vector<double>& term_val,
    const int branching_variable, const int branching_bound, const double branching_value,
    const std::string& term_type, const int orig_node_id) {

  // we assume tmpSolverParent already has the term_ variable bounds applied
  checkBounds(tmpSolverParent, term_var, term_bound, term_val);
  verify(-1 <= orig_node_id && orig_node_id < static_cast<int>(stats_.size()),
         "node_id must be -1 if not in stats_.");
  verify(term_type == "pruned" || term_type == "complement" || term_type == "leaf",
         "term_type must be pruned, complement, or leaf.");
  verify((0 <= orig_node_id) == (term_type == "leaf"),
         "orig_node_id must be nonnegative if and only if term_type is leaf.");
  if (branching_variable >= 0) {
    checkColumnAndBound(tmpSolverParent, branching_variable, branching_bound);
  }
  
  // create a copy of the solver so we can add the branching decision and recover the basis
  SolverInterface* tmpSolverNode =
      dynamic_cast<SolverInterface*>(tmpSolverParent->clone());
  if (branching_variable >= 0){
    addVarBound(tmpSolverNode, branching_variable, branching_bound, branching_value);
  }
  tmpSolverNode->resolve();
  enableFactorization(tmpSolverNode, owner->params.get(doubleParam::EPS)); // this may change the solution slightly

  // create the term
  this->owner->num_terms++;
  DisjunctiveTerm term;
  term.basis = dynamic_cast<CoinWarmStartBasis*>(tmpSolverNode->getWarmStart());
  term.obj = tmpSolverNode->getObjValue();
  term.is_feasible = checkSolverOptimality(tmpSolverNode, true);
  if (!term.is_feasible){
    std::cerr << "Infeasible term!!!" << std::endl;
  }
  term.type = term_type;
  term.changed_var = term_var;
  term.changed_bound = term_bound;
  term.changed_value = term_val;
  term.changed_var.push_back(branching_variable);
  term.changed_bound.push_back(branching_bound);
  term.changed_value.push_back(branching_value);

  // add the term
  owner->terms.push_back(term);

  // delete the solver
  delete tmpSolverNode;

  if (term.type == "leaf"){
    if (term.is_feasible){
      // if this was a feasible leaf, update the dual bound
      updateDualBound(orig_node_id, term);
      return true;
    } else {
      // decrement the number of leaf nodes if the leaf was infeasible
      this->numLeafNodes_--;
    }
  } else {
    owner->num_pruned_terms++;
  }

  // return whether or not this was a feasible leaf
  return false;
} /* setupDisjunctiveTerm */

/**
 * @details Given the bounds tightened at the child node that occurred due to strong
 * branching, complete the corresponding subtree rooted at the parent node, and
 * create the disjunctive terms for all leaves that are not the child node.
 *
 * @param child_exclusive_var: the indices of variables with bounds tightened by
 * strong branching at the child node
 * @param child_exclusive_bound: the bounds tightened by strong branching at the
 * child node
 * @param child_exclusive_val: the values that variables with bounds tightened by
 * strong branching at the child node are bound by
 * @param tmpSolver: a solver with the parent node's bounds already applied
 * @param parent_var: the indices of variables with bounds tightened at the
 * parent node
 * @param parent_bound: the bounds tightened at the parent node
 * @param parent_val: the values that variables with bounds tightened at the
 * parent node are bound by
 *
 * @return void
 */
void VPCEventHandler::createStrongBranchingTerms(
    std::vector<int> child_exclusive_var, std::vector<int> child_exclusive_bound,
    std::vector<double> child_exclusive_val, const SolverInterface* const tmpSolver,
    std::vector<int> parent_var, std::vector<int> parent_bound,
    std::vector<double> parent_val){

  // we assume tmpSolver already has the parent_ variable bounds applied
  checkBounds(tmpSolver, parent_var, parent_bound, parent_val);
  // we assume tmpSolver should not have any of the bounds exclusive to the child node applied
  checkBounds(tmpSolver, child_exclusive_var, child_exclusive_bound,
              child_exclusive_val, false);

  // create a solver to track the parent node of each infeasible term we create
  SolverInterface* tmpSolverParent =
      dynamic_cast<SolverInterface*>(tmpSolver->clone());

  // create the nodes that were ignored between parent and child due to strong branching
  for (int idx = 0; idx < child_exclusive_var.size(); idx++) {
    // flip the branching decision to get the parameters for creating the pruned term
    int bound = child_exclusive_bound[idx] == 0 ? 1 : 0;
    double val = child_exclusive_val[idx] + (child_exclusive_bound[idx] == 0 ? -1 : 1);
    // create the disjunctive term for the node pruned by strong branching
    setupDisjunctiveTerm(tmpSolverParent, parent_var, parent_bound, parent_val,
                         child_exclusive_var[idx], bound, val, "complement");
    // augment the solver and common branching decisions for the next term
    addVarBound(tmpSolverParent, child_exclusive_var[idx], child_exclusive_bound[idx],
                child_exclusive_val[idx], parent_var, parent_bound, parent_val);
  }
  delete tmpSolverParent;
}

/**
 * @details Recursively find the variable bounds tightened by strong branching between
 * node <node_id> and the root node, and create the corresponding disjunctive terms
 *
 * @param node_id: the id of the node to start from
 * @param common_var: the variables with bounds tightened at the root node
 * @param common_bound: the bounds tightened at the root node
 * @param common_value: the values that variables are bound by at the root node
 * @param tmpSolverBase: the solver to use to create the disjunctive terms. It
 * is assumed to have the common_ variable bounds applied
 *
 * @return void
 */
void VPCEventHandler::recursivelyCreateStrongBranchingTerms(
    const int node_id, const std::vector<int>& common_var,
    const std::vector<int>& common_bound, const std::vector<double>& common_value,
    const SolverInterface* const tmpSolverBase){

  // we assume tmpSolverBase already has the common_ variable bounds applied
  checkBounds(tmpSolverBase, common_var, common_bound, common_value);
  verify(0 <= node_id && node_id < stats_.size(), "node_id must be in stats_.");

  // declare variables - need orig_id's for the branching decisions
  int orig_node_id = stats_[node_id].orig_id; // to get tightenings at this node
  // both siblings have same parent so querying with orig_id is fine
  // we sometimes get the wrong node for parent id if we use node_id because parent_id
  // reflects the state of the parent node when this node was created and the parent may
  // have been branched on already leading to a wrong id for our purposes
  int parent_id = stats_[orig_node_id].parent_id; // to get parent branching decision
  int orig_parent_id = stats_[parent_id].orig_id; // to get tightenings at the parent

  // if the node has a parent and has not already been checked
  if (parent_id >= 0 && checked_nodes_.find(orig_node_id) == checked_nodes_.end()){

    // recurse on the parent node
    recursivelyCreateStrongBranchingTerms(parent_id, common_var, common_bound,
                                          common_value, tmpSolverBase);

    // if the node has more tightened bounds relative to its parent than just its parent's branching choice.
    // first condition needs to ignore child's of root node bounds b/c it includes root node fixes
    // which are excluded from all other nodes and accounted for in common_x vectors
    if ((stats_[orig_node_id].changed_var.size() > 1 && orig_parent_id == 0) ||
        (stats_[orig_node_id].changed_var.size() > stats_[orig_parent_id].changed_var.size() + 1)){
      SolverInterface* tmpSolver = dynamic_cast<SolverInterface*>(tmpSolverBase->clone());

      // create vectors to hold the terms tightened bounds
      std::vector<int> term_bound = common_bound;
      std::vector<double> term_value = common_value;
      std::vector<int> term_var = common_var;

      // augment term_x vector and the solver with parent's bounds prior to when it branched
      // again skip children of root b/c they include root node fixes which were captured in common_x
      if (orig_parent_id > 0){
        for (int idx = 0; idx < stats_[orig_parent_id].changed_var.size(); idx++){
          addVarBound(tmpSolver, stats_[orig_parent_id].changed_var[idx],
                      stats_[orig_parent_id].changed_bound[idx],
                      stats_[orig_parent_id].changed_value[idx],
                      term_var, term_bound, term_value);
        }
      }

      // add parent's bound from branching to each term_x vector and the solver
      addVarBound(tmpSolver, stats_[parent_id].variable, stats_[parent_id].way <= 0 ? 1 : 0,
                  stats_[parent_id].way <= 0 ? stats_[parent_id].ub : stats_[parent_id].lb,
                  term_var, term_bound, term_value);

      // determine which bounds were set by strong branching on the current node
      std::vector<int> differing_indices =
          findIndicesOfDifference(stats_[orig_node_id].changed_var, term_var,
                                  stats_[orig_node_id].changed_bound, term_bound,
                                  stats_[orig_node_id].changed_value, term_value);
      std::vector<int> child_exclusive_var =
          subselectVector(stats_[orig_node_id].changed_var, differing_indices);
      std::vector<int> child_exclusive_bound =
          subselectVector(stats_[orig_node_id].changed_bound, differing_indices);
      std::vector<double> child_exclusive_val =
          subselectVector(stats_[orig_node_id].changed_value, differing_indices);

      // create disjunctive terms for those fixed by strong branching at current node
      createStrongBranchingTerms(child_exclusive_var, child_exclusive_bound,
                                 child_exclusive_val, tmpSolver, term_var,
                                 term_bound, term_value);
    }

    // add the current node to the set of checked nodes
    checked_nodes_.insert(orig_node_id);
  }
}

/**
 * @details Save all relevant information before it gets deleted by BB finishing
 * This version of saveInformation also keeps information related to pruned nodes
 * (including variable fixes made by strong branching)
 *
 * @return status: 0 if everything is okay, otherwise 1 (e.g., if all but one of the terms is infeasible)
 */
int VPCEventHandler::saveInformationWithPrunes() {
  const bool hitTimeLimit = model_->getCurrentSeconds() >= maxTime_;
  const bool hitHardNodeLimit = false;
  //const bool hitHardNodeLimit = model_->getNodeCount2() > maxNumLeafNodes_ * 10;

  clearInformation();
  this->owner->data.num_nodes_on_tree = this->getNumNodesOnTree();
  this->owner->data.num_partial_bb_nodes = model_->getNodeCount(); // save number of nodes looked at
  this->owner->data.num_pruned_nodes = this->getPrunedStatsVector().size() - this->isIntegerSolutionFound();
  this->owner->data.num_fixed_vars = model_->strongInfo()[1]; // number fixed during b&b

  // clear out tightened bounds on continuous variables since we can't create terms with them
  removeContinuousVariableTightenedBounds(stats_);
  removeContinuousVariableTightenedBounds(pruned_stats_);

  // Save variables with bounds that were changed at the root
  // don't save as common to disjunction and instead prepend to each disjunctive term
  std::vector<int> common_bound;
  std::vector<double> common_value;
  std::vector<int> common_var;
  if (this->stats_.size() > 0) {
    common_bound = this->stats_[0].changed_bound;
    common_value = this->stats_[0].changed_value;
    common_var = this->stats_[0].changed_var;
    this->owner->root_obj = this->stats_[0].obj;
    this->owner->root.var = this->stats_[0].variable;
    this->owner->root.val = this->stats_[0].value;
  }

  // create a solver to represent the original LP relaxation
  SolverInterface* tmpSolverRoot =
      dynamic_cast<SolverInterface*>(originalSolver_->clone());
  setLPSolverParameters(tmpSolverRoot, owner->params.get(VERBOSITY), owner->params.get(TIMELIMIT));
  CoinWarmStartBasis* original_basis = dynamic_cast<CoinWarmStartBasis*>(originalSolver_->getWarmStart());
  setWarmStart(tmpSolverRoot, original_basis);

  // add disjunctive terms for the pruned branching decisions found at the root
  createStrongBranchingTerms(common_var, common_bound, common_value, tmpSolverRoot);

  // add disjunctive terms for the nodes pruned by branch and bound
  for (int tmp_ind = 0; tmp_ind < (int) pruned_stats_.size(); tmp_ind++) {

    // create a new solver (solve from LP optimal basis since parent basis may be gone)
    SolverInterface* tmpSolverPruned = dynamic_cast<SolverInterface*>(originalSolver_->clone());
    setWarmStart(tmpSolverPruned, original_basis);

    // Change bounds in the solver for the root node variable fixes
    const int init_num_changed_bounds = common_var.size();
    for (int i = 0; i < init_num_changed_bounds; i++) {
      addVarBound(tmpSolverPruned, common_var[i], common_bound[i], common_value[i]);
    }

    // For safety, rebuild the term from parent's information
    // do I need to fix the id usage here too?
    const int parent_id = pruned_stats_[tmp_ind].parent_id;
    const int node_id = parent_id;
    const int orig_node_id = stats_[node_id].orig_id;
    const int branching_variable = stats_[node_id].variable;
    const int branching_way = (stats_[node_id].way == 1);
    const double branching_value = (branching_way <= 0) ?
        stats_[node_id].ub : stats_[node_id].lb;

    // order branching decisions leading to node so they're ordered in disjunctive terms
    sortBranchingDecisions(node_id);

    // find and create any disjunctive terms pruned due to strong branching between
    // the root and the node we're about to create a disjunctive term for
    recursivelyCreateStrongBranchingTerms(node_id, common_var, common_bound,
                                          common_value, tmpSolverPruned); // stop here on tmp_idx = 3

    // Collect the parent changed node bounds
    const int curr_num_changed_bounds =
        (orig_node_id == 0) ? 0 : stats_[orig_node_id].changed_var.size();
    std::vector<int> term_bound = common_bound;
    std::vector<double> term_value = common_value;
    std::vector<int> term_var = common_var;
    for (int i = 0; i < curr_num_changed_bounds; i++) {
      addVarBound(tmpSolverPruned, stats_[orig_node_id].changed_var[i],
                  stats_[orig_node_id].changed_bound[i],
                  stats_[orig_node_id].changed_value[i], term_var,
                  term_bound, term_value);
    } // loop over changed bounds

    // create the term
    setupDisjunctiveTerm(tmpSolverPruned, term_var, term_bound, term_value,
                         branching_variable, branching_way == 1 ? 0 : 1,
                         branching_value, "pruned");
  } // loop over pruned nodes

  // If an integer solution was found, save it
  if (isIntegerSolutionFound()) {
    this->owner->num_terms++;
    this->owner->integer_sol = savedSolution_;
    std::string feasname = "feasSol0";
    Disjunction::setCgsName(this->owner->name, feasname);
  }

  // For each node on the tree, use the warm start basis and branching directions
  // to create the corresponding disjunctive term
  for (int tmp_ind = 0; tmp_ind < numNodesOnTree_; tmp_ind++) {
    // create a new copy of the solver for the term/s we're about to create
    SolverInterface* tmpSolverBase = dynamic_cast<SolverInterface*>(originalSolver_->clone());
    setLPSolverParameters(tmpSolverBase, owner->params.get(VERBOSITY), owner->params.get(TIMELIMIT));

    // get this node's (the parent of disjunctive term/s about to be made) final basis
    CoinWarmStartBasis* parent_basis = dynamic_cast<CoinWarmStartBasis*>(original_basis->clone());
    CbcNode* node = model_->tree()->nodePointer(tmp_ind);
    CbcNodeInfo* nodeInfo = node->nodeInfo();
    nodeInfo->buildRowBasis(*parent_basis);
    finalNodeIndices_.push_back(stats_[nodeInfo->nodeNumber()].id);

    // orig_node_id is node's id the first time we saw it (this differs from
    // node_id if one branching direction of the node has already been explored)
    const int node_id = stats_[nodeInfo->nodeNumber()].id;
    const int orig_node_id = stats_[node_id].orig_id;

    // we have the basis - now just need to recreate the bounds
    // Change bounds in the solver for the root node variable fixes
    const int init_num_changed_bounds = common_var.size();
    for (int i = 0; i < init_num_changed_bounds; i++) {
      addVarBound(tmpSolverBase, common_var[i], common_bound[i],
                  common_value[i]);
    }

    // order branching decisions leading to node so they're ordered in disjunctive terms
    sortBranchingDecisions(node_id);

    // find and create any disjunctive terms pruned due to strong branching between
    // the root and the nodes we're about to create disjunctive terms for
    recursivelyCreateStrongBranchingTerms(node_id, common_var, common_bound,
                                          common_value, tmpSolverBase);

    // Change bounds in the solver for variable bounds tightened between the root
    // node and this node's children (i.e. the nodes for the terms we're about to create)
    const int curr_num_changed_bounds =
        (orig_node_id == 0) ? 0 : stats_[orig_node_id].changed_var.size();
    std::vector<int> term_bound = common_bound;
    std::vector<double> term_value = common_value;
    std::vector<int> term_var = common_var;
    for (int i = 0; i < curr_num_changed_bounds; i++) {
      addVarBound(tmpSolverBase, stats_[orig_node_id].changed_var[i],
                  stats_[orig_node_id].changed_bound[i],
                  stats_[orig_node_id].changed_value[i], term_var,
                  term_bound, term_value);
    }

    // set the warm start - tmpSolverBase should now have bounds such that parent_basis is optimal for it
    setWarmStart(tmpSolverBase, parent_basis);

    // Now we check the branch or branches left for this node
    bool hasFeasibleChild = false;
    const int branching_index = stats_[node_id].branch_index;
    const int branching_variable = stats_[node_id].variable;
    int branching_way = stats_[node_id].way;
    double branching_value = (branching_way <= 0) ? stats_[node_id].ub : stats_[node_id].lb;

#ifdef TRACE
    printf("\n## Solving first child of parent %d/%d (term %d). ##\n",
        tmp_ind + 1, this->numNodesOnTree_, this->owner->num_terms);
#endif
    hasFeasibleChild = setupDisjunctiveTerm(
        tmpSolverBase, term_var, term_bound, term_value, branching_variable,
        branching_way == 1 ? 0 : 1, branching_value, "leaf", node_id);

    if (branching_index == 0) { // should we compute a second branch?
#ifdef TRACE
      printf("\n## Solving second child of parent %d/%d (term %d). ##\n",
          tmp_ind + 1, this->numNodesOnTree_, this->owner->num_terms);
#endif
      branching_way = (stats_[node_id].way <= 0) ? 1 : 0;
      branching_value = (branching_way <= 0) ? branching_value - 1 : branching_value + 1;
      const bool childIsFeasible = setupDisjunctiveTerm(
          tmpSolverBase, term_var, term_bound, term_value, branching_variable,
          branching_way == 1 ? 0 : 1, branching_value, "leaf", node_id);
      hasFeasibleChild = hasFeasibleChild || childIsFeasible;
    } // end second branch computation

    if (hasFeasibleChild) {
      if (stats_[node_id].depth + 1 < owner->data.min_node_depth) {
        owner->data.min_node_depth = stats_[node_id].depth + 1;
      }
      if (stats_[node_id].depth + 1 > owner->data.max_node_depth) {
        owner->data.max_node_depth = stats_[node_id].depth + 1;
      }
    }

    if (tmpSolverBase)
      delete tmpSolverBase;
    if (parent_basis)
      delete parent_basis;
  } // loop over num nodes on tree // DONE

  if (original_basis) {
    delete original_basis;
  }

  // validate the disjunction represents a full binary tree if we don't have a solution
  // and didn't error out. We'll discard the disjunction later in case of solution or error
  // (discarding for solution is just because we haven't put the thought into how to handle yet)
  if (this->owner->integer_sol.size() == 0){
#ifdef TRACE
    // print number of terms with type "pruned", "complement", and "leaf"
    int p = 0;
    int l = 0;
    int c = 0;
    for (const DisjunctiveTerm& term : owner->terms){
      if (term.type == "pruned"){
        p++;
      } else if (term.type == "leaf"){
        l++;
      } else {
        verify(term.type == "complement", "term type must be pruned, complement, or leaf");
        c++;
      }
    }
    printf("Number of pruned terms: %d\n", p);
    printf("Number of complement terms: %d\n", c);
    printf("Number of leaf terms: %d\n", l);
    printf("Number of expected pruned and complement terms: %d\n", numPrunedNodes_);
    printf("Number of expected leaf terms: %d\n", numLeafNodes_);
#endif
    isFullBinaryTree();
    verify(numLeafNodes_ + numPrunedNodes_ == owner->num_terms,
           "the size of our disjunction is not what we expected it to be");
  }
  return 0;
} /* saveInformationWithPrunes */

/**
 * @details check if the disjunction represents the leaves of a full binary tree.
 * This function assumes that the branching decisions in each disjunctive term
 * are sorted in the order they occurred.
 *
 * @return true if the tree of the disjunction is complete, false otherwise
 */
void VPCEventHandler::isFullBinaryTree(){

  std::vector<DisjunctiveTerm> terms = owner->terms;

  // check that term does not contain term2
  // we can get away with just checking for ancestral relationship because it
  // isn't possible to end up with the same branching decisions just in different order
  for (const DisjunctiveTerm& term : terms){
    for (const DisjunctiveTerm& term2 : terms){

      // term can only contain term2 if term2 is at least as deep in the tree
      if (&term != &term2 && term.changed_var.size() <= term2.changed_var.size()){

        // make sure that term is not an ancestor of term2 by checking they branch
        // on different variables or at least in different directions
        verify(!std::equal(term2.changed_var.begin(), term2.changed_var.begin() +
                           term.changed_var.size(), term.changed_var.begin()) ||
               !std::equal(term2.changed_bound.begin(), term2.changed_bound.begin() +
                           term.changed_bound.size(), term.changed_bound.begin()),
               "the LP relaxation of term contains the LP relaxation of term2");
      }
    }
  }

  // get the maximum depth
  int max_depth = 0;
  for (DisjunctiveTerm term : terms){
    max_depth = term.changed_var.size() > max_depth ? term.changed_var.size() : max_depth;
  }

  // check each leaf has a sibling
  for (int depth = max_depth; depth > 0; depth--){

    // get the terms at this depth
    std::vector<DisjunctiveTerm> depth_terms;
    for (DisjunctiveTerm term : terms){
      if (term.changed_var.size() == depth){
        depth_terms.push_back(term);
      }
    }

    // keep a running list of paired terms
    std::set<const DisjunctiveTerm*> paired_terms;

    // find a sibling for each term
    for (const DisjunctiveTerm& term : depth_terms){

      // Check if the term was found to be another's sibling earlier
      if (paired_terms.find(&term) != paired_terms.end()){
        continue;
      }

      // check term against all other terms at this depth
      for (const DisjunctiveTerm& term2 : depth_terms){
        std::vector<int> differing_idx;
        
        // we can only be siblings if we share the same variables that were
        // branched on and we're not the same term
        if (term.changed_var == term2.changed_var && &term != &term2){
          
          // record branching directions we differ on
          for (int i = 0; i < depth; i++) {
            if (term.changed_bound[i] != term2.changed_bound[i]) {
              differing_idx.push_back(i);
            }
          }

          // if we differ by only the last branching decision, we're siblings
          if (differing_idx.size() == 1 && differing_idx[0] == depth - 1){
            // check that we have the reciprocal branching decision
            double expected_val = term.changed_value[differing_idx[0]] +
                                  (term.changed_bound[differing_idx[0]] == 0 ? -1 : 1);
            verify(term2.changed_value[differing_idx[0]] == expected_val,
                   "term2 branch value doesnt meet expectation");

            // record siblings
            paired_terms.insert(&term);
            paired_terms.insert(&term2);

            // create a parent node to leave in the next level above
            DisjunctiveTerm parent_term = term;
            parent_term.changed_var.erase(parent_term.changed_var.begin() + differing_idx[0]);
            parent_term.changed_bound.erase(parent_term.changed_bound.begin() + differing_idx[0]);
            parent_term.changed_value.erase(parent_term.changed_value.begin() + differing_idx[0]);
            parent_term.type = "parent";
            terms.push_back(parent_term);
            break;
          }
        }
      }

      // if we didn't find a sibling, the disjunction is not complete
      if (paired_terms.find(&term) == paired_terms.end()){
        verify(false, "Disjunction does not represent a full binary tree.");
      }
      
    } // find a sibling for each term
  }
}

