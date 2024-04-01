/**
 * @file VPCEventHandler.hpp
 * @brief Custom event handler for Cbc for partial b&b tree generation (requires use of USE_CBC)
 *
 * @author A. M. Kazachkov
 * @date 2018-12-25
 */
#pragma once

#include <set>
#include <string>
#include <vector>
#include <memory> // unique_ptr

class PartialBBDisjunction;
class OsiCuts;
#include "Disjunction.hpp"
#include "SolverInterface.hpp"

#ifdef USE_CBC
#include <CbcEventHandler.hpp>
#include <CbcNode.hpp>

class VPCEventHandler;
enum class PruneNodeOption;

/// @brief Useful data structure for COIN-OR tracking
struct NodeStatistics {
  int id = -1; ///< in stats vector
  int parent_id = -1; ///< in stats vector
  int variable = -1; ///< variable we will branch on if this node is popped
  int branch_index = -1; ///< how many times we have branched from this node
  int way = -2; ///< which way we will branch if this node is popped (+1 = up, -1 = down)
  double obj = std::numeric_limits<double>::max(); ///< objective at this node
  double value = -1; ///< value of the variable at the solution of the parent, I believe
  double lb = -1; ///< new lb
  double ub = -1; ///< new ub
  int orig_id = -1; ///< first time this node appeared in stats
  int number = -1; ///< this is the number in the tree, not stats vector (differs from Cbc count)
  int depth = -1; ///< depth of the node in the tree
  bool found_integer_solution = false; ///< was an integer-feasible solution found during this node exploration? (may have been during strong branching)
  double integer_obj = std::numeric_limits<double>::max(); ///< integer objective at this node
  std::string changed_bounds; ///< string concatenating bounds changed at this node
  std::vector<int> changed_var; ///< vector of variable indices changed to get to this node
  std::vector<int> changed_bound; ///< vector with which bound (0: lower, 1: upper) changed for the changed variables
  std::vector<double> changed_value; ///< vector of new lower/upper bound for each variable that has been changed
  std::string status; ///< status of the node (normal or pruned by integrality / bound / infeasibility)
}; /* Node Statistics */

/// @brief Define the statistics for a particular node
void setNodeStatistics(NodeStatistics& stats, const CbcNode* const node,
    const CbcModel* const model, const std::vector<NodeStatistics>& stats_vec,
    const std::vector<double>& originalLB,
    const std::vector<double>& originalUB, const bool will_create_child,
    const bool solution_is_found = false);
/// @brief Define the statistics to be saved for a pruned node
void setPrunedNodeStatistics(NodeStatistics& stats,
    const CbcNode* const parent_node, const CbcModel* const model,
    /*const std::vector<NodeStatistics>& stats_vec,*/const bool prune_by_integrality,
    const std::vector<double>& originalLB,
    const std::vector<double>& originalUB, const double obj,
    const bool solution_is_found,
    const PruneNodeOption& pruneReason);

/// @brief Call printNodeStatistics(const std::vector<NodeStatistics>&, const bool, FILE*)
void printNodeStatistics(const std::vector<NodeStatistics>& stats);
/// @brief Print node statistics for a vector of NodeStatistics
void printNodeStatistics(const std::vector<NodeStatistics>& stats,
    const bool print_bounds, FILE* myfile = stdout);
/// @brief Print a given node's statistics
void printNodeStatistics(const NodeStatistics& stats, const bool print_bounds,
    FILE* myfile = stdout);

/// @brief Save which bounds have been changed with respect to the original solver
void changedBounds(NodeStatistics& stats,
    const OsiSolverInterface* const solver,
    const std::vector<double>& originalLB,
    const std::vector<double>& originalUB);

/// @brief This is so user can trap events and do useful stuff.
/// CbcModel model_ is available as well as anything else you care to pass in.
class VPCEventHandler: public CbcEventHandler {
public:
  PartialBBDisjunction* owner;
  /// @brief Reason a node is pruned
  enum class PruneNodeOption {
    NO = 0, PRUNE_BY_INFEASIBILITY, PRUNE_BY_BOUND, PRUNE_BY_INTEGRALITY
  }; /* PruneNodeOption */

  /**@name Overrides */
  //@{
  /// @brief Custom event handler that keeps history of the tree
  virtual CbcAction event(CbcEvent whichEvent);
  //@}

  /**@name Constructors, destructor, etc. */
  //@{
  /** @brief Default constructor */
  VPCEventHandler();
  /** @brief VPC special constructors */
  VPCEventHandler(PartialBBDisjunction* const disj, const int maxNumLeafNodes,
      const double maxTime);
  /// @brief Constructor with pointer to model (redundant as setEventHandler does)
  VPCEventHandler(CbcModel* model);
  /** Destructor */
  virtual ~VPCEventHandler();
  /** @brief The copy constructor. */
  VPCEventHandler(const VPCEventHandler& rhs);
  /// @brief Assignment
  VPCEventHandler& operator=(const VPCEventHandler& rhs);
  /// @brief Clone
  virtual CbcEventHandler* clone() const;
  //@}

  /**@name Access methods */
  //@{
  // @brief Number of nodes on the tree at the end
  inline int getNumNodesOnTree() const {
    return numNodesOnTree_;
  }
  inline int getNumNodes() const {
    return numNodes_;
  }
  /// @brief For the nodes on the tree at the end, where do we find the information in #stats_
  inline int getNodeIndex(const int i) const {
    return finalNodeIndices_[i];
  }
  inline int getNumLeafNodes() const {
    return numLeafNodes_;
  }
  inline int getMaxNumLeafNodes() const {
    return maxNumLeafNodes_;
  }
  inline int getMaxTime() const {
    return maxTime_;
  }
  inline SolverInterface* getOriginalSolver() const {
    return originalSolver_;
  }
  /// @brief set #originalSolver_ to a clone of copySolver, and call #setOriginalLB and #setOriginalUB
  SolverInterface* setOriginalSolver(const OsiSolverInterface* const copySolver,
      const bool return_old = false);
  inline int getOriginalLB(const int col) const {
    return originalLB_[col];
  }
  inline int getOriginalUB(const int col) const {
    return originalUB_[col];
  }
  void setOriginalLB(const int num_cols, const double* const vec);
  void setOriginalUB(const int num_cols, const double* const vec);
//  inline const CoinWarmStartBasis* getBasisForNode(const int node_ind) const {
//    return bases_[node_ind];
//  }
  inline std::vector<NodeStatistics> getStatsVector() const {
    return stats_;
  }
  inline std::vector<NodeStatistics> getPrunedStatsVector() const {
    return pruned_stats_;
  }
  inline const double* getIntegerFeasibleSolution() const {
    return savedSolution_.data();
  }
  inline int getNumStats() const {
    return stats_.size();
  }
  inline int getNumPrunedStats() const {
    return pruned_stats_.size();
  }
  inline bool isIntegerSolutionFound() const {
    return !savedSolution_.empty();
  }
  //@}

protected:
  int maxNumLeafNodes_; ///< stop when we reach this many leaf nodes
  int numCuts_; ///< number of cuts generated
  double maxTime_; ///< maximum allowable time

  ///@{
  /// @name Things that will be saved at the end
  int numNodesOnTree_; ///< number of active nodes on the tree
  int numLeafNodes_; ///< number of leaf nodes
  int numNodes_; ///< total number of nodes explored
  double obj_; ///< integer-feasible solution value
  SolverInterface* originalSolver_; ///< pointer to original solver to compare bounds
  std::vector<double> originalLB_; ///< vector of original lower bounds at the *end* of the root node processing
  std::vector<double> originalUB_; ///< vector of original upper bounds at the *end* of the root node processing
//    std::vector<CoinWarmStartBasis*> bases_; ///< bases of nodes
  std::vector<NodeStatistics> stats_; ///< all stats that we need to recreate tree
  std::vector<NodeStatistics> pruned_stats_; ///< info for all children that were pruned
  std::vector<int> finalNodeIndices_; ///< node numbers for the nodes on the final tree
  std::vector<double> savedSolution_; ///< when pruneNode_ = 3, the saved solution might have been deleted somehow
  std::unique_ptr<OsiCuts> cuts_; ///< cuts we generate
  ///@}

  ///@{
  /// @name Temporary information we want to keep during the B&B process
  std::vector<CbcNode*> currentNodes_; ///< current set of nodes on the tree
  CbcNode* parent_; ///< pointer to parent of current node being explored
  CbcNodeInfo* parentInfo_; ///< full node info for parent
  CbcNode* child_; ///< current child being explored
  PruneNodeOption pruneNode_ = PruneNodeOption::NO; ///< reason a node is pruned
  bool reachedEnd_ = false; ///< are we done?
  bool foundSolution_ = false; ///< was an integer-feasible solution found?
  std::set<int> checked_nodes_; ///< nodes we have already checked for strong branching fixes
  std::set<int> sorted_nodes_; ///< nodes we have already sorted the branching decisions for
  int numPrunedNodes_; ///< number of nodes implicitly created and pruned during strong branching
  ///@}

  ///@{
  /// @name Helper methods

  /// @brief Copy our stuff
  void initialize(const VPCEventHandler* const rhs);

  /// @brief Update the dual bound of the disjunction if the term worsens it
  void updateDualBound(const int orig_node_id, const DisjunctiveTerm& term);

  /// @brief Solve for a disjunctive term
  bool setupDisjunctiveTerm(const int node_id, const int branching_variable,
      const int branching_way, const double branching_value,
      const SolverInterface* const tmpSolverBase,
      const int curr_num_changed_bounds,
      std::vector<std::vector<int> >& commonTermIndices,
      std::vector<std::vector<double> >& commonTermCoeff,
      std::vector<double>& commonTermRHS);

  /// @brief Set the warm start information for the parent node
  void setWarmStart(SolverInterface* tmpSolverBase, CoinWarmStartBasis* parent_basis,
                    int tmp_ind = -1, int node_id = -1, const int curr_num_changed_bounds = 0,
                    const std::vector<std::vector<int> >& commonTermIndices = std::vector<std::vector<int> >(),
                    const std::vector<std::vector<double> >& commonTermCoeff = std::vector<std::vector<double> >(),
                    const std::vector<double>& commonTermRHS = std::vector<double>());

  /// @brief Clear saved information and free memory
  void clearInformation();

  /// @brief Save node information before we reach endSearch_
  int saveInformation();

  /// @brief drop all tightened bounds occurring on continuous variables.
  void removeContinuousVariableTightenedBounds(std::vector<NodeStatistics>& stat_vec);

  /// @brief reorder the branching decisions leading to each node in the tree in the order that they occurred
  void sortBranchingDecisions(const int node_id);

  /// @brief Creates a disjunctive term in a generalized enough manner to allow for pruned terms
  bool setupDisjunctiveTerm(
      const SolverInterface* const tmpSolverParent, const std::vector<int>& term_var,
      const std::vector<int>& term_bound, const std::vector<double>& term_val,
      const int branching_variable, const int branching_bound, const double branching_value,
      const std::string& term_type, const int orig_node_id = -1);

  /// @brief Create disjunctive terms due to strong branching on the child node
  void createStrongBranchingTerms(
      std::vector<int> child_pre_branch_var, std::vector<int> child_pre_branch_bound,
      std::vector<double> child_pre_branch_val, const SolverInterface* const tmpSolver,
      std::vector<int> parent_var = std::vector<int>(),
      std::vector<int> parent_bound = std::vector<int>(),
      std::vector<double> parent_val = std::vector<double>());

  /// @brief Create disjunctive terms pruned from strong branching between node
  // <node_id> and the root node
  void recursivelyCreateStrongBranchingTerms(
    const int node_id, const std::vector<int>& common_var,
    const std::vector<int>& common_bound, const std::vector<double>& commo_value,
    const SolverInterface* const tmpSolverBase);

  /// @brief Save node information (including those pruned) before we reach endSearch_
  int saveInformationWithPrunes();

  /// @brief check if the disjunction represents the leaves of a full binary tree
  void isFullBinaryTree();
  //@}
};
/* VPCEventHandler definition */
#endif // USE_CBC
