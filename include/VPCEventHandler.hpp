//============================================================================
// Name        : VPCEventHandler.hpp
// Author      : A. M. Kazachkov
// Version     : 2018-Dec-25
// Description : Custom event handler for Cbc
//
// This code is licensed under the terms of the Eclipse Public License (EPL).
//============================================================================

#pragma once

#include <string>
#include <vector>

#include <CbcEventHandler.hpp>
#include <CbcNode.hpp>

#include "VPCParameters.hpp"

// Useful data structure for COIN-OR tracking
struct NodeStatistics {
  int id = -1; // in stats vector
  int parent_id = -1; // in stats vector
  int variable = -1; // variable we will branch on if this node is popped
  int branch_index = -1; // how many times we have branched from this node
  int way = -2; // which way we will branch if this node is popped
  double obj = std::numeric_limits<double>::max(); // objective at this node
  double value = -1; // value of the variable at the solution of the parent, I believe
  double lb = -1; // new lb
  double ub = -1; // new ub
  int orig_id = -1; // first time this node appeared in stats
  int number = -1; // this is the number in the tree, not stats vector (differs from Cbc count)
  int depth = -1; // depth of the node in the tree
  bool found_integer_solution = false; // was an integer-feasible solution found during this node exploration? (may have been during strong branching)
  double integer_obj = std::numeric_limits<double>::max();
  std::string changed_bounds;
  std::vector<int> changed_var;
  std::vector<int> changed_bound;
  std::vector<double> changed_value;
};

void setNodeStatistics(NodeStatistics& stats, const CbcNode* const node, 
    const CbcModel* const model, const std::vector<NodeStatistics>& stats_vec,
    const std::vector<double>& originalLB, const std::vector<double>& originalUB,
    const bool will_create_child, const bool solution_is_found = false);
void setPrunedNodeStatistics(NodeStatistics& stats,
    const CbcNode* const parent_node, const CbcModel* const model,
    /*const std::vector<NodeStatistics>& stats_vec,*/const bool prune_by_integrality,
    const std::vector<double>& originalLB,
    const std::vector<double>& originalUB, const double obj,
    const bool solution_is_found);
void printNodeStatistics(const std::vector<NodeStatistics>& stats);
void printNodeStatistics(const std::vector<NodeStatistics>& stats,
    const bool print_bounds, FILE* myfile = stdout);
void printNodeStatistics(const NodeStatistics& stats, const bool print_bounds,
    FILE* myfile = stdout);

//std::string changedBoundsString(const CbcNodeInfo* const nodeInfo,
//    const std::vector<double>& originalLB, const std::vector<double>& originalUB);
void changedBounds(NodeStatistics& stats, const OsiSolverInterface* const solver,
    const std::vector<double>& originalLB, const std::vector<double>& originalUB);

/************************************************************/
/** This is so user can trap events and do useful stuff.  
 *
 *     CbcModel model_ is available as well as anything else you care 
 *         to pass in
 */
class VPCEventHandler : public CbcEventHandler {
  public:
    /**@name Overrides */
    //@{
    virtual CbcAction event(CbcEvent whichEvent);
    //@}
    
    /**@name Constructors, destructor, etc. */
    //@{
    /** Default constructor */
    VPCEventHandler();
    /** VPC special constructors */
    VPCEventHandler(const int maxNumLeafNodes, const double maxTime, const VPCParameters* params);
    /// Constructor with pointer to model (redundant as setEventHandler does)
    VPCEventHandler(CbcModel* model);
    /** Destructor */
    virtual ~VPCEventHandler();
    /** The copy constructor. */
    VPCEventHandler(const VPCEventHandler& rhs);
    /// Assignment
    VPCEventHandler& operator=(const VPCEventHandler& rhs);
    /// Clone
    virtual CbcEventHandler* clone() const;
    /// Copy our stuff
    virtual void copyOurStuff(const VPCEventHandler* const rhs);
    //@}

    /**@name Helper methods */
    //@{
    virtual void saveInformation();
    //@}

    /**@name Access methods */
    //@{
    /*
    inline double getObjLower(const int i) { return lb[i]; }
    inline double getObjUpper(const int i) { return ub[i]; }
    */
    // Number of nodes on the tree at the end
    inline int getNumNodesOnTree() const { return numNodesOnTree_; }
    inline int getNumNodes() const { return numNodes_; }
    // For the nodes on the tree at the end, where do we find the information in stats_
    inline int getNodeIndex(const int i) const { return finalNodeIndices_[i]; }
    inline int getNumLeafNodes() const { return numLeafNodes_; }
    inline int getMaxNumLeafNodes() const { return maxNumLeafNodes_; }
    inline int getMaxTime() const { return maxTime_; }
//    inline CoinWarmStartBasis* getOriginalBasis() const { return originalBasis_; };
//    CoinWarmStartBasis* setOriginalBasis(const CoinWarmStart* const copyBasis, const bool return_old = false);
    inline SolverInterface* getOriginalSolver() const { return originalSolver_; }
    SolverInterface* setOriginalSolver(
      const OsiSolverInterface* const copySolver,
      const bool return_old = false);
    inline int getOriginalLB(const int col) const { return originalLB_[col]; }
    inline int getOriginalUB(const int col) const { return originalUB_[col]; }
    void setOriginalLB(const int num_cols, const double* const vec);
    void setOriginalUB(const int num_cols, const double* const vec);
    inline const CoinWarmStartBasis* getBasisForNode(const int node_ind) const {
      return bases_[node_ind];
    }
    inline std::vector<NodeStatistics> getStatsVector() const { return stats_; }
    inline std::vector<NodeStatistics> getPrunedStatsVector() const { return pruned_stats_; }
    inline const double* getIntegerFeasibleSolution() const { return savedSolution_.data(); }
    inline int getNumStats() const { return stats_.size(); }
    inline int getNumPrunedStats() const { return pruned_stats_.size(); }
    inline bool isIntegerSolutionFound() const { return !savedSolution_.empty(); }
    //@}
    
  protected:
    const VPCParameters* params;
    int maxNumLeafNodes_;
    double maxTime_;
    
    // Things that will be saved at the end
    int numNodesOnTree_, numLeafNodes_;
    int numNodes_;
    //CoinWarmStartBasis* originalBasis_;
    SolverInterface* originalSolver_;
    std::vector<double> originalLB_, originalUB_;
    std::vector<CoinWarmStartBasis*> bases_; // bases of node
    std::vector<NodeStatistics> stats_; // all stats that we need to recreate tree
    std::vector<NodeStatistics> pruned_stats_; // info for all children that were pruned
    std::vector<int> finalNodeIndices_; // node numbers for the nodes on the final tree
    //std::vector<std::vector<double>> integerFeasibleSolutions_; // for those pruned nodes that were actually feasible
    std::vector<double> savedSolution_; // when pruneNode_ = 3, the saved solution might have been deleted somehow

    // Temporary information we want to keep during the B&B process
    std::vector<CbcNode*> currentNodes_;
    CbcNode* parent_;
    CbcNodeInfo* parentInfo_;
    CbcNode* child_;
    int pruneNode_ = 0; // 0: no, 1: prune by infeasibility, 2: prune by bound, 3: prune by integrality
    bool reachedEnd_ = false;
    bool foundSolution_ = false;
}; /* VPCEventHandler definition */
