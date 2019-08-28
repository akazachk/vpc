// Name:     PartialBBDisjunction.hpp
// Author:   A. M. Kazachkov
// Date:     2018-02-22
//-----------------------------------------------------------------------------
#pragma once

#include <limits> // numeric_limits

#include "Disjunction.hpp"
#include "VPCParameters.hpp"

class OsiSolverInterface;
class PartialBBDisjunction;

#ifdef USE_CBC
class CbcModel;
class CbcEventHandler;

/**
 * Set parameters for Cbc used for VPCs, as well as the custom branching decision
 */
void setCbcParametersForPartialBB(const VPCParameters& param,
    CbcModel* const cbc_model, CbcEventHandler* eventHandler = NULL,
    const int numStrong = 5, const int numBeforeTrusted = 10,
    const double max_time = std::numeric_limits<double>::max());

/**
 * Generate a partial branch-and-bound tree with at most max_leaf_nodes leaf nodes
 */
void generatePartialBBTree(PartialBBDisjunction* const owner,
    CbcModel* cbc_model, const OsiSolverInterface* const solver,
    const int max_nodes, const int num_strong, const int num_before_trusted);
#endif // USE_CBC

/**
 * Keeps split information at root, and best bound on each side
 */
struct RootTerm {
  int var;       ///< var branched on
  double val;    ///< value of var when it is branched on
  double boundD; ///< best bound in down branch
  double boundU; ///< best bound in up branch

  /** Constructor (calls initialize) */
  RootTerm() {
    initialize();
  } /* constructor */

  /** Reset root term information and deletes any allocated memory */
  void clear() {
    var = -1;
    val = 0.;
    boundD = std::numeric_limits<double>::max();
    boundU = std::numeric_limits<double>::max();
  } /* clear */

  /** Initalize root term (calls clear) */
  void initialize() {
    clear();
  } /* initialize */
}; /* RootTerm */

/****************************************************************/
/*  Disjunction generated from a partial branch-and-bound tree  */
/*  Currently requires the use of Cbc                           */
/****************************************************************/
class PartialBBDisjunction : public Disjunction {
public:
  RootTerm root; /// Keeps split information at root and best bound on each side

  struct PartialBBDisjunctionData {
    int num_nodes_on_tree = 0, num_partial_bb_nodes = 0, num_pruned_nodes = 0;
    int min_node_depth = std::numeric_limits<int>::max();
    int max_node_depth = 0;
    int num_fixed_vars = 0;
  } data;

  /** Params */
  VPCParameters params;
  /** Param constructor */
  PartialBBDisjunction(const VPCParameters& params);
  /** setParams based on VPCParameters */
  void setParams(const VPCParameters& params);

  /** Default constructor */
  PartialBBDisjunction();

  /** Copy constructor */
  PartialBBDisjunction(const PartialBBDisjunction& source);

  /** Destructor */
  ~PartialBBDisjunction();

  /** Assignment operator */
  PartialBBDisjunction& operator=(const PartialBBDisjunction& source);

  /** Clone */
  virtual PartialBBDisjunction* clone() const;

  /** For clearing things and setting up the disjunction as new */
  virtual void setupAsNew();

  /** Get disjunction */
  virtual ExitReason prepareDisjunction(const OsiSolverInterface* const si);
protected:
  void initialize(const PartialBBDisjunction* const source = NULL, const VPCParameters* const params = NULL);
}; /* PartialBBDisjunction */
