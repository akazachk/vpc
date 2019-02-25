// Name:     PartialBBDisjunction.hpp
// Author:   A. M. Kazachkov
// Date:     2018-02-22
//-----------------------------------------------------------------------------
#pragma once

/****************************************************************/
/*  Disjunction generated from a partial branch-and-bound tree  */ 
/****************************************************************/

#include "Disjunction.hpp"
#include "VPCParameters.hpp"

#include <limits> // numeric_limits

class PartialBBDisjunction : public Disjunction {
public:
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
