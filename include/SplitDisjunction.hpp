/**
 * @file SplitDisjunction.hpp
 * @author A. M. Kazachkov
 * @date 2018-02-24
 */
#pragma once

/****************************************/
/*  Disjunction generated from a split  */
/****************************************/

#include "VPCDisjunction.hpp"

class SplitDisjunction : public VPCDisjunction {
public:
  int var;

  /** Param constructor */
  SplitDisjunction(const VPCParametersNamespace::VPCParameters& params);

  /** Copy and param constructor */
  SplitDisjunction(const SplitDisjunction& source, const VPCParametersNamespace::VPCParameters& params);

  /** Default constructor */
  SplitDisjunction();

  /** Copy constructor */
  SplitDisjunction(const SplitDisjunction& source);

  /** Destructor */
  ~SplitDisjunction();

  /** Assignment operator */
  SplitDisjunction& operator=(const SplitDisjunction& source);

  /** Clone */
  virtual SplitDisjunction* clone() const;

  /** For clearing things and setting up the disjunction as new */
  virtual void setupAsNew();

  /** Get disjunction */
  virtual DisjExitReason prepareDisjunction(const OsiSolverInterface* const si);

protected:
  void initialize(const SplitDisjunction* const source = NULL, const VPCParametersNamespace::VPCParameters* const params = NULL);
  bool checkVar(OsiSolverInterface* si, int col);
  void setCgsName(const int var, const double val);
  void addTerm(const int branching_variable,
      const int branching_way, const double branching_value,
      const OsiSolverInterface* const solver);
}; /* SplitDisjunction */
