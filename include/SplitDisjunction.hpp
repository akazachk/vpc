// Name:     SplitDisjunction.hpp
// Author:   A. M. Kazachkov
// Date:     2018-02-24
//-----------------------------------------------------------------------------
#pragma once

/****************************************/
/*  Disjunction generated from a split  */
/****************************************/

#include "Disjunction.hpp"
#include "TimeStats.hpp"
#include "VPCParameters.hpp"

class SplitDisjunction : public Disjunction {
public:
  int var;

  /** Params */
  VPCParameters params;
  /** Param constructor */
  SplitDisjunction(const VPCParameters& params);
  /** setParams based on VPCParameters */
  void setParams(const VPCParameters& params);

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
  virtual ExitReason prepareDisjunction(const OsiSolverInterface* const si);
protected:
  void initialize(const SplitDisjunction* const source = NULL, const VPCParameters* const params = NULL);
  bool checkVar(OsiSolverInterface* si, int col);
  void setCgsName(const int var, const double val);
  void addTerm(const int branching_variable,
      const int branching_way, const double branching_value,
      const OsiSolverInterface* const solver);
}; /* SplitDisjunction */
