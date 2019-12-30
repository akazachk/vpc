/**
 * @file VPCDisjunction.hpp
 * @author A. M. Kazachkov
 * @date 2019-12-26
 */
#pragma once

#include <string>
#include <vector>

// Project files
#include "Disjunction.hpp"
#include "TimeStats.hpp"
#include "VPCParameters.hpp"

/**********************************************/
/*  Generic abstract VPCDisjunction class     */
/*  from which a point-ray collection can be  */
/*  constructed and VPCs can be generated     */
/**********************************************/
class VPCDisjunction : public Disjunction {
public:
  friend class CglVPC; // to access protected members

  VPCParametersNamespace::VPCParameters params; // parameters for the class

  // Optional members
  TimeStats* timer; // not owned by VPCDisjunction

  /** setParams based on VPCParameters */
  void setParams(const VPCParametersNamespace::VPCParameters& params);

  /** Param constructor */
  VPCDisjunction(const VPCParametersNamespace::VPCParameters& params);

  /** Copy and param constructor */
  VPCDisjunction(const VPCDisjunction& source, const VPCParametersNamespace::VPCParameters& params);

  /** Default constructor */
  VPCDisjunction();

  /** Copy constructor */
  VPCDisjunction(const VPCDisjunction& source);

  /** Destructor */
  virtual ~VPCDisjunction();

  /** Assignment operator */
  VPCDisjunction& operator=(const VPCDisjunction& source);

  /** Clone */
  virtual VPCDisjunction* clone() const = 0;

  /** For clearing things and setting up the disjunction as new */
  virtual void setupAsNew();

  /** Get disjunction */
#ifdef USE_COIN
  virtual DisjExitReason prepareDisjunction(const OsiSolverInterface* const si) = 0;
#else
  virtual DisjExitReason prepareDisjunction() = 0;
#endif

protected:
  void initialize(const VPCDisjunction* const source = NULL, const VPCParametersNamespace::VPCParameters* const params = NULL);
}; /* VPCDisjunction */
