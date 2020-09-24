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

/**
 * @brief Generic abstract VPCDisjunction class
 * 
 * Allows a point-ray collection to be constructed,
 * from which we can generate VPCs
 */
class VPCDisjunction : public Disjunction {
public:
  friend class CglVPC; ///< to access protected members

  VPCParametersNamespace::VPCParameters params; ///< parameters for the class

  /// @name Optional members
  ///@{
  TimeStats* timer; ///< not owned by VPCDisjunction
  ///@}

  /// @brief setParams based on VPCParameters
  void setParams(const VPCParametersNamespace::VPCParameters& params);

  /// @brief Param constructor
  VPCDisjunction(const VPCParametersNamespace::VPCParameters& params);

  /// @brief Copy and param constructor
  VPCDisjunction(const VPCDisjunction& source, const VPCParametersNamespace::VPCParameters& params);

  /// @brief Default constructor
  VPCDisjunction();

  /// @brief Copy constructor
  VPCDisjunction(const VPCDisjunction& source);

  /// @brief Destructor
  virtual ~VPCDisjunction();

  /// @brief Assignment operator
  VPCDisjunction& operator=(const VPCDisjunction& source);

  /// @brief Clone
  virtual VPCDisjunction* clone() const = 0;

  /// @brief For clearing things and setting up the disjunction as new
  virtual void setupAsNew();

  /// @brief Get disjunction
#ifdef USE_COIN
  virtual DisjExitReason prepareDisjunction(const OsiSolverInterface* const si) = 0;
#else
  virtual DisjExitReason prepareDisjunction() = 0;
#endif

protected:
  /// @brief Initialize class members (copy from \p source if provided)
  void initialize(const VPCDisjunction* const source = NULL, const VPCParametersNamespace::VPCParameters* const params = NULL);
}; /* VPCDisjunction */
