/**
 * @file SplitDisjunction.hpp
 * @author A. M. Kazachkov
 * @date 2018-02-24
 */
#pragma once

#include "VPCDisjunction.hpp"

/// @brief Disjunction generated from a (variable) split
class SplitDisjunction : public VPCDisjunction {
public:
  int var; ///< variable on which we are taking split

  /// @brief Param constructor
  SplitDisjunction(const VPCParametersNamespace::VPCParameters& params);

  /// @brief Copy and param constructor
  SplitDisjunction(const SplitDisjunction& source, const VPCParametersNamespace::VPCParameters& params);

  /// @brief Default constructor
  SplitDisjunction();

  /// @brief Copy constructor
  SplitDisjunction(const SplitDisjunction& source);

  /// @brief Destructor
  ~SplitDisjunction();

  /// @brief Assignment operator
  SplitDisjunction& operator=(const SplitDisjunction& source);

  /// @brief Clone
  virtual SplitDisjunction* clone() const;

  /// @brief For clearing things and setting up the disjunction as new
  virtual void setupAsNew();

  /// @brief Prepare a new disjunction
  virtual DisjExitReason prepareDisjunction(const OsiSolverInterface* const si);

protected:
  /// @brief Setup class members (copy from \p source if provided)
  void initialize(const SplitDisjunction* const source = NULL, const VPCParametersNamespace::VPCParameters* const params = NULL);
  /// @brief Check if var is fractional
  bool checkVar(OsiSolverInterface* si, int col);
  /// @brief Set #name
  void setCgsName(const int var, const double val);
  /// @brief Add disjunctive term
  void addTerm(const int branching_variable,
      const int branching_way, const double branching_value,
      const OsiSolverInterface* const solver);
}; /* SplitDisjunction */

/// @brief Generate SplitDisjunction for all variables that are fractional, using hot starts
int generateSplitDisjunctions(
    DisjunctionSet* const disjSet,
    const OsiSolverInterface* const si,
    const VPCParametersNamespace::VPCParameters& params);