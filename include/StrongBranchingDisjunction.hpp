/**
 * @file StrongBranchingDisjunction.hpp
 * @author Sean Kelley
 * @date 2023-04-20
 */
#pragma once

// standard library
#include "vector"

// VPC modules
#include "VPCDisjunction.hpp"

/// @brief Disjunction generated from a (variable) split
class StrongBranchingDisjunction : public VPCDisjunction {
public:

  /// @brief Param constructor
  StrongBranchingDisjunction(const VPCParametersNamespace::VPCParameters& params);

  /// @brief Copy and param constructor
  StrongBranchingDisjunction(const StrongBranchingDisjunction& source, const VPCParametersNamespace::VPCParameters& params);

  /// @brief Default constructor
  StrongBranchingDisjunction();

  /// @brief Copy constructor
  StrongBranchingDisjunction(const StrongBranchingDisjunction& source);

  /// @brief Destructor
  ~StrongBranchingDisjunction();

  /// @brief Assignment operator
  StrongBranchingDisjunction& operator=(const StrongBranchingDisjunction& source);

  /// @brief Clone
  virtual StrongBranchingDisjunction* clone() const;

  /// @brief For clearing things and setting up the disjunction as new
  virtual void setupAsNew();

  /// @brief Get disjunction
  virtual DisjExitReason prepareDisjunction(const OsiSolverInterface* const si);

protected:
  /// @brief Setup class members (copy from \p source if provided)
  void initialize(const StrongBranchingDisjunction* const source = NULL,
                  const VPCParametersNamespace::VPCParameters* const params = NULL);
  /// @brief Set #name
  void setCgsName(const int var, const double val);
  /// @brief Add disjunctive term
  void addTerm(const std::vector<int>& branching_variables,
      const std::vector<int>& branching_ways, const std::vector<double>& branching_values,
      OsiSolverInterface* solver);
  /// @brief Get the cartesian product of a vector of vectors
  template <typename T>
  void cartesianProduct(
      /// [out] empty vector of vectors that will return as cartesian product
      std::vector< std::vector<T> >& finalResult,
      /// [out] empty vector that will have work done in it
      std::vector<T>& currentResult,
      /// [in] begin iterator for input vector of vectors to cartesian product
      typename std::vector< std::vector<T> >::const_iterator currentInput,
      /// [in] end iterator for input vector of vectors to cartesian product
      typename std::vector< std::vector<T> >::const_iterator lastInput);
}; /* StrongBranchingDisjunction */
