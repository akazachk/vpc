/**
 * @file CompleteDisjunction.hpp
 * @author Sean Kelley
 * @date 2023-05-03
 */
#pragma once

// standard library
#include "vector"

// VPC modules
#include "SolverInterface.hpp" // SolverInterface
#include "VPCDisjunction.hpp"

/// @brief A Disjunction with all variable fixes in a partial disjunction included as additional terms
class CompleteDisjunction : public VPCDisjunction {
public:

  /// @brief Disjunction to be completed
  Disjunction * partialDisj;

  /// @brief Default constructor
  CompleteDisjunction();

  /// @brief Copy constructor
  CompleteDisjunction(const CompleteDisjunction& source);

  /// @brief Destructor
  ~CompleteDisjunction();

  /// @brief Assignment operator
  CompleteDisjunction& operator=(const CompleteDisjunction& source);

  /// @brief Clone
  CompleteDisjunction* clone() const;

  /// @brief For clearing things and setting up the disjunction as new
  void setupAsNew();

  /// @brief Get disjunction
  DisjExitReason prepareDisjunction(const OsiSolverInterface* const si,
                                    const VPCDisjunction* const partialDisj);

  /// @brief add the terms arising from the tree that created the fixed variables
  /// common for all disjunctive terms
  void addCommonFixedTerms(OsiSolverInterface* solver);

  /// @brief add the terms arising from the tree that created fixed variables for each term
  void addDisjunctionFixedTerms(OsiSolverInterface* solver);

  /// @brief finds and removes redundant terms from the disjunction
  void removeRedundantTerms(const OsiSolverInterface* const si);

protected:
  /// @brief Setup class members (copy from \p source if provided)
  void initialize(const CompleteDisjunction* const source = NULL);
  /// @brief Get disjunction
  virtual DisjExitReason prepareDisjunction(const OsiSolverInterface* const si);
  /// @brief Complete the tree encoded in the disjunction
  void getTerms(
      /// [out] the variables to fix for each term
      std::vector<std::vector<int> >& variablesByTerm,
      /// [out] the bounds to fix for each term
      std::vector<std::vector<int> >& boundsByTerm,
      /// [out] the values to fix respective bounds for each term
      std::vector<std::vector<double> >& valuesByTerm,
      /// [in] which variable had a bound changed
      const std::vector<int>& changedVariables,
      /// [in] which bound was fixed (0 - lower; 1 - upper)
      const std::vector<int>& changedBounds,
      /// [in] which value was the bound fixed to
      const std::vector<double>& changedValues,
      /// [in] which variables are fixed for these terms
      const std::vector<int>* fixedVariables=NULL,
      /// [in] which bounds are fixed for these terms
      const std::vector<int>* fixedBounds=NULL,
      /// [in] which values are the fixed bounds for these terms set to
      const std::vector<double>* fixedValues=NULL);
  /// @brief Add disjunctive term
  void addTerm(const std::vector<int>& branching_variables,
               const std::vector<int>& branching_ways, const std::vector<double>& branching_values,
               OsiSolverInterface* solver);

}; /* CompleteDisjunction */