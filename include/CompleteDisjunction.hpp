/**
 * @file CompleteDisjunction.hpp
 * @author Sean Kelley
 * @date 2023-05-03
 */
#pragma once

// standard library
#include "vector"

// VPC modules
#include "Disjunction.hpp"

/// @brief A Disjunction with all variable fixes in a partial disjunction included as additional terms
class CompleteDisjunction : public Disjunction {
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
                                    const Disjunction* const partialDisj);

protected:
  /// @brief Setup class members (copy from \p source if provided)
  void initialize(const CompleteDisjunction* const source = NULL);
  /// @brief Get disjunction
  virtual DisjExitReason prepareDisjunction(const OsiSolverInterface* const si);
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
}; /* CompleteDisjunction */