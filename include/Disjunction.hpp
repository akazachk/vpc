/**
 * @file Disjunction.hpp
 * @author A. M. Kazachkov
 * @date 2018-02-22
 */
#pragma once

#include <string>
#include <vector>

#ifdef USE_COIN
// COIN-OR
#include <OsiSolverInterface.hpp>
#include <OsiRowCut.hpp>
#endif

enum class DisjExitReason {
  SUCCESS_EXIT = 0,
  OPTIMAL_SOLUTION_FOUND_EXIT, /// integer-optimal solution found during creation of disjunction
  TOO_FEW_TERMS_EXIT, /// disjunction has only one term
  NO_DISJUNCTION_EXIT, /// no disjunction was able to be created
  UNKNOWN, /// unknown exit reason
  NUM_DISJ_EXIT_REASONS
}; /* DisjExitReason */

/**
 * Struct DisjunctiveTerm
 * Each disjunctive term is specified by a basis
 * as well as which bounds are changed to arrive at that term
 * (not counting the changed values at the root node that are common to all nodes)
 * We can also specify a set of inequalities to add
 */
class DisjunctiveTerm {
public:
  double obj = std::numeric_limits<double>::max();
  std::vector<int> changed_var;
  std::vector<int> changed_bound; // <= 0: lower bound, 1: upper bound
  std::vector<double> changed_value;
#ifdef USE_COIN
  CoinWarmStart* basis = NULL;
  std::vector<OsiRowCut> ineqs; // optional: inequalities to add aside from changed bounds
#endif

  /** Default constructor */
  DisjunctiveTerm();

  /** Copy constructor */
  DisjunctiveTerm(const DisjunctiveTerm& source);

  /** Destructor */
  virtual ~DisjunctiveTerm();

  /** Assignment operator */
  DisjunctiveTerm& operator=(const DisjunctiveTerm& source);

  /** Clone */
  virtual DisjunctiveTerm* clone() const;

  /** Initialize class members */
  virtual void initialize(const DisjunctiveTerm* const source);

  /** Clear memory */
  void clear();
}; /* DisjunctiveTerm */

/**********************************************/
/*  Generic abstract Disjunction class        */
/**********************************************/
class Disjunction {
public:
  // Required members
  int num_terms;
  std::vector<DisjunctiveTerm> terms; // optimal bases of parents of each of the disjunctive terms

  // Optional members
  std::string name;
  double best_obj, worst_obj; // value of term with best and worst objective
  double integer_obj;  // value of best integer-feasible solution
  std::vector<double> integer_sol; // integer-feasible solution

  /// Save changed variable bounds at root node (bound <= 0 is LB, bound = 1 is UB)
  std::vector<int> common_changed_var;
  std::vector<int> common_changed_bound;
  std::vector<double> common_changed_value;
#ifdef USE_COIN
  std::vector<OsiRowCut> common_ineqs;
#endif

  /** Default constructor */
  Disjunction();

  /** Copy constructor */
  Disjunction(const Disjunction& source);

  /** Destructor */
  virtual ~Disjunction();

  /** Assignment operator */
  Disjunction& operator=(const Disjunction& source);

  /** Clone */
  virtual Disjunction* clone() const = 0;

  /** For clearing things and setting up the disjunction as new */
  virtual void setupAsNew();

  /** Get disjunction */
#ifdef USE_COIN
  virtual DisjExitReason prepareDisjunction(const OsiSolverInterface* const si) = 0;
#else
  virtual DisjExitReason prepareDisjunction() = 0;
#endif

  /** Set/update name of cut generating set (disjunction) */
  static void setCgsName(std::string& cgsName, const std::string& disjTermName);
  static void setCgsName(std::string& cgsName, const int num_coeff,
      const int* const termIndices, const double* const termCoeff,
      const double termRHS, const bool append = false);
  static void setCgsName(std::string& cgsName, const int num_ineq_per_term,
      const std::vector<std::vector<int> >& termIndices,
      const std::vector<std::vector<double> >& termCoeff,
      const std::vector<double>& termRHS, const bool append = false);

  void updateObjValue(const double obj);

protected:
  void initialize(const Disjunction* const source = NULL);
}; /* Disjunction */
