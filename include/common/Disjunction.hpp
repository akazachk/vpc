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

/// Possible reasons we will exit disjunction generation
enum class DisjExitReason {
  SUCCESS_EXIT = 0, ///< everything looks okay
  OPTIMAL_SOLUTION_FOUND_EXIT, ///< integer-optimal solution found during creation of disjunction
  TOO_FEW_TERMS_EXIT, ///< disjunction has only one term
  NO_DISJUNCTION_EXIT, ///< no disjunction was able to be created
  UNKNOWN, ///< unknown exit reason
  NUM_DISJ_EXIT_REASONS ///< number of exit reasons
}; /* DisjExitReason */

/**
 * @brief Information on one disjunctive term
 *
 * Each disjunctive term is specified by a basis
 * as well as which bounds are changed to arrive at that term
 * (not counting the changed values at the root node that are common to all nodes)
 * We can also specify a set of inequalities to add
 */
class DisjunctiveTerm {
public:
  double obj = std::numeric_limits<double>::max(); ///< objective value
  std::vector<int> changed_var; ///< list of indices of variables with changed bounds
  std::vector<int> changed_bound; ///< for each var in #changed_var, which bound was changed: <= 0: lower bound, 1: upper bound
  std::vector<double> changed_value; ///< new value of the variable
  bool is_feasible; ///< is this term LP-feasible?
#ifdef USE_COIN
  CoinWarmStart* basis = NULL; ///< optional: saved basis for this term (to enable quick warm start)
  std::vector<OsiRowCut> ineqs; ///< optional: inequalities to add aside from changed bounds
#endif

  /// @brief Default constructor
  DisjunctiveTerm();

  /// @brief Copy constructor
  DisjunctiveTerm(const DisjunctiveTerm& source);

  /// @brief Destructor
  virtual ~DisjunctiveTerm();

  /// @brief Assignment operator
  DisjunctiveTerm& operator=(const DisjunctiveTerm& source);

  /// @brief Clone
  virtual DisjunctiveTerm* clone() const;

  /// @brief Initialize class members
  virtual void initialize(const DisjunctiveTerm* const source);

  /// @brief Clear memory
  void clear();
}; /* DisjunctiveTerm */

//<!--------------------------------------------->
/// @brief Generic abstract Disjunction class        
class Disjunction {
public:
  /// @name Required members
  ///@{
  int num_terms; ///< number of terms in disjunction (some may be infeasible)
  std::vector<DisjunctiveTerm> terms; ///< information about disjunctive terms, including optimal bases of each feasible term when available
  ///@}

  /// @name Optional members
  ///@{
  std::string name; ///< name (for printing)
  double best_obj; ///< value of term with best objective
  double worst_obj; ///< value of term with worst objective
  double integer_obj; ///< value of best integer-feasible solution
  std::vector<double> integer_sol; ///< integer-feasible solution

  /// indices of variables with changed bounds at root node
  std::vector<int> common_changed_var;
  /// which bound is changed: bound <= 0 is LB, bound = 1 is UB
  std::vector<int> common_changed_bound;
  /// value of changed variables at root node
  std::vector<double> common_changed_value;
  /// objective after applying common changes at root node
  double root_obj;
#ifdef USE_COIN
  /// Instead of specifying bounds that are changed at the root, you can provide inequalities valid throughout the tree
  std::vector<OsiRowCut> common_ineqs;
#endif
  ///@} // optional members

  /// @brief Default constructor
  Disjunction();

  /// @brief Copy constructor
  Disjunction(const Disjunction& source);

  /// @brief Destructor
  virtual ~Disjunction();

  /// @brief Assignment operator
  Disjunction& operator=(const Disjunction& source);

  /// @brief Clone
  virtual Disjunction* clone() const = 0;

  /// @brief For clearing things and setting up the disjunction as new
  virtual void setupAsNew();

  /// @brief Get disjunction
#ifdef USE_COIN
  virtual DisjExitReason prepareDisjunction(const OsiSolverInterface* const si) = 0;
#else
  virtual DisjExitReason prepareDisjunction() = 0;
#endif

  /// @brief Retrieve disjunctive term solver
#ifdef USE_COIN
  void getSolverForTerm(
    OsiSolverInterface*& termSolver,
    const int term_ind,
    const OsiSolverInterface* const solver,
    const bool shouldChangeBounds,
    const double DIFFEPS,
    FILE* logfile) const;
#else
  void getSolverForTerm(const int term_ind) const;
#endif

  /// @brief Set/update name of cut generating set (disjunction)
  static void setCgsName(std::string& cgsName, const std::string& disjTermName);
  /// @brief Set/update name of cut generating set (disjunction)
  static void setCgsName(std::string& cgsName, const int num_coeff,
      const int* const termIndices, const double* const termCoeff,
      const double termRHS, const bool append = false);
  /// @brief Set/update name of cut generating set (disjunction)
  static void setCgsName(std::string& cgsName, const int num_ineq_per_term,
      const std::vector<std::vector<int> >& termIndices,
      const std::vector<std::vector<double> >& termCoeff,
      const std::vector<double>& termRHS, const bool append = false);

  /// @brief Update best/worst obj
  void updateObjValue(const double obj);

protected:
  /// @brief Initialize class members (copy from \p source if provided)
  void initialize(const Disjunction* const source = NULL);
}; /* Disjunction */

/// @brief Set of disjunctions to test
class DisjunctionSet {
  public:
  /// @name Required members
  ///@{
    std::vector<Disjunction*> disjunctions; ///< disjunctions to test
  ///@}

  /// @name Optional members
  ///@{
    double best_obj; ///< value of term with best objective
    double worst_obj; ///< value of term with worst objective
    double integer_obj; ///< value of best integer-feasible solution
    std::vector<double> integer_sol; ///< integer-feasible solution
  ///@} // optional members

  /// @brief Default constructor
  DisjunctionSet();

  /// @brief Copy constructor
  DisjunctionSet(const DisjunctionSet& source);

  /// @brief Destructor
  virtual ~DisjunctionSet();

  /// @brief Assignment operator
  DisjunctionSet& operator=(const DisjunctionSet& source);

  /// @brief Clone
  virtual DisjunctionSet* clone() const;

  /// @brief For clearing things and setting up the disjunction as new
  virtual void setupAsNew();

  /// @brief Add disjunction to the set
  inline void addDisjunction(const Disjunction* const disj) {
    this->disjunctions.push_back(disj->clone());
  }

  /// @brief Number of disjunctions
  inline int size() const {
    return this->disjunctions.size();
  }

  /// @brief Return disjunction with index \p ind
  inline const Disjunction* const getDisjunction(const int ind) const {
    return (ind < this->size()) ? this->disjunctions[ind] : NULL;
  }

  /// @brief Prepare set of disjunctions
  #ifdef USE_COIN
  virtual DisjExitReason prepareDisjunction(const OsiSolverInterface* const si);
#else
  virtual DisjExitReason prepareDisjunction() = 0;
#endif

  /// @brief Update best/worst obj
  void updateObjValue(const double obj);

  protected:
  /// @brief Initialize class members (copy from \p source if provided)
  void initialize(const DisjunctionSet* const source = NULL);  
}; /* DisjunctionSet */
