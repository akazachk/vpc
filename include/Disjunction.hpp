// Name:     Disjunction.hpp
// Author:   A. M. Kazachkov
// Date:     2018-02-22
//-----------------------------------------------------------------------------
#pragma once

/**********************************************/
/*  Generic abstract Disjunction class        */
/*  from which a point-ray collection can be  */
/*  constructed and VPCs can be generated     */
/**********************************************/

#include <string>
#include <vector>

// COIN-OR
#include <OsiSolverInterface.hpp>

// Project files
#include "TimeStats.hpp"

// Project files
enum class ExitReason; // defined in CglVPC.hpp, which is included in the source file

/**
 * Struct DisjunctiveTerm
 * Each disjunctive term is specified by a basis
 * as well as which bounds are changed to arrive at that term
 * (not counting the changed values at the root node that are common to all nodes)
 * We can also specify a set of inequalities to add
 */
struct DisjunctiveTerm {
  CoinWarmStart* basis = NULL;
  double obj = std::numeric_limits<double>::max();
  std::vector<int> changed_var;
  std::vector<int> changed_bound; // <= 0: lower bound, 1: upper bound
  std::vector<double> changed_value;
  std::vector<OsiRowCut> ineqs; // optional: inequalities to add aside from changed bounds
};

/**
 * Abstract Class Disjunction
 */
class Disjunction {
public:
  friend class CglVPC;

  // Required members
  int num_terms;
  std::vector<DisjunctiveTerm> terms; // optimal bases of parents of each of the disjunctive terms

  // Optional members
  std::string name;
  double best_obj, worst_obj; // value of term with best and worst objective
  double integer_obj;  // value of best integer-feasible solution
  std::vector<double> integer_sol; // integer-feasible solution
  TimeStats* timer; // not owned by Disjunction

  // Save changed variable bounds at root node (bound <= 0 is LB, bound = 1 is UB)
  std::vector<int> common_changed_var;
  std::vector<int> common_changed_bound;
  std::vector<double> common_changed_value;
  std::vector<OsiRowCut> common_ineqs;

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
  virtual ExitReason prepareDisjunction(const OsiSolverInterface* const si) = 0;

  /** Set/update name of cut generating set (disjunction) */
  inline static void setCgsName(std::string& cgsName, const std::string& disjTermName) {
    if (disjTermName.empty()) {
      return;
    }
    if (!cgsName.empty()) {
      cgsName += " V ";
    }
    cgsName += "(";
    cgsName += disjTermName;
    cgsName += ")";
  } /* setCgsName (given disj term name) */

  inline static void setCgsName(std::string& cgsName, const int num_coeff,
      const int* const termIndices, const double* const termCoeff,
      const double termRHS, const bool append = false) {
    if (num_coeff == 0) {
      return;
    }
    if (!cgsName.empty()) {
      if (!append) {
        cgsName += " V ";
      } else {
        cgsName.resize(cgsName.size() - 1); // remove last ")"
        cgsName += "; ";
      }
    }
    cgsName += append ? "" : "(";
    for (int coeff_ind = 0; coeff_ind < num_coeff; coeff_ind++) {
      cgsName += (termCoeff[coeff_ind] > 0) ? "+" : "-";
      cgsName += "x";
      cgsName += std::to_string(termIndices[coeff_ind]);
    }
    cgsName += " >= ";
    cgsName += std::to_string((int) termRHS);
    cgsName += ")";
  } /* setCgsName (one ineq per term) */

  inline static void setCgsName(std::string& cgsName, const int num_ineq_per_term,
      const std::vector<std::vector<int> >& termIndices,
      const std::vector<std::vector<double> >& termCoeff,
      const std::vector<double>& termRHS, const bool append = false) {
    if (num_ineq_per_term == 0) {
      return;
    }
    if (!cgsName.empty()) {
      if (!append) {
        cgsName += " V ";
      } else {
        cgsName.resize(cgsName.size() - 1); // remove last ")"
        cgsName += "; ";
      }
    }
    cgsName += append ? "" : "(";
    for (int i = 0; i < num_ineq_per_term; i++) {
      setCgsName(cgsName, termIndices[i].size(), termIndices[i].data(),
          termCoeff[i].data(), termRHS[i], (i > 0));
    }
    cgsName += ")";
  } /* setCgsName */
protected:
  void initialize(const Disjunction* const source = NULL);
  void updateObjValue(const double obj);
//  void updateNBObjValue(const double curr_nb_obj_val);
}; /* Disjunction */
