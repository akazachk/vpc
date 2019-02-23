// Name:     Disjunction.hpp
// Author:   A. M. Kazachkov
// Date:     2018-02-22
//-----------------------------------------------------------------------------
#pragma once

/**********************************************/
/*  Generic abstract Disjunction class from   */
/*  which a point-ray collection can be       */
/*  constructed and VPCs can be generated     */
/**********************************************/

#include <string>
#include <vector>

// COIN-OR
#include <OsiSolverInterface.hpp>

// Project files
enum class ExitReason; // defined in CglVPC.hpp, which is included in the source file

// Abstract class
class Disjunction {
public:
  friend class CglVPC;
  int num_terms;
  double best_obj, worst_obj;
  double min_nb_obj_val;
  double integer_obj; // value of term with best and worst objective, and integer obj (if found)
  std::vector<double> integer_sol; // integer-feasible solution
  std::vector<CoinWarmStart*> bases; // optimal bases of each of the disjunctive terms
  std::string name;

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

  /** Get disjunction */
  ExitReason setBases(const OsiSolverInterface* const si,
      std::vector<int>& changed_var, std::vector<int>& changed_bound,
      std::vector<double>& changed_value);
  virtual ExitReason prepareDisjunction(OsiSolverInterface* const si) = 0;
protected:
  void initialize(const Disjunction* const source = NULL);
  void updateObjValue(const double obj);
  void updateNBObjValue(const double curr_nb_obj_val);
}; /* Disjunction */
