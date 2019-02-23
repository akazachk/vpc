// Name:     Disjunction.cpp
// Author:   A. M. Kazachkov
// Date:     2018-02-22
//-----------------------------------------------------------------------------
#include "Disjunction.hpp"
#include <limits>

// Project files
#include "CglVPC.hpp" // ExitReason

/****************** PUBLIC  **********************/
/** Default constructor */
Disjunction::Disjunction() {
  initialize();
} /* default constructor */

/** Copy constructor */
Disjunction::Disjunction(const Disjunction& source) {
  initialize(&source);
} /* copy constructor */

/** Destructor */
Disjunction::~Disjunction() {
} /* destructor */

/** Assignment operator */
Disjunction& Disjunction::operator=(const Disjunction& source) {
  if (this != &source) {
    initialize(&source);
  }
  return *this;
} /* assignment operator */

ExitReason Disjunction::setBases(const OsiSolverInterface* const si,
    std::vector<int>& changed_var, std::vector<int>& changed_bound,
    std::vector<double>& changed_value) {
  ExitReason retval = ExitReason::UNKNOWN;
  return retval;
} /* setBases */

/****************** PROTECTED **********************/
void Disjunction::initialize(const Disjunction* const source) {
  if (source != NULL) {
    this->num_terms = source->num_terms;
    this->best_obj = source->best_obj;
    this->worst_obj = source->worst_obj;
    this->min_nb_obj_val = source->min_nb_obj_val;
    this->integer_obj = source->integer_obj;
    this->integer_sol = source->integer_sol;
    this->bases = source->bases;
    this->name = source->name;
  } else {
    this->num_terms = 0;
    this->best_obj = std::numeric_limits<double>::lowest();
    this->worst_obj = std::numeric_limits<double>::max();
    this->min_nb_obj_val = std::numeric_limits<double>::max();
    this->integer_obj = std::numeric_limits<double>::max();
    this->integer_sol.resize(0);
    this->bases.resize(0);
    this->name = "";
  }
} /* initialize */

void Disjunction::updateObjValue(const double objVal) {
  if (objVal < this->best_obj) {
    this->best_obj = objVal;
  }
  if (objVal > this->worst_obj) {
    this->worst_obj = objVal;
  }
} /* updateObjValue */

void Disjunction::updateNBObjValue(const double curr_nb_obj_val) {
  if (curr_nb_obj_val < this->min_nb_obj_val) {
    this->min_nb_obj_val = curr_nb_obj_val;
  }
} /* updateNBObjValue */
