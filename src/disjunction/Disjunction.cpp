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

/** Set up the disjunction class as new (except the timer pointer) */
void Disjunction::setupAsNew() {
  this->name = "";
  this->best_obj = std::numeric_limits<double>::max();
  this->worst_obj = std::numeric_limits<double>::lowest();
  this->min_nb_obj_val = std::numeric_limits<double>::max();
  this->integer_obj = std::numeric_limits<double>::max();
  this->integer_sol.resize(0);
  this->common_changed_bound.resize(0);
  this->common_changed_value.resize(0);
  this->common_changed_var.resize(0);
  this->num_terms = 0;
  this->terms.resize(0);
} /* setupAsNew */

/****************** PROTECTED **********************/
void Disjunction::initialize(const Disjunction* const source) {
  if (source != NULL) {
    this->name = source->name;
    this->best_obj = source->best_obj;
    this->worst_obj = source->worst_obj;
    this->min_nb_obj_val = source->min_nb_obj_val;
    this->integer_obj = source->integer_obj;
    this->integer_sol = source->integer_sol;
    this->timer = source->timer;
    this->common_changed_bound = source->common_changed_bound;
    this->common_changed_value = source->common_changed_value;
    this->common_changed_var = source->common_changed_var;
    this->num_terms = source->num_terms;
    this->terms = source->terms;
  } else {
    this->timer = NULL;
    setupAsNew();
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
