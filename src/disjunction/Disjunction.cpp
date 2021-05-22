/**
 * @file Disjunction.cpp
 * @author A. M. Kazachkov
 * @date 2018-02-22
 */
#include "Disjunction.hpp"
#include <limits>
#include <cmath> // abs

#include "SolverHelper.hpp" // checkSolverOptimality
#include "utility.hpp" // error_msg, warning_msg, writeErrorToLogfile, stringValue

// Defining DisjunctiveTerm class
/** Default constructor */
DisjunctiveTerm::DisjunctiveTerm() {
  initialize(NULL);
} /* default constructor */

/** Copy constructor */
DisjunctiveTerm::DisjunctiveTerm(const DisjunctiveTerm& source) {
  initialize(&source);
} /* copy constructor */

/** Destructor */
DisjunctiveTerm::~DisjunctiveTerm() {
  clear();
} /* destructor */

/** Assignment operator */
DisjunctiveTerm& DisjunctiveTerm::operator=(const DisjunctiveTerm& source) {
  if (this != &source) {
    initialize(&source);
  }
  return *this;
} /* assignment operator */

/** Clone */
DisjunctiveTerm* DisjunctiveTerm::clone() const {
  return new DisjunctiveTerm(*this);
}

/**
 * Initialize members of DisjunctiveTerm
 * (copy from source when that is provided)
 */
void DisjunctiveTerm::initialize(const DisjunctiveTerm* const source) {
  if (source != NULL) {
    obj = source->obj;
    changed_var = source->changed_var;
    changed_bound = source->changed_bound;
    changed_value = source->changed_value;
#ifdef USE_COIN
    basis = source->basis->clone();
    ineqs = source->ineqs;
#endif
  } else {
    clear();
    obj = std::numeric_limits<double>::max();
    changed_var.resize(0);
    changed_bound.resize(0);
    changed_value.resize(0);
#ifdef USE_COIN
    ineqs.resize(0);
#endif
  }
} /* initialize */

/**
 * Clear memory
 */
void DisjunctiveTerm::clear() {
#ifdef USE_COIN
  if (basis) {
    delete basis;
    basis = NULL;
  }
#endif
} /* clear */

// Defining Disjunction class
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
  for (auto& term : this->terms) {
    term.clear();
  }
} /* destructor */

/** Assignment operator */
Disjunction& Disjunction::operator=(const Disjunction& source) {
  if (this != &source) {
    initialize(&source);
  }
  return *this;
} /* assignment operator */

/** Set up the disjunction class as new */
void Disjunction::setupAsNew() {
  this->name = "";
  this->best_obj = std::numeric_limits<double>::max();
  this->worst_obj = std::numeric_limits<double>::lowest();
  this->integer_obj = std::numeric_limits<double>::max();
  this->integer_sol.resize(0);
  this->common_changed_bound.resize(0);
  this->common_changed_value.resize(0);
  this->common_changed_var.resize(0);
#ifdef USE_COIN
  this->common_ineqs.resize(0);
#endif
  this->num_terms = 0;
  for (auto& term : this->terms)
    term.clear();
  this->terms.resize(0);
} /* setupAsNew */

#ifdef USE_COIN
/// Retrieve solver set up for disjunctive term (user responsibility to free this solver)
void Disjunction::getSolverForTerm(
    /// [out] solver that is created by this method, corresponding to the specified disjunctive term
    OsiSolverInterface*& termSolver,
    /// [in] term index
    const int term_ind,
    /// [in] original solver
    const OsiSolverInterface* const solver,
    /// [in] whether bounds should be changed or explicit constraints should be added
    const bool shouldChangeBounds,
    /// [in] tolerance for reconstructing disjunctive term solution
    const double DIFFEPS,
    /// [in] logfile for error printing
    FILE* logfile) const {
  termSolver = solver->clone();
  const DisjunctiveTerm* const term = &(this->terms[term_ind]);

  for (int i = 0; i < (int) this->common_changed_var.size(); i++) {
    const int col = this->common_changed_var[i];
    const double coeff = (this->common_changed_bound[i] <= 0) ? 1. : -1.;
    const double val = this->common_changed_value[i];
    if (shouldChangeBounds) {
      if (this->common_changed_bound[i] <= 0) {
        termSolver->setColLower(col, val);
      } else {
        termSolver->setColUpper(col, val);
      }
    } else {
      termSolver->addRow(1, &col, &coeff, coeff * val, solver->getInfinity());
    }
  } // loop over changed vars in common across terms

  for (int i = 0; i < (int) term->changed_var.size(); i++) {
    const int col = term->changed_var[i];
    const double coeff = (term->changed_bound[i] <= 0) ? 1. : -1.;
    const double val = term->changed_value[i];
    if (shouldChangeBounds) {
      if (term->changed_bound[i] <= 0) {
        termSolver->setColLower(col, val);
      } else {
        termSolver->setColUpper(col, val);
      }
    } else {
      termSolver->addRow(1, &col, &coeff, coeff * val, solver->getInfinity());
    }
  } // loop over changed vars in term

  const int curr_num_added_ineqs = term->ineqs.size();
  for (int i = 0; i < curr_num_added_ineqs; i++) {
    const OsiRowCut* currCut = &term->ineqs[i];
    termSolver->applyRowCuts(1, currCut); // hopefully this works
  }

  // Set the warm start
  if (term->basis && !(termSolver->setWarmStart(term->basis))) {
    error_msg(errorstring,
        "Warm start information not accepted for term %d/%d.\n", term_ind+1, this->num_terms);
    writeErrorToLog(errorstring, logfile);
    exit(1);
  }

  // Resolve and check the objective matches
#ifdef TRACE
  printf("\n## Solving for term %d/%d. ##\n", term_ind+1, this->num_terms);
#endif
  termSolver->resolve();
  const bool calcAndFeasTerm = checkSolverOptimality(termSolver, true);

  if (!calcAndFeasTerm) {
    printf("\n## Term %d/%d is not proven optimal. Exiting from this term. ##\n", term_ind+1, this->num_terms);
    delete termSolver;;
    return;
  }

  // Sometimes we run into a few issues getting the ``right'' value
  if (!isVal(termSolver->getObjValue(), term->obj, DIFFEPS)) {
    termSolver->resolve();
  }
  if (!isVal(termSolver->getObjValue(), term->obj, DIFFEPS)) {
    double ratio = termSolver->getObjValue() / term->obj;
    if (ratio < 1.) {
      ratio = 1. / ratio;
    }
    // Allow it to be up to 3% off without causing an error
    if (greaterThanVal(ratio, 1.03)) {
      error_msg(errorstring,
          "Objective at disjunctive term %d/%d is incorrect. Before, it was %s, now it is %s.\n",
          term_ind+1, this->num_terms, stringValue(term->obj, "%1.3f").c_str(),
          stringValue(termSolver->getObjValue(), "%1.3f").c_str());
      writeErrorToLog(errorstring, logfile);
      exit(1);
    } else {
      warning_msg(warnstring,
          "Objective at disjunctive term %d/%d is incorrect. Before, it was %s, now it is %s.\n",
          term_ind+1, this->num_terms, stringValue(term->obj, "%1.3f").c_str(),
          stringValue(termSolver->getObjValue(), "%1.3f").c_str());
    }
#ifdef TRACE
    std::string commonName;
    const int curr_num_changed_bounds = term->changed_var.size();
    std::vector < std::vector<int> > termIndices(curr_num_changed_bounds);
    std::vector < std::vector<double> > termCoeff(curr_num_changed_bounds);
    std::vector<double> termRHS(curr_num_changed_bounds);
    for (int i = 0; i < curr_num_changed_bounds; i++) {
      const int col = term->changed_var[i];
      const double coeff = (term->changed_bound[i] <= 0) ? 1. : -1.;
      const double val = term->changed_value[i];
      termIndices[i].resize(1, col);
      termCoeff[i].resize(1, coeff);
      termRHS[i] = coeff * val;
    }
    Disjunction::setCgsName(commonName, curr_num_changed_bounds, termIndices,
        termCoeff, termRHS, false);
    printf("Bounds changed: %s.\n", commonName.c_str());
#endif
  } // check that objective value matches
} /* getSolverForTerm */
#else
void Disjunction::getSolverForTerm(
    /// [in] term index
    const int term_ind) const {
  // This function left intentionally blank for now
} /* getSolverForTerm */
#endif // USE_COIN

/**
 * @details Append \p disjTermName to \p cgsName, adding a "V" if \p cgsName is nonempty (i.e., other terms are already present).
 */
void Disjunction::setCgsName(std::string& cgsName,
    const std::string& disjTermName) {
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

/**
 * @details If \p append is true, then we are adding an inequality to a particular term,
 * so the last ")" is replaced by a ";" and the inequality given by \p termIndices, \p termCoeff, \p termRHS
 * is added, followed by a ")".
 *
 * If \p append is false, then we are adding a new disjunctive term, so if \p cgsName is nonempty,
 * then "V" is inserted before the inequality is added within a new set of parentheses.
 */
void Disjunction::setCgsName(std::string& cgsName, const int num_coeff,
    const int* const termIndices, const double* const termCoeff,
    const double termRHS, const bool append) {
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
  const double EPS = 1e-7;
  for (int coeff_ind = 0; coeff_ind < num_coeff; coeff_ind++) {
    const double absCurrCoeff = std::abs(termCoeff[coeff_ind]);
    cgsName += (termCoeff[coeff_ind] > 0) ? "+" : "-";
    if (!(std::abs(absCurrCoeff - 1.) < EPS)) {
      if ((std::abs(absCurrCoeff - std::floor(absCurrCoeff)) < EPS)
          || (std::abs(absCurrCoeff - std::ceil(absCurrCoeff)) < EPS)) {
        cgsName += std::to_string(static_cast<int>(absCurrCoeff));
      } else {
        cgsName += std::to_string(absCurrCoeff);
      }
      cgsName += " ";
    }
    cgsName += "x";
    cgsName += std::to_string(termIndices[coeff_ind]);
  }
  cgsName += " >= ";

  if ((std::abs(termRHS - std::floor(termRHS)) < EPS)
      || (std::abs(termRHS - std::ceil(termRHS)) < EPS))
    cgsName += std::to_string(static_cast<int>(termRHS));
  else
    cgsName += std::to_string(termRHS);
  cgsName += ")";
} /* setCgsName (one ineq per term) */

/**
 * @details Adds a set of inequalties using #setCgsName(std::string&, const int, const int* const, const double* const, const double, const bool) called on each of the inequalities given by \p termIndices, \p termCoeff, and \p termRHS.
 */
void Disjunction::setCgsName(std::string& cgsName, const int num_ineq_per_term,
    const std::vector<std::vector<int> >& termIndices,
    const std::vector<std::vector<double> >& termCoeff,
    const std::vector<double>& termRHS, const bool append) {
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

void Disjunction::updateObjValue(const double objVal) {
  if (objVal < this->best_obj) {
    this->best_obj = objVal;
  }
  if (objVal > this->worst_obj) {
    this->worst_obj = objVal;
  }
} /* updateObjValue */

/****************** PROTECTED **********************/
void Disjunction::initialize(const Disjunction* const source) {
  if (source != NULL) {
    this->name = source->name;
    this->best_obj = source->best_obj;
    this->worst_obj = source->worst_obj;
    this->integer_obj = source->integer_obj;
    this->integer_sol = source->integer_sol;
    this->common_changed_bound = source->common_changed_bound;
    this->common_changed_value = source->common_changed_value;
    this->common_changed_var = source->common_changed_var;
#ifdef USE_COIN
    this->common_ineqs = source->common_ineqs;
#endif
    this->num_terms = source->num_terms;
    for (auto& term : this->terms)
      term.clear();
    this->terms = source->terms;
  } else {
    setupAsNew();
  }
} /* initialize */

