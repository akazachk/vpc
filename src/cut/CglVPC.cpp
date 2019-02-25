// Name:     CglVPC.cpp
// Author:   A. M. Kazachkov
// Date:     2018-Dec-24
//-----------------------------------------------------------------------------
#include "CglVPC.hpp"

#include <cmath> // abs, floor, ceil
#include <limits> // numeric_limits
#include <algorithm> // std::min_element, std::max_element

// COIN-OR files
#include <CbcModel.hpp>
#include <CglGMI.hpp>

// Project files
#include "BBHelper.hpp"
#include "VPCEventHandler.hpp"
#include "PRLP.hpp"
#include "SolverHelper.hpp"
#include "nbspace.hpp"
#include "utility.hpp"

// Various pre-built disjunction options
#include "Disjunction.hpp"
#include "PartialBBDisjunction.hpp"
#include "SplitDisjunction.hpp"

#ifdef TRACE
#include "debug.hpp"
#endif

const std::vector<std::string> CglVPC::VPCModeName {
  "PARTIAL_BB",
  "SPLITS",
  "CROSSES",
  "CUSTOM"
}; /* VPCModeName */
const std::vector<std::string> CglVPC::VPCTimeStatsName {
  "TOTAL_TIME",
  "INIT_SOLVE_TIME",
  "DISJ_SETUP_TIME",
  "GEN_DISJ_TIME",
  "PRLP_SETUP_TIME",
  "PRLP_SOLVE_TIME",
  "GEN_CUTS_TIME"
}; /* VPCTimeStatsName */
const std::vector<std::string> CglVPC::CutTypeName {
  "ONE_SIDED_CUT", "OPTIMALITY_CUT", "VPC"
}; /* CutTypeName */
const std::vector<std::string> CglVPC::CutHeuristicsName {
  "DUMMY_OBJ",
  "ALL_ONES",
  "CUT_VERTICES",
  "ITER_BILINEAR",
  "UNIT_VECTORS",
  "STRONG_LB",
  "TIGHT_POINTS",
  "TIGHT_RAYS",
  "TIGHT_POINTS2",
  "TIGHT_RAYS2",
  "ONE_SIDED"
}; /* CutHeuristicsName */
const std::vector<std::string> CglVPC::FailureTypeName {
  "ABANDONED",
  "BAD_DYNAMISM",
  "BAD_SUPPORT",
  "BAD_VIOLATION",
  "CUT_LIMIT",
  "DUAL_INFEASIBLE",
  "DUPLICATE_SIC",
  "DUPLICATE_VPC",
  "ITERATION_LIMIT",
  "ORTHOGONALITY_SIC",
  "ORTHOGONALITY_VPC",
  "PRIMAL_INFEASIBLE",
  "TIME_LIMIT",
  "NUMERICAL_ISSUES_WARNING",
  "PRIMAL_INFEASIBLE_NO_OBJ",
  "NUMERICAL_ISSUES_NO_OBJ",
  "UNKNOWN"
}; /* FailureTypeName */
const std::string CglVPC::time_T1 = "TIME_TYPE1_";
const std::string CglVPC::time_T2 = "TIME_TYPE2_";

/** Default constructor */
CglVPC::CglVPC() {
  initialize();
} /* default constructor */

/** Param constructor */
CglVPC::CglVPC(const VPCParameters& param) {
  initialize(NULL, &param);
} /* param constructor */

/** Copy constructor */
CglVPC::CglVPC(const CglVPC& source) : CglCutGenerator(source) {
  initialize(&source);
} /* copy constructor */

/** Destructor */
CglVPC::~CglVPC() {
  if (disj) { delete disj; }
} /* destructor */

/** Assignment operator */
CglVPC& CglVPC::operator=(const CglVPC& source) {
  if (this != &source) {
    initialize(&source);
  }
  return *this;
} /* assignment operator */

/** Clone */
CglCutGenerator* CglVPC::clone() const {
  return new CglVPC(*this);
} /* clone */

/** setParams */
void CglVPC::setParams(const VPCParameters& param) {
  this->params = param;
} /* setParams */

/**
 * @brief Generate VPCs from a partial branch-and-bound tree
 */
void CglVPC::generateCuts(const OsiSolverInterface& si, OsiCuts& cuts, const CglTreeInfo info) {
  ExitReason status = ExitReason::UNKNOWN;
  if (params.get(intParam::DISJ_TERMS) == 0) {
    status = ExitReason::NO_DISJUNCTION_EXIT;
    finish(status);
    return;
  }
  if (reachedTimeLimit(VPCTimeStatsName[static_cast<int>(VPCTimeStats::TOTAL_TIME)],
      params.get(TIMELIMIT))) {
    status = ExitReason::TIME_LIMIT_EXIT;
    finish(status);
    return;
  }

  mode = static_cast<VPCMode>(params.get(MODE));
  if (mode == VPCMode::PARTIAL_BB) {
    if (params.get(intParam::DISJ_TERMS) < 2) {
      status = ExitReason::NO_DISJUNCTION_EXIT;
      finish(status);
      return;
    }
    printf("\n## Starting VPC generation from partial branch-and-bound tree with up to %d disjunctive terms. ##\n", params.get(intParam::DISJ_TERMS));
    disj = new PartialBBDisjunction(this->params);
  }
  else if (mode == VPCMode::SPLITS) {
    disj = new SplitDisjunction(this->params);
  }
  else if (mode == VPCMode::CUSTOM) {
    if (!disj) {
      error_msg(errorstring, "Mode chosen is CUSTOM but no disjunction is set.\n");
      writeErrorToLog(errorstring, params.logfile);
      exit(1);
    }
  }
  else {
    error_msg(errorstring,
        "Mode that is chosen has not yet been implemented for VPC generation: %s.\n",
        VPCModeName[static_cast<int>(mode)].c_str());
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }
  disj->timer = &timer;

  // Solver can't be const because custom enableFactorization function might do initialSolve or resolve
  SolverInterface* solver;
  try {
    solver = const_cast<SolverInterface*>(dynamic_cast<const SolverInterface*>(&si));
  } catch (std::exception& e) {
    error_msg(errorstring, "Unable to cast as SolverInterface.\n");
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }
  if (!solver->isProvenOptimal()) {
    error_msg(errorstring, "CglVPC::generateCuts: Solver not proven optimal.\n");
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }

  // Set cut limit if needed
  if (params.get(CUTLIMIT) < 0) {
    params.set(CUTLIMIT,
        -1 * params.get(CUTLIMIT) * solver->getFractionalIndices().size());
  } else if (params.get(CUTLIMIT) == 0) {
    params.set(CUTLIMIT, std::numeric_limits<int>::max());
  }

  // Read opt value
  ip_opt = std::numeric_limits<double>::max();
  if (!params.get(stringParam::OPTFILE).empty()) {
#ifdef TRACE
    std::cout << "Reading objective information from \"" + params.get(stringParam::OPTFILE) + "\"" << std::endl;
#endif
    ip_opt = getObjValueFromFile(params.get(stringParam::OPTFILE), params.get(stringParam::FILENAME), params.logfile);
#ifdef TRACE
    std::cout << "Best known objective value is " << ip_opt << std::endl;
#endif
    if (isInfinity(ip_opt)) {
      warning_msg(warnstring, "Did not find objective value.\n");
    }
  }

  // Get disjunctive terms and obtain their optimal bases
  timer.start_timer(VPCTimeStatsName[static_cast<int>(VPCTimeStats::TOTAL_TIME)]);
  status = disj->prepareDisjunction(solver);
  if (status == ExitReason::PARTIAL_BB_OPTIMAL_SOLUTION_FOUND_EXIT) {
    warning_msg(warnstr,
        "An integer (optimal) solution was found prior while getting disjunction. "
        "We will generate between n and 2n cuts, restricting the value of each variable.\n");
    const double* solution = disj->integer_sol.data();
    for (int col = 0; col < static_cast<int>(disj->integer_sol.size()); col++) {
      const double val = solution[col];

      // Check which of the bounds needs to be fixed
      for (int b = 0; b < 2; b++) {
        if ((b == 0 && greaterThanVal(val, solver->getColLower()[col]))
            || (b == 1 && lessThanVal(val, solver->getColUpper()[col]))) {
          const double mult = (b == 0) ? 1. : -1.;
          const double el = mult * 1.;

          OsiRowCut currCut;
          currCut.setLb(mult * val);
          currCut.setRow(1, &col, &el, false);
          addCut(currCut, CutType::OPTIMALITY_CUT, cuts);
        }
      }
    } // iterate over columns and add optimality cut if needed
  } // exit out early if integer-optimal solution found
  if (status != ExitReason::SUCCESS_EXIT) {
    finish(status);
    return;
  }

  // Save the V-polyhedral relaxations of each optimal basis in the terms of the PRLP
  status = setupConstraints(solver, cuts);
  if (status != ExitReason::SUCCESS_EXIT) {
    finish(status);
    return;
  }

  // Generate cuts from the PRLP
  status = tryObjectives(cuts, solver, NULL, time_T1);
  if (status != ExitReason::SUCCESS_EXIT) {
    finish(status);
    return;
  }

  finish(status);
} /* generateCuts */

/****************** PROTECTED **********************/

void CglVPC::initialize(const CglVPC* const source, const VPCParameters* const param) {
  if (param != NULL)
    setParams(*param);
  if (source != NULL) {
    if (param == NULL)
      setParams(source->params);
    this->mode = source->mode;
    this->disj = source->disj;
    this->exitReason = source->exitReason;
    this->timer = source->timer;
    this->cutType = source->cutType;
    this->cutHeurVec = source->cutHeurVec;
    this->numCutsOfType = source->numCutsOfType;
    this->numObjFromHeur = source->numObjFromHeur;
    this->numCutsFromHeur = source->numCutsFromHeur;
    this->numFails = source->numFails;
    this->ip_opt = source->ip_opt;
    this->num_cgs = source->num_cgs;
    this->num_cgs_actually_used = source->num_cgs_actually_used;
    this->num_cgs_leading_to_cuts = source->num_cgs_leading_to_cuts;
    this->num_cuts = source->num_cuts;
    this->num_obj_tried = source->num_obj_tried;
    this->probData = source->probData;
    this->prlpData = source->prlpData;
  }
  else {
    this->mode = VPCMode::PARTIAL_BB;
    this->disj = NULL;
    this->exitReason = ExitReason::UNKNOWN;
    for (int t = 0; t < static_cast<int>(VPCTimeStats::NUM_TIME_STATS); t++) {
      timer.register_name(VPCTimeStatsName[t]);
    }
    for (int t = 0; t < static_cast<int>(CutHeuristics::NUM_CUT_HEUR); t++) {
      timer.register_name(CutHeuristicsName[t] + "_TIME");
    }
    this->cutType.resize(0);
    this->cutHeurVec.resize(0);
    this->numCutsOfType.resize(static_cast<int>(CutType::NUM_CUT_TYPES), 0);
    this->numObjFromHeur.resize(static_cast<int>(CutHeuristics::NUM_CUT_HEUR), 0);
    this->numCutsFromHeur.resize(static_cast<int>(CutHeuristics::NUM_CUT_HEUR), 0);
    this->numFails.resize(static_cast<int>(FailureType::NUM_FAILURES), 0);
    this->ip_opt = std::numeric_limits<double>::max();
    this->num_cgs = 1; // for now, always 1, unless we go back to doing multiple cgs
    this->num_cgs_actually_used = 0;
    this->num_cgs_leading_to_cuts = 0;
    this->num_cuts = 0;
    this->num_obj_tried = 0;
    this->probData.EPS = this->params.get(EPS);
  }
} /* initialize */

/**
 * Get problem data such as min/max coeff, problem-specific epsilon,
 * nonbasic variables, row in which each variable is basic, etc.
 */
void CglVPC::getProblemData(SolverInterface* const solver,
    ProblemData& probData, const ProblemData* const origProbData,
    const bool enable_factorization) {
  const int numCols = solver->getNumCols();
  const int numRows = solver->getNumRows();

  if (enable_factorization)
    enableFactorization(solver, params.get(doubleParam::EPS)); // this may change the solution slightly

  // Set min/max reference values based on problem data
  double minReferenceValue = 1, maxReferenceValue = 1;
  const double* elements = solver->getMatrixByCol()->getElements();
  probData.minAbsCoeff = *std::min_element(elements,
      elements + solver->getMatrixByCol()->getNumElements(),
      [](double i, double j) {return std::abs(i) < std::abs(j);});
  probData.maxAbsCoeff = *std::max_element(elements,
        elements + solver->getMatrixByCol()->getNumElements(),
        [](double i, double j) {return std::abs(i) < std::abs(j);});
  probData.minAbsCoeff = std::abs(probData.minAbsCoeff);
  probData.maxAbsCoeff = std::abs(probData.maxAbsCoeff);
  minReferenceValue =
      (probData.minAbsCoeff > 0) ?
          CoinMin(minReferenceValue, probData.minAbsCoeff) : minReferenceValue;
  maxReferenceValue = CoinMax(maxReferenceValue, probData.maxAbsCoeff);

  double lp_opt = solver->getObjValue();
  minReferenceValue =
      (lp_opt != 0) ?
          CoinMin(minReferenceValue, std::abs(lp_opt)) : minReferenceValue;
  maxReferenceValue = CoinMax(maxReferenceValue, std::abs(lp_opt));

  // Prepare data structures
  probData.num_cols = numCols;
  probData.lp_opt = lp_opt;
  probData.NBVarIndex.clear();
  probData.NBVarIndex.reserve(numCols);
  probData.rowOfVar.clear();
  probData.rowOfVar.resize(numCols + numRows, -1);
  probData.varBasicInRow.resize(numRows);
  solver->getBasics(&probData.varBasicInRow[0]); // get which variable is basic in each row
  if (origProbData) {
    probData.rowOfOrigNBVar.clear();
    probData.rowOfOrigNBVar.resize(numCols, -1);
  }

  // Get the basic and nonbasic original variable info
  // Note that fixed nonbasic variables might not need to be added to the nonbasic index vector
  // since they do not correspond to any rays...
  // but right now we typically assume we get n rays in some parts of the code
  for (int var = 0; var < numCols + numRows; var++) {
//  for (int col = 0; col < numCols; col++) {
    // Update min/max reference values
    if (var < numCols) {
      // Count how many non-inf lower and upper bounds
      double colLB, colUB;
      colLB = solver->getColLower()[var];
      colUB = solver->getColUpper()[var];
      if (colLB > -1 * solver->getInfinity() + params.get(doubleParam::EPS)) {
        minReferenceValue =
            (colLB != 0.) ?
                CoinMin(minReferenceValue, std::abs(colLB)) : minReferenceValue;
        maxReferenceValue = CoinMax(maxReferenceValue, std::abs(colLB));
      }

      if (colUB < solver->getInfinity() - params.get(doubleParam::EPS)) {
        minReferenceValue =
            (colUB != 0.) ?
                CoinMin(minReferenceValue, std::abs(colUB)) : minReferenceValue;
        maxReferenceValue = CoinMax(maxReferenceValue, std::abs(colUB));
      }
    } else {
      // Should only store the nonbasic slacks corresponding to inequalities since
      // the equality slack columns don't correspond to any rays
      // But, just as with fixed variables, we don't currently do this
      const int row = var - numCols;
      const double absrhs = std::abs(solver->getRightHandSide()[row]);
      minReferenceValue =
          (absrhs > 0) ? CoinMin(minReferenceValue, absrhs) : minReferenceValue;
      maxReferenceValue = CoinMax(maxReferenceValue, absrhs);

      // Assign var basic in this row to rowOfVar
      const int basic_var = probData.varBasicInRow[row];
      probData.rowOfVar[basic_var] = row;

      // If this basic_var was non-basic in the basis at v,
      // we need to add this row in the right spot in rowOfOrigNBVar
      if (origProbData) {
        const int NBIndexOfBasicVar = origProbData->getVarNBIndex(basic_var);
        if (NBIndexOfBasicVar >= 0) {
          probData.rowOfOrigNBVar[NBIndexOfBasicVar] = row;
        }
      }
    }

    if (!isBasicVar(solver, var)) {
      // Recall that rowOfVar stores -1 - nb index for nb variables
      // The -1 is to prevent the conflict of the 0th nb var and var basic in row 0
      probData.rowOfVar[var] -= probData.NBVarIndex.size();
      probData.NBVarIndex.push_back(var);
    } else {
      // Quick check that basic slack vars are basic in their row
      const int row = var - numCols;
      if (row >= 0) {
        if (probData.varBasicInRow[row] != var) {
          // Elsewhere we are using that each slack is basic in its own row,
          // so if this is not true, we will have to adjust
          // We use this, for example, in PCut, to calculate the right order
          // for the packed NB rays in genCornerNB
          error_msg(errstr,
              "Basic variable in row %d is variable %d, but it should be the slack on this row (%d).\n",
              row, probData.varBasicInRow[row], var);
          writeErrorToLog(errstr, params.logfile);
          exit(1);
        }
      }
    }
  } // loop over variables

  // TODO May also need to save rays of C1, where the coefficients of inv(B) * A
  // are sometimes negated because we have
  // \bar x = inv(B) * b - inv(B) * A * x_N
  const int numNB = probData.NBVarIndex.size();
  int tempIndex = 0;
//  const char* rowSense = solver->getRowSense();
  probData.NBReducedCost.resize(numNB);
  for (int j = 0; j < numNB; j++) {
    const int NBVar = probData.NBVarIndex[j];

    if (NBVar < numCols) {
      if (isNonBasicUBCol(solver, NBVar))
        probData.NBReducedCost[j] = -1.0 * solver->getReducedCost()[NBVar];
      else
        probData.NBReducedCost[j] = solver->getReducedCost()[NBVar];
    } else {
      tempIndex = NBVar - numCols;
      if (isNonBasicUBSlack(solver, tempIndex))
        probData.NBReducedCost[j] = solver->getRowPrice()[tempIndex];
      else
        probData.NBReducedCost[j] = -1.0 * solver->getRowPrice()[tempIndex];
    }
    if (lessThanVal(probData.NBReducedCost[j], 0.)) {
      if (lessThanVal(probData.NBReducedCost[j], 0., params.get(doubleConst::DIFFEPS))) {
        error_msg(errorstring,
            "Nonbasic reduced cost should be >= 0 in the complemented nonbasic space. "
            "However, for nb var %d (real var %d), it is %e.\n",
            j, NBVar, probData.NBReducedCost[j]);
        writeErrorToLog(errorstring, params.logfile);
        exit(1);
      } else {
        warning_msg(warnstring,
            "Nonbasic reduced cost should be >= 0 in the complemented nonbasic space. "
            "However, for nb var %d (real var %d), it is %f. Small enough error that we only send warning.\n",
            j, NBVar, probData.NBReducedCost[j]);
        this->numFails[static_cast<int>(FailureType::NUMERICAL_ISSUES_WARNING)]++;
      }
    }

//    if (isZero(nonBasicReducedCost[j], EPS)) {
//      numDualDegeneratePivots++;
//    }
  } // loop over nonbasic vars to collect reduced costs

  // Set data-specific epsilon
  probData.EPS = CoinMin(params.get(doubleParam::EPS), minReferenceValue / maxReferenceValue);

  if (enable_factorization)
    solver->disableFactorization();
} /* getProblemData */

ExitReason CglVPC::setupConstraints(const SolverInterface* const si, OsiCuts& cuts) {
  /***********************************************************************************
   * Change initial bounds
   ***********************************************************************************/
  // Make a copy of the solver to allow for fixing variables and changed bounds at root
  SolverInterface* vpcsolver = dynamic_cast<SolverInterface*>(si->clone());
  const int num_changed_bounds = this->disj->common_changed_var.size();
  int num_fixed = 0;
  for (int i = 0; i < num_changed_bounds; i++) {
    const int col = this->disj->common_changed_var[i];
    if (this->disj->common_changed_bound[i] <= 0) {
      vpcsolver->setColLower(col, this->disj->common_changed_value[i]);
    } else {
      vpcsolver->setColUpper(col, this->disj->common_changed_value[i]);
    }

    const double mult = (this->disj->common_changed_bound[i] <= 0) ? 1. : -1.;
    const double el = mult * 1.;
    const double val = this->disj->common_changed_value[i];
    OsiRowCut currCut;
    currCut.setLb(mult * val);
    currCut.setRow(1, &col, &el, false);
    addCut(currCut, CutType::ONE_SIDED_CUT, cuts);

    if (isVal(vpcsolver->getColLower()[col], vpcsolver->getColUpper()[col])) {
      num_fixed++;
    }
  } // account for bounds changed at root
#ifdef TRACE
  printf(
      "\n## Total number changed bounds: %d. Number fixed: %d. ##\n",
      num_changed_bounds, num_fixed);
#endif

  if (num_changed_bounds + num_fixed > 0) {
    vpcsolver->resolve();
    if (!checkSolverOptimality(vpcsolver, false)) {
      error_msg(errorstring,
          "Solver not proven optimal after updating bounds from root node.\n");
      writeErrorToLog(errorstring, params.logfile);
      exit(1);
    }
  }

  // Save problem data for nonbasic space usage and decide on problem-specific epsilon
#ifdef TRACE
  printf("\n## Saving root node optimal solution for nonbasic space usage ##\n");
#endif
  getProblemData(vpcsolver, this->probData);
#ifdef TRACE
  printf("Problem-specific epsilon set to %.10e\n", this->probData.EPS);
#endif

//  /***********************************************************************************
//   * Save objective
//   ***********************************************************************************/
//  // In order to calculate how good the points we find are
//  // with respect to the objective function (which should be in NB space when appropriate)
//  // Clearly a valid inequality is c^T x \ge c^\T v,
//  // where v is the LP optimum, since we are minimizing
//  // Either we have the objective, or we have the reduced costs
//  // of the nonbasic variables (which is exactly the objective
//  // in the nonbasic space).
//  // Recall that the reduced cost of slack variables is exactly
//  // the negative of the row price (dual variable) for that row.
//  OsiRowCut objCut;
//  if (inNBSpace) {
//    std::vector<int> indices;
//    std::vector<double> vals;
//    indices.reserve(dim);
//    vals.reserve(dim);
//    for (int i = 0; i < dim; i++) {
//      if (!isVal(probData.NBReducedCost[i], probData.EPS)) {
//        indices.push_back(i);
//        vals.push_back(probData.NBReducedCost[i]);
//      }
//    }
//    objCut.setRow(indices.size(), indices.data(), vals.data(), false);
//    objCut.setLb(0.0);
//  } else {
//    std::vector<int> indices;
//    std::vector<double> vals;
//    indices.reserve(dim);
//    vals.reserve(dim);
//    for (int i = 0; i < dim; i++) {
//      if (!isVal(vpcsolver->getObjCoefficients()[i], probData.EPS)) {
//        indices.push_back(i);
//        vals.push_back(vpcsolver->getObjCoefficients()[i]);
//      }
//    }
//    objCut.setRow(indices.size(), indices.data(), vals.data(), false);
//    objCut.setLb(vpcsolver->getObjValue());
//  }

  /***********************************************************************************
   * Set up intersection point and ray storage
   ***********************************************************************************/
  const int num_disj_terms = disj->num_terms;
  const int dim = si->getNumCols(); // NB: treating fixed vars as at bound
  std::vector < std::vector<std::vector<int> > > termIndices(num_disj_terms);
  std::vector < std::vector<std::vector<double> > > termCoeff(num_disj_terms);
  std::vector < std::vector<double> > termRHS(num_disj_terms);

  /***********************************************************************************
   * Get bases and generate VPCs
   ***********************************************************************************/
//  std::vector<ProblemData> disjProbData(num_disj_terms);

  int terms_added = -1;
  std::vector<bool> calcAndFeasTerm(num_disj_terms);

  // Start with the integer-feasible solution
  if (!(disj->integer_sol.empty())) {
    terms_added++;

    // Get solution and calculate slacks
    const double* sol = disj->integer_sol.data();
    std::vector<double> slack(vpcsolver->getNumRows(), 0.);
    const CoinPackedMatrix* mx = vpcsolver->getMatrixByRow();
    mx->times(sol, &slack[0]);
    for (int row_ind = 0; row_ind < vpcsolver->getNumRows(); row_ind++) {
      slack[row_ind] = vpcsolver->getRightHandSide()[row_ind] - slack[row_ind];
    }

    // Update the point
    CoinPackedVector point;
    double curr_nb_obj_val = 0.;
    setCompNBCoor(point, curr_nb_obj_val, params, sol, slack.data(), vpcsolver,
        probData.NBVarIndex, probData.NBReducedCost);
    const double beta = 1.0;
//        params.get(FLIP_BETA) >= 0 ? 1.0 : -1.0;
    prlpData.addConstraint(point, beta, terms_added, curr_nb_obj_val);
#ifdef TRACE
    printf("\n## Saving integer feasible solution as term %d. ##\n", terms_added);
#endif
    calcAndFeasTerm[terms_added] = true;

    // Update disjunctive bound info
    // This is here rather than in the Disjunction class,
    // because it is unclear whether, in that class,
    // the user computes with the variables changed at the root
    double objOffset = 0.;
    vpcsolver->getDblParam(OsiDblParam::OsiObjOffset, objOffset);
    const double objVal =
        dotProduct(vpcsolver->getObjCoefficients(), sol, dim) - objOffset;
    this->disj->updateObjValue(objVal);
    this->disj->updateNBObjValue(curr_nb_obj_val);
  } // integer-feasible solution

  // Now we handle the normal terms
  const int num_normal_terms = disj->terms.size();
  for (int tmp_ind = 0; tmp_ind < num_normal_terms; tmp_ind++) {
    terms_added++;
    SolverInterface* tmpSolver = dynamic_cast<SolverInterface*>(vpcsolver->clone());
    tmpSolver->disableFactorization();

    // Change bounds in the solver
    DisjunctiveTerm* term = &this->disj->terms[tmp_ind];
    const int curr_num_changed_bounds = term->changed_var.size();
    std::vector < std::vector<int> > commonTermIndices(curr_num_changed_bounds);
    std::vector < std::vector<double> > commonTermCoeff(curr_num_changed_bounds);
    std::vector<double> commonTermRHS(curr_num_changed_bounds);
    for (int i = 0; i < curr_num_changed_bounds; i++) {
      const int col = term->changed_var[i];
      const double coeff = (term->changed_bound[i] <= 0) ? 1. : -1.;
      const double val = term->changed_value[i];
      commonTermIndices[i].resize(1, col);
      commonTermCoeff[i].resize(1, coeff);
      commonTermRHS[i] = coeff * val;
      if (term->changed_bound[i] <= 0) {
        tmpSolver->setColLower(col, val);
      } else {
        tmpSolver->setColUpper(col, val);
      }
    }

    // Set the warm start
    if (!(tmpSolver->setWarmStart(term->basis))) {
      error_msg(errorstring,
          "Warm start information not accepted for leaf node %d/%d.\n",
          tmp_ind + 1, num_normal_terms);
      writeErrorToLog(errorstring, params.logfile);
      exit(1);
    }

    // Resolve and check the objective matches
#ifdef TRACE
    printf("\n## Solving for leaf node %d/%d. ##\n",
        tmp_ind + 1, num_normal_terms);
#endif
    tmpSolver->resolve();
    enableFactorization(tmpSolver, params.get(doubleParam::EPS)); // this may change the solution slightly
    calcAndFeasTerm[terms_added] = checkSolverOptimality(tmpSolver, true);

    // Sometimes we run into a few issues getting the ``right'' value
    if (!isVal(tmpSolver->getObjValue(), term->obj, params.get(doubleConst::DIFFEPS))) {
      tmpSolver->resolve();
    }
    if (!isVal(tmpSolver->getObjValue(), term->obj, params.get(doubleConst::DIFFEPS))) {
      double ratio = tmpSolver->getObjValue() / term->obj;
      if (ratio < 1.) {
        ratio = 1. / ratio;
      }
      // Allow it to be up to 3% off without causing an error
      if (greaterThanVal(ratio, 1.03)) {
        error_msg(errorstring,
            "Objective at leaf node %d/%d is incorrect. During BB, it was %s, now it is %s.\n",
            tmp_ind, num_normal_terms, stringValue(term->obj, "%1.3f").c_str(),
            stringValue(tmpSolver->getObjValue(), "%1.3f").c_str());
        writeErrorToLog(errorstring, params.logfile);
        exit(1);
      } else {
        warning_msg(warnstring,
            "Objective at disjunctive term %d/%d is incorrect. During BB, it was %s, now it is %s.\n",
            tmp_ind, num_normal_terms, stringValue(term->obj, "%1.3f").c_str(),
            stringValue(tmpSolver->getObjValue(), "%1.3f").c_str());
      }
#ifdef TRACE
      std::string commonName;
      Disjunction::setCgsName(commonName, curr_num_changed_bounds, commonTermIndices,
          commonTermCoeff, commonTermRHS, false);
      printf("Bounds changed: %s.\n", commonName.c_str());
#endif
    } // check that objective value matches

    // Possibly strengthen the term
    if (calcAndFeasTerm[terms_added] && params.get(intParam::STRENGTHEN) == 2) {
    // Add Gomory cuts on disjunctive term and resolve
      OsiCuts GMICs;
      CglGMI GMIGen;
      GMIGen.generateCuts(*tmpSolver, GMICs);
      tmpSolver->applyCuts(GMICs);
      tmpSolver->resolve();
      calcAndFeasTerm[terms_added] = checkSolverOptimality(tmpSolver, true);
    }
    if (calcAndFeasTerm[terms_added]) {
      // Update disjunctive bound info
      // This is here rather than in the Disjunction class,
      // because it is unclear whether, in that class,
      // the user computes with the variables changed at the root
      this->disj->updateObjValue(tmpSolver->getObjValue());
      this->disj->updateNBObjValue(tmpSolver->getObjValue() - vpcsolver->getObjValue());

      timer.register_name(CglVPC::time_T1 + std::to_string(terms_added));
      timer.register_name(CglVPC::time_T2 + std::to_string(terms_added));

      if (!reachedTimeLimit(CglVPC::time_T1 + "TOTAL", params.get(TIMELIMIT))) {
        timer.start_timer(CglVPC::time_T1 + "TOTAL");
        timer.start_timer(CglVPC::time_T1 + std::to_string(terms_added));

        // Get cobasis information and PR collection
        enableFactorization(tmpSolver, probData.EPS);
        ProblemData tmpData;
        getProblemData(tmpSolver, tmpData, &probData, false);
        genDepth1PRCollection(vpcsolver, tmpSolver,
            probData, tmpData, terms_added);

        timer.end_timer(CglVPC::time_T1 + "TOTAL");
        timer.end_timer(CglVPC::time_T1 + std::to_string(terms_added));
      }
    } // compute PR collection if the disj term is feasible
    if (tmpSolver) {
      delete tmpSolver;
    }
  } /* end loop over nodes on tree */

#ifdef TRACE
  printf(
      "\nFinished generating cgs from partial BB. Min obj val: Structural: %1.3f, NB: %1.3f.\n",
      this->disj->best_obj, this->disj->min_nb_obj_val);
#endif
  return ExitReason::SUCCESS_EXIT;
} /* setupConstraints */

/**
 * IN NON-BASIC SPACE
 * Get Point and Rays from corner polyhedron defined by current optimum at solver
 * Assumed to be optimal already
 */
void CglVPC::genDepth1PRCollection(const SolverInterface* const vpcsolver,
    const SolverInterface* const tmpSolver, const ProblemData& origProbData,
    const ProblemData& tmpProbData, const int term_ind) {
  // Ensure solver is optimal
  if (!tmpSolver->isProvenOptimal()) {
    error_msg(errstr, "Solver is not proven optimal.\n");
    writeErrorToLog(errstr, params.logfile);
    exit(1);
  }

  /***********************************************************************************
   * The first point is simply the first vertex
   ***********************************************************************************/
  // We need to calculate its coefficients in the NB space given in origProbData
  // That is, the NB space defined by v, the solution to the initial LP relaxation
  // *However* we need to work in the complemented NB space, where v is the origin
  // The point will be stored wrt origProbData cobasis,
  // but the cobasis of the point is given by nonBasicVarIndex
  CoinPackedVector optPoint;
  double curr_nb_obj_val = 0.;
  setCompNBCoorPoint(optPoint, curr_nb_obj_val, params, tmpSolver, vpcsolver,
      probData.NBVarIndex, probData.NBReducedCost);
  if (!isVal(curr_nb_obj_val, tmpSolver->getObjValue() - origProbData.lp_opt,
      params.get(doubleConst::DIFFEPS))) {
    error_msg(errorstring,
        "NB obj value is somehow incorrect. Calculated: %s, should be: %s.\n",
        stringValue(curr_nb_obj_val).c_str(),
        stringValue(tmpSolver->getObjValue() - origProbData.lp_opt).c_str());
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }
  const double beta = params.get(PRLP_BETA) >= 0 ? 1.0 : -1.0;
  if (beta > 0 && !isZero(curr_nb_obj_val)) { // check things are still okay after we scale; only applies when beta > 0, as otherwise the obj cut is not valid
    const double activity = curr_nb_obj_val
        / (tmpSolver->getObjValue() - origProbData.lp_opt);
    if (lessThanVal(activity, beta, params.get(EPS))) {
      numFails[static_cast<int>(FailureType::NUMERICAL_ISSUES_WARNING)]++;
      // Is it less by a little or by a lot?
      if (lessThanVal(activity, beta, params.get(doubleConst::DIFFEPS))) {
        error_msg(errorstring,
            "Term %d (point) does not satisfy the objective cut. Activity: %.8f. RHS: %.8f\n",
            term_ind, activity, beta);
        writeErrorToLog(errorstring, params.logfile);
        exit(1);
      } else {
        warning_msg(warnstring,
            "Term %d (point) does not satisfy the objective cut.\n\tActivity: %.8f\n\tRHS: %.8f.\n" "\tSmall enough violation that we only send a warning.\n",
            term_ind, activity, beta);
      }
    }
  }
  prlpData.addConstraint(optPoint, beta, term_ind, curr_nb_obj_val);

  /***********************************************************************************
   * We get the rays corresponding to non-basic variables
   ***********************************************************************************/
  // We now store the packed vector for the corner polyhedron rays.
  // For each non-basic var in nonBasicVarIndex, we get Binv * A.
  // Since we are looking at x_B = \bar a_0 - \bar A_N x_N,
  // the rays should be negative of Binv * A.
  // Exception is when the variable is non-basic at its upper bound.
  // In that case, ray is negated, which is equivalent to taking Binv * A as-is.
  // For free variables, we need both the ray from Binv * A and its negation.
  // For fixed variables, there is no need to get a ray, since we would not be able
  // to increase along that direction anyway.
  // Finally, considering the slack variables, if the row is
  //  ax <= b
  // then Clp stores a slack of the form
  //  ax + s = b, s >= 0
  // in which case we need to again take -Binv * A.
  // If the row is
  //  ax >= b
  // then Clp stores a slack of the form
  //  ax + s = b, s <= 0
  // in which case we already have the negative by using Binv * A.
  // This is because the Clp version is an upper-bounded variable,
  // so we have to negate.
  //
  // For the rows that are >=, if the variable basic in the row is slack,
  // then we have it be a_i - s = b_i, so we move everything to the other side.
  // I.e., we negate the *row* of Binv * A.
  // We iterate through each of the rays
  std::vector<double> currRay(tmpSolver->getNumRows());
  for (int ray_ind = 0; ray_ind < (int) tmpProbData.NBVarIndex.size(); ray_ind++) {
    const int NBVar = tmpProbData.NBVarIndex[ray_ind]; // This is the current ray

      tmpSolver->getBInvACol(NBVar, &(currRay[0]));
      if (isNonBasicLBVar(tmpSolver, NBVar)) {
        for (int row = 0; row < tmpSolver->getNumRows(); row++) {
          currRay[row] *= -1.;
        }
      }

    // We want to store this ray in the complemented original NB space
    CoinPackedVector currRayNB;
    curr_nb_obj_val = 0.;
    double scale = 1.;
    setCompNBCoorRay(currRayNB, currRay.data(), curr_nb_obj_val, scale, params,
        tmpSolver, vpcsolver, tmpProbData.rowOfOrigNBVar, probData.NBVarIndex,
        probData.NBReducedCost, NBVar);
    if (currRayNB.getNumElements() == 0) {
      continue;
    }
    // For debugging purposes, check obj dot product is correct
    const int numElem = currRayNB.getNumElements();
    const int* ind = currRayNB.getIndices();
    const double* vals = currRayNB.getElements();
    double calcRedCost = 0.;
    for (int el = 0; el < numElem; el++) {
      calcRedCost += vals[el] * origProbData.NBReducedCost[ind[el]];
    }
    calcRedCost /= scale;
    if (!isVal(curr_nb_obj_val, calcRedCost, params.get(doubleConst::DIFFEPS))) {
      error_msg(errorstring,
          "Calculated reduced cost is somehow incorrect. Calculated: %s, should be: %s.\n",
          stringValue(calcRedCost).c_str(),
          stringValue(curr_nb_obj_val).c_str());
      writeErrorToLog(errorstring, params.logfile);
      exit(1);
    }
    if (lessThanVal(calcRedCost, 0., params.get(EPS))) {
      numFails[static_cast<int>(FailureType::NUMERICAL_ISSUES_WARNING)]++;
      // Is it less by a little or by a lot?
      if (lessThanVal(calcRedCost, 0., params.get(doubleConst::DIFFEPS))) {
        error_msg(errorstring,
            "Term %d ray %d does not satisfy the objective cut. Activity: %.8f. RHS: %.8f\n",
            term_ind, ray_ind, calcRedCost, 0.);
        writeErrorToLog(errorstring, params.logfile);
        exit(1);
      } else {
        warning_msg(warnstring,
            "Term %d ray %d does not satisfy the objective cut.\n\tActivity: %.8f\n\tRHS: %.8f.\n" "\tSmall enough violation that we only send a warning.\n",
            term_ind, ray_ind, calcRedCost, 0.);
      }
    }
    prlpData.addConstraint(currRayNB, 0.0, term_ind, curr_nb_obj_val);
    // double objDepth = curr_nb_obj_val / currRayNB.twoNorm();
  } // iterate over nonbasic vars
} /* genDepth1PRCollection */

ExitReason CglVPC::tryObjectives(OsiCuts& cuts,
    const OsiSolverInterface* const origSolver, const OsiCuts* const structSICs,
    const std::string& timeName) {
  if (reachedTimeLimit(time_T1 + "TOTAL", params.get(TIMELIMIT))) {
    return ExitReason::TIME_LIMIT_EXIT;
  }

#ifdef TRACE
  printf("\n## CglVPC: Trying objectives. ##\n");
#endif

  PRLP* prlp = new PRLP(this);
  setLPSolverParameters(prlp, params.get(VERBOSITY));

  // We can scale the rhs for points by min_nb_obj_val
  const bool useScale = true && !isInfinity(std::abs(this->disj->min_nb_obj_val));
  const double scale = (!useScale || lessThanVal(this->disj->min_nb_obj_val, 1.)) ? 1. : this->disj->min_nb_obj_val;
  double beta = params.get(intParam::PRLP_BETA) >= 0 ? scale : -1. * scale;

//  std::vector<double> ortho;
  const bool isCutSolverPrimalFeas = prlp->setup(scale, false);
//  printf("# rows: %d\t # cols: %d\n", prlp->getNumRows(), prlp->getNumCols());
//  printf("# points: %d\t # rays: %d\n", prlp->numPoints, prlp->numRays);

  int curr_num_cuts = 0, init_num_obj = this->num_obj_tried;
  if (isCutSolverPrimalFeas) {
    curr_num_cuts += prlp->targetStrongAndDifferentCuts(beta, cuts, this->num_cuts,
        this->num_obj_tried, origSolver, structSICs, timeName,
        params.get(intConst::NB_SPACE));
  }

  const bool flipBeta = params.get(PRLP_BETA) > 0; // if the PRLP is infeasible, we will try flipping the beta
  if (flipBeta || greaterThanVal(this->disj->worst_obj, this->disj->best_obj)
      || greaterThanVal(this->disj->min_nb_obj_val, 0.)) {
//    for (int i = 0; i < this->num_cgs; i++) {
//      this->num_disj_terms += calcAndFeasFacet[i];
//    }
//    const int old_num_obj_tried = num_obj_tried;
//    const int old_num_cuts = cuts.sizeCuts();
//    const int old_num_primal_infeas = this->numFails[FailureType::PRIMAL_INFEASIBLE];
//    int num_generated = 0;
//
//    generateVPCsT1(structVPCs, probData, vpcTimeStats, structSICs,
//        interPtsAndRays, num_vpcs_per_cgs_T1[0], num_generated_vpcsT1, objCut,
//        calcAndFeasFacet, 0, 0, cgsName[0], max_num_cgs, dim, inNBSpace);
//
//    if (num_generated > old_num_cuts) {
//      num_cgs_actually_used++;
//      num_cgs_leading_to_cuts++;
//    } else if ((num_obj_tried > old_num_obj_tried)
//        && (old_num_primal_infeas == numFails[FailureType::PRIMAL_INFEASIBLE])) {
//      this->num_cgs_actually_used++;
//    } else {
//      // Clear out stuff because this cgs was not used?
//    }
  }

  if (prlp)
    delete prlp;

#ifdef TRACE
  printf("\n## CglVPC: Finished trying %d objectives. Generated %d cuts. ##\n",
      this->num_obj_tried - init_num_obj, curr_num_cuts);
#endif
  return ExitReason::SUCCESS_EXIT;
} /* tryObjectives */

void CglVPC::addCut(const OsiRowCut& cut, const CutType& type, OsiCuts& cuts) {
  cuts.insert(cut);
  cutType.push_back(type);
  numCutsOfType[static_cast<int>(type)]++;
  num_cuts++;
} /* addCut */

/**
 * @brief Universal way to check whether we reached the limit for the number of cuts for each split
 * This allows us to change between restricting number of cuts per split and total number of cuts easily
 */
int CglVPC::getCutLimit() const {
  // The cut limit is either across all cut-generating sets
  // or it is divided among the cut-generating sets (either as the limit / num_cgs, or as a fixed number per cgs)
  // If CUT_LIMIT = 0 => no cut limit
  // If CUT_LIMIT > 0 => absolute cut limit
  // If CUT_LIMIT < 0 => cut limit per cgs
  if (params.get(intParam::CUTLIMIT) == 0) {
    return std::numeric_limits<int>::max();
  } else if (params.get(intParam::CUTLIMIT) > 0) {
    return params.get(intParam::CUTLIMIT);
  } else {
    const int num_cgs = 1; //param.get(NUM_CGS);
    return (-1. * params.get(intParam::CUTLIMIT) / num_cgs);
  }
} /* getCutLimit */

/**
 * @brief Checks whether too many unsuccessful objective attempts have been made
 *
 * There are four types of checks. The first three have non-decreasing success requirements.
 * 1. Few cuts have been generated:
 *      default is FEW_CUTS = 1 cut, and the success threshold is at least 1 cut every 20 obj (fail ratio = .95).
 * 2. Many cuts have been generated:
 *      default is MANY_CUTS = .25 * CUT_LIMIT, and the success threshold is at least 1 cut every 10 obj (fail ratio = .90).
 * 3. Many obj have been tried, i.e., we have been trying for a long time so we better be successful super often:
 *      MANY_OBJ = max(FEW_CUTS / (1-few_cuts_fail_threshold), MANY_CUTS / (1-many_cuts_fail_threshold));
 *      default = max(20, 2.5 * CUT_LIMIT) and the success threshold is at least 1 cut every 5 obj (fail ratio = .80).
 * 4. Time is too long and we are not too successful:
 *       # obj tried >= MANY_OBJ && time >= 10 && average time / obj >= CUTSOLVER_TIMELIMIT + 1
 *       the success threshold is at least 1 cut every 3 obj
 *
 * Examples:
 * If the failure ratio is .85, then this will return true only if # obj >= MANY_OBJ.
 * This means that, as long as we have not tried too many times unsuccessfully,
 * then even if we have MANY_CUTS cuts, we would like more, and we feel that we are doing pretty well generating them.
 * (Of course, we have not hit the cut limit yet.)
 *
 * If your current failure ratio (# fails / # obj) is .93, then this will return true if # cuts >= MANY_CUTS, as .93 > .90.
 * Note that if # cuts < MANY_CUTS, but # obj >= MANY_OBJ,
 * then we would have failed earlier based on the more stringent limit for the MANY_OBJ case.
 *
 * If the failure ratio is .99, then we got here with # cuts < MANY_CUTS and # obj < MANY_OBJ.
 * Otherwise, if # cuts >= MANY_CUTS, then we would have hit the failure limit earlier (the first time it went above .90).
 * Similarly, if # obj >= MANY_OBJ, then we would have hit the limit even earlier.
 *
 * If the failure ratio is 1 (all failures), then we will reject if # obj > FEW_CUTS / (1.-few_cuts_fail_threshold).
 */
bool CglVPC::reachedFailureLimit(const int num_cuts, const int num_fails, //const double time,
    const double few_cuts_fail_threshold, const double many_cuts_fail_threshold,
    const double many_obj_fail_threshold, const double time_fail_threshold) const {
  const int num_obj_tried = num_cuts + num_fails;
  if (num_obj_tried == 0) {
    return false;
  }
  const int CUT_LIMIT = getCutLimit();
  const int FEW_CUTS = 1;
  const int NO_CUTS_OBJ_LIMIT = std::ceil(FEW_CUTS / (1. - few_cuts_fail_threshold));
  const int MANY_CUTS = std::ceil(.25 * CUT_LIMIT);
  const int MANY_OBJ = CoinMax(
      std::ceil(FEW_CUTS / (1. - few_cuts_fail_threshold)),
      std::ceil(MANY_CUTS / (1. - many_cuts_fail_threshold)));
  const double fail_ratio = (double) num_fails / num_obj_tried;
  bool reached_limit = false;
  if (num_obj_tried >= MANY_OBJ && greaterThanVal(fail_ratio, many_obj_fail_threshold)) {
    reached_limit = true;
  } else if (num_cuts >= MANY_CUTS && greaterThanVal(fail_ratio, many_cuts_fail_threshold)) {
    reached_limit = true;
  } else if (num_cuts >= FEW_CUTS && greaterThanVal(fail_ratio, few_cuts_fail_threshold)) {
    reached_limit = true;
  } else if (num_cuts < FEW_CUTS && num_obj_tried >= NO_CUTS_OBJ_LIMIT) {
    reached_limit = true;
  }
  const double time = timer.get_total_time(VPCTimeStatsName[static_cast<int>(VPCTimeStats::PRLP_SOLVE_TIME)]);
  if (!reached_limit && num_obj_tried >= MANY_OBJ && time > 10.
      && time / num_obj_tried >= params.get(PRLP_TIMELIMIT)) { // checks if average PRLP solve time is too high
    reached_limit = true;
  }
//  if (reached_limit) {
//    this->exitReason = ExitReason::FAIL_LIMIT_EXIT;
//  }
  return reached_limit;
} /* reachedFailureLimit */
