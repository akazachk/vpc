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

// Project files
#include "utility.hpp"
#include "BBHelper.hpp"
#include "VPCEventHandler.hpp"

#ifdef TRACE
#include "debug.hpp"
#endif

/** Default constructor */
CglVPC::CglVPC() {
  initialize();
} /* default constructor */

/** Param constructor */
CglVPC::CglVPC(const VPCParameters& param) {
  initialize(NULL, &param);
} /* param constructor */

/** Copy constructor */
CglVPC::CglVPC(const CglVPC& source) {
  initialize(&source);
} /* copy constructor */

/** Destructor */
CglVPC::~CglVPC() {
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
  if (params.get(intParam::NUM_DISJ_TERMS) <= 1) {
    finish(ExitReason::SUCCESS_EXIT);
  }
#ifdef TRACE
  printf("\n############### Starting VPC generation from partial branch-and-bound tree. ###############\n");
#endif

//  auto start_chrono = std::chrono::high_resolution_clock::now();
  timer.start_timer(VPCTimeStatsName[VPCTimeStats::TOTAL_TIME]);
  if (reachedTimeLimit(VPCTimeStatsName[VPCTimeStats::TOTAL_TIME], params.get(TIMELIMIT))) {
    finish(ExitReason::TIME_LIMIT_EXIT);
    return;
  }

  OsiClpSolverInterface* solver = dynamic_cast<OsiClpSolverInterface*>(si.clone());
  if (!solver->isProvenOptimal()) {
    error_msg(errorstring, "CglVPC::generateCuts: Solver not proven optimal.\n");
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }


  // Read opt value
  ip_opt = std::numeric_limits<double>::max();
  if (!params.get(stringParam::OPTFILE).empty()) {
#ifdef TRACE
    std::cout << "Reading objective information from \"" + params.get(stringParam::OPTFILE) + "\"" << std::endl;
#endif
    ip_opt = getObjValueFromFile(params);
#ifdef TRACE
    std::cout << "Best known objective value is " << ip_opt << std::endl;
#endif
    if (isInfinity(ip_opt)) {
      warning_msg(warnstring, "Did not find objective value.\n");
    }
  }

  getProblemData(solver, this->probData);

  ExitReason status;

  // Choose disjunctive terms and obtain their optimal bases
  status = prepareDisjunction(solver, cuts);
  if (status != ExitReason::SUCCESS_EXIT) {
    if (solver)
      delete solver;
    finish(status);
  }

  // Save the V-polyhedral relaxations of each optimal basis in the terms of the PRLP
  setupPRLP(si);

  // Generate cuts from the PRLP
  tryObjectives(si);

  if (solver)
    delete solver;
  finish(status);
} /* generateCuts */

/****************** PROTECTED **********************/

void CglVPC::getProblemData(OsiSolverInterface* const solver, ProblemData& probData) {
  OsiClpSolverInterface* clpsolver = NULL;
  try {
    clpsolver = dynamic_cast<OsiClpSolverInterface*>(solver);
  } catch (std::exception& e) {
    error_msg(errorstring,
        "Unable to cast as OsiClpSolverInterface.\n");
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }
  enableFactorization(solver); // this may change the solution slightly

  int numCols = solver->getNumCols();
  int numRows = solver->getNumRows();

  // Set min/max reference values based on problem data
  double minReferenceValue = 1, maxReferenceValue = 1;
  const double* elements = solver->getMatrixByCol()->getElements();
  double minAbsCoeff = *std::min_element(elements,
      elements + solver->getMatrixByCol()->getNumElements(),
      [](double i, double j) {return std::abs(i) < std::abs(j);});
  double maxAbsCoeff = *std::max_element(elements,
        elements + solver->getMatrixByCol()->getNumElements(),
        [](double i, double j) {return std::abs(i) < std::abs(j);});
  minAbsCoeff = std::abs(minAbsCoeff);
  maxAbsCoeff = std::abs(maxAbsCoeff);
  minReferenceValue =
      (minAbsCoeff > 0) ?
          CoinMin(minReferenceValue, minAbsCoeff) : minReferenceValue;
  maxReferenceValue = CoinMax(maxReferenceValue, maxAbsCoeff);

  double lp_opt = solver->getObjValue();
  minReferenceValue =
      (lp_opt != 0) ?
          CoinMin(minReferenceValue, std::abs(lp_opt)) : minReferenceValue;
  maxReferenceValue = CoinMax(maxReferenceValue, std::abs(lp_opt));

//  // Get cstat and rstat
//  probData.cstat.resize(solver->getNumCols());
//  probData.rstat.resize(solver->getNumRows());
//  solver->getBasisStatus(&(probData.cstat[0]), &(probData.rstat[0]));

  // Prepare data structures
  probData.NBVarIndex.clear();
  probData.NBVarIndex.reserve(numCols);
  probData.rowOfVar.clear();
  probData.rowOfVar.resize(numCols + numRows, -1);
  probData.varBasicInRow.resize(numRows);
  solver->getBasics(&probData.varBasicInRow[0]); // get which variable is basic in each row
  probData.fractionalCore.clear();
//  basicOrigVarIndex.clear();
//  basicSlackVarIndex.clear();
//  basicOrigVarIndex.reserve(numCols);
//  basicSlackVarIndex.reserve(numRows);

  // Get the basic and nonbasic original variable info
  // Note that fixed nonbasic variables might not need to be added to the nonbasic index vector
  // since they do not correspond to any rays...
  // but right now we typically assume we get n rays in some parts of the code
  for (int col = 0; col < numCols; col++) {
    // Count how many non-inf lower and upper bounds
    double colLB, colUB;
    colLB = solver->getColLower()[col];
    colUB = solver->getColUpper()[col];
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

    if (isBasicVar(clpsolver, col)) {
//      basicOrigVarIndex.push_back(col);
      if (solver->isInteger(col)) {
        const double val = solver->getColSolution()[col];
        if (!isVal(val, std::floor(val),
            params.get(AWAY))
            && !isVal(val, std::ceil(val),
                params.get(AWAY))) {
          probData.fractionalCore.push_back(col);
        }
      }

//      // If it is at its lower or upper bound, it is primal degenerate
//      if (isVal(solver->getColSolution()[col], colLB)
//          || isVal(solver->getColSolution()[col], colUB)) {
//        numPrimalDegeneratePivots++;
//      }
    } else {
      // Recall that rowOfVar stores -1 - nb index for nb variables
      // The -1 is to prevent the conflict of the 0th nb var and var basic in row 0
      probData.rowOfVar[col] = -1 - (int) probData.NBVarIndex.size();
      probData.NBVarIndex.push_back(col);
    }
  } // loop over columns

  // Should only store the nonbasic slacks corresponding to inequalities since
  // the equality slack columns don't correspond to any rays
  const CoinPackedMatrix* mx = solver->getMatrixByRow();
//  const char* rowSense = solver->getRowSense();
  for (int row = 0; row < numRows; row++) {
    const double absrhs = std::abs(solver->getRightHandSide()[row]);
    minReferenceValue =
        (absrhs > 0) ? CoinMin(minReferenceValue, absrhs) : minReferenceValue;
    maxReferenceValue = CoinMax(maxReferenceValue, absrhs);
    //    slackVar.push_back(rhs[row] - rowActivity[row]);
    //    if (rowSense[row] == 'G')
    //      slackVar[row] *= -1.0;
    //    dualSoln.push_back(tempDualSoln[row]);
    //
    //if (basicSlack(row)) { // rstat[row] == 0 should also be basic? Seemingly not from ClpSimplex model status definitions.
    if (isBasicVar(clpsolver, row + numCols)) {
//      basicSlackVarIndex.push_back(row + numCols);

//      // If the slack in the row is zero, then this is a primal degenerate pivot
//      if (isVal(solver->getRowActivity()[row],
//          solver->getRightHandSide()[row])) {
//        numPrimalDegeneratePivots++;
//      }
    } else {
      probData.rowOfVar[row + numCols] = -1 - (int) probData.NBVarIndex.size();
      probData.NBVarIndex.push_back(row + numCols);
    }

    const int var = probData.varBasicInRow[row];
    const int currrow = var - numCols;
    if (isBasicVar(clpsolver, row + numCols) && currrow != row) {
      // Elsewhere we are using that each slack is basic in its own row,
      // so if this is not true, we will have to adjust
      // We use this, for example, in PCut, to calculate the right order
      // for the packed NB rays in genCornerNB
      error_msg(errstr,
          "Basic variable in row %d is variable %d, but it should be the slack on this row (%d).\n",
          row, var, row + numCols);
      writeErrorToLog(errstr, params.logfile);
      exit(1);
    }
    if (var < numCols) {
      probData.rowOfVar[var] = row;
    } else {
      probData.rowOfVar[var] = row;
    }
  } // loop over rows

  // TODO May also need ot save rays of C1, where the coefficients of inv(B) * A
  // are sometimes negated because we have
  // \bar x = inv(B) * b - inv(B) * A * x_N
  const int numNB = probData.NBVarIndex.size();
  int tempIndex = 0;
//  const char* rowSense = solver->getRowSense();
  probData.NBReducedCost.resize(numNB);
  for (int j = 0; j < numNB; j++) {
    const int NBVar = probData.NBVarIndex[j];

    if (NBVar < numCols) {
      if (isNonBasicUBCol(clpsolver, NBVar))
        probData.NBReducedCost[j] = -1.0 * solver->getReducedCost()[NBVar];
      else
        probData.NBReducedCost[j] = solver->getReducedCost()[NBVar];
    } else {
      tempIndex = NBVar - numCols;
      if (isNonBasicUBSlack(clpsolver, tempIndex))
        probData.NBReducedCost[j] = solver->getRowPrice()[tempIndex];
      else
        probData.NBReducedCost[j] = -1.0 * solver->getRowPrice()[tempIndex];
    }
    if (lessThanVal(probData.NBReducedCost[j], 0.)) {
      if (lessThanVal(probData.NBReducedCost[j], 0., 1e-3)) {
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

//  solver->disableFactorization();
} /* getProblemData */

CglVPC::ExitReason CglVPC::prepareDisjunction(OsiClpSolverInterface* solver,
    OsiCuts& cuts) {
  if (params.get(intParam::NUM_DISJ_TERMS) == 0) {
    return ExitReason::UNKNOWN;
  }

  // Set up solver
  OsiClpSolverInterface* BBSolver;
  try {
    BBSolver = dynamic_cast<OsiClpSolverInterface*>(solver->clone()); // this will fail if anything except Clp is used
  } catch (std::exception& e) {
    error_msg(errorstring, "Unable to clone solver as OsiClpSolverInterface.\n");
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }
  setupClpForCbc(BBSolver);

#ifdef TRACE
  printf("\n## Generating partial branch-and-bound tree. ##\n");
#endif
  CbcModel* cbc_model = new CbcModel;
  cbc_model->swapSolver(BBSolver);
  cbc_model->setModelOwnsSolver(true); // solver will be deleted with cbc object
  setIPSolverParameters(cbc_model);

  timer.start_timer(VPCTimeStatsName[VPCTimeStats::GEN_TREE_TIME]);
  int num_strong = params.get(intParam::PARTIAL_BB_NUM_STRONG);
  if (num_strong == -1) {
    num_strong = solver->getNumCols();
  }
  if (num_strong == -2) {
    num_strong = static_cast<int>(std::ceil(std::sqrt(solver->getNumCols())));
  }
  const int num_before_trusted = std::numeric_limits<int>::max(); // 10;
  generatePartialBBTree(params, cbc_model, solver,
      params.get(intParam::NUM_DISJ_TERMS), num_strong,
      num_before_trusted);
  timer.end_timer(VPCTimeStatsName[VPCTimeStats::GEN_TREE_TIME]);

  num_partial_bb_nodes = cbc_model->getNodeCount(); // save number of nodes looked at

  // Do stuff with the branch and bound tree
  VPCEventHandler* eventHandler;
  try {
    eventHandler = dynamic_cast<VPCEventHandler*>(cbc_model->getEventHandler());
  } catch (std::exception& e) {
    error_msg(errstr, "Could not get event handler.\n");
    writeErrorToLog(errstr, params.logfile);
    exit(1);
  }

  std::vector<NodeStatistics> stats = eventHandler->getStatsVector();
#ifdef TRACE
  printNodeStatistics(stats, false);
  if (eventHandler->getPrunedStatsVector().size() > 0) {
    printf("\n");
    printNodeStatistics(eventHandler->getPrunedStatsVector(), false);
  }
#endif

#ifdef TRACE
  if (std::abs(params.get(intParam::TEMP)) >= 10 && std::abs(params.get(intParam::TEMP)) <= 15) {
    generateTikzTreeString(eventHandler, params,
        params.get(intParam::PARTIAL_BB_STRATEGY),
        solver->getObjValue(), true);
    if (std::abs(params.get(intParam::TEMP)) == 14) {
      // Free
      if (cbc_model) {
        delete cbc_model;
        cbc_model = NULL;
      }
      return ExitReason::SUCCESS_EXIT;
    }
    if (std::abs(params.get(intParam::TEMP)) == 15) {
      exit(1); // this is during debug and does not free memory
    }
  }
#endif

  num_disj_terms = eventHandler->getNumLeafNodes()
      + eventHandler->isIntegerSolutionFound();
  num_pruned_nodes = eventHandler->getPrunedStatsVector().size()
      - eventHandler->isIntegerSolutionFound();
  if (cbc_model->status() == 0 || cbc_model->status() == 1
      || cbc_model->status() == 2 || num_disj_terms <= 1) {
    if (cbc_model->isProvenOptimal()) {
      warning_msg(warnstr,
          "An integer (optimal) solution was found prior while branching. " "We will generate between n and 2n cuts, restricting the value of each variable.\n");
      const double* solution = cbc_model->getColSolution();
      for (int col = 0; col < cbc_model->getNumCols(); col++) {
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
            addCut(currCut, CutType::ONE_SIDED_CUT, cuts);
          }
        }
      }

      // Set strong branching lb and ub to both be this value
      this->branching_lb = cbc_model->getObjValue();
      this->branching_ub = cbc_model->getObjValue();
    } else {
      warning_msg(warnstr,
          "Giving up on getting cuts from the partial branch-and-bound tree (bad status or too few terms). Model status is %d.\n",
          cbc_model->status());
    }

    // Free
    if (cbc_model) {
      delete cbc_model;
      cbc_model = NULL;
    }
    return ExitReason::PARTIAL_BB_INTEGER_SOLUTION_FOUND_EXIT;
  } /* exit out early if cbc_model status is 0 or insufficiently many disjunctive terms */

  /***********************************************************************************
   * Change initial bounds
   ***********************************************************************************/
  const int num_changed_bounds = stats[0].changed_var.size();
  int num_fixed = 0;
  for (int i = 0; i < num_changed_bounds; i++) {
    const int col = stats[0].changed_var[i];
    if (stats[0].changed_bound[i] <= 0) {
      solver->setColLower(col, stats[0].changed_value[i]);
    } else {
      solver->setColUpper(col, stats[0].changed_value[i]);
    }

    const double mult = (stats[0].changed_bound[i] <= 0) ? 1. : -1.;
    const double el = mult * 1.;
    const double val = stats[0].changed_value[i];
    OsiRowCut currCut;
    currCut.setLb(mult * val);
    currCut.setRow(1, &col, &el, false);
    addCut(currCut, CutType::ONE_SIDED_CUT, cuts);

    if (isVal(solver->getColLower()[col], solver->getColUpper()[col])) {
      num_fixed++;
    }
  }
#ifdef TRACE
  printf(
      "\n## Total number changed bounds: %d. Number fixed: %d. Number fixed in course of BB: %d. ##\n",
      num_changed_bounds, num_fixed, cbc_model->strongInfo()[1]);
#endif

  solver->resolve();

  if (!checkSolverOptimality(solver, false)) {
    error_msg(errorstring,
        "Solver not proven optimal after updating bounds from root node.\n");
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }

  // Decide if we are in the nonbasic space, based on option
  const bool inNBSpace = true; // TODO later
  const int dim = solver->getNumCols(); // TODO treating fixed vars as at bound

  /***********************************************************************************
   * Save objective
   ***********************************************************************************/
  // In order to calculate how good the points we find are
  // with respect to the objective function (which should be in NB space when appropriate)
  // Clearly a valid inequality is c^T x \ge c^\T v,
  // where v is the LP optimum, since we are minimizing
  // Either we have the objective, or we have the reduced costs
  // of the nonbasic variables (which is exactly the objective
  // in the nonbasic space).
  // Recall that the reduced cost of slack variables is exactly
  // the negative of the row price (dual variable) for that row.
  OsiRowCut objCut;
  if (inNBSpace) {
    std::vector<int> indices;
    std::vector<double> vals;
    indices.reserve(dim);
    vals.reserve(dim);
    for (int i = 0; i < dim; i++) {
      if (!isVal(probData.NBReducedCost[i], probData.EPS)) {
        indices.push_back(i);
        vals.push_back(probData.NBReducedCost[i]);
      }
    }
    objCut.setRow(indices.size(), indices.data(), vals.data(), false);
    objCut.setLb(0.0);
  } else {
    std::vector<int> indices;
    std::vector<double> vals;
    indices.reserve(dim);
    vals.reserve(dim);
    for (int i = 0; i < dim; i++) {
      if (!isVal(solver->getObjCoefficients()[i], probData.EPS)) {
        indices.push_back(i);
        vals.push_back(solver->getObjCoefficients()[i]);
      }
    }
    objCut.setRow(indices.size(), indices.data(), vals.data(), false);
    objCut.setLb(solver->getObjValue());
  }

  /***********************************************************************************
   * Set up intersection point and ray storage
   ***********************************************************************************/
//  std::vector < std::vector<std::vector<int> > > termIndices(num_disj_terms);
//  std::vector < std::vector<std::vector<double> >
//      > termCoeff(num_disj_terms);
//  std::vector < std::vector<double> > termRHS(num_disj_terms);
  std::vector<CoinPackedVector> PRLP_constraints;
  std::vector<double> PRLP_rhs;

  /***********************************************************************************
   * Get bases and generate VPCs
   ***********************************************************************************/
  std::vector<ProblemData> disjProbData(num_disj_terms);
//  std::vector < std::vector<int> > tmpVarBasicInRow(num_disj_terms);
//  std::vector < std::vector<int> > tmpBasicStructVar(num_disj_terms);
//  std::vector < std::vector<int> > tmpBasicSlackVar(num_disj_terms);
//  std::vector < std::vector<int> > rowOfTmpVar(num_disj_terms);
//  std::vector < std::vector<int> > rowOfOrigNBVar(num_disj_terms);
//  std::vector < std::vector<int> > tmpNonBasicVarIndex(num_disj_terms);
//  std::vector < std::vector<std::vector<double> >
//      > currBInvACol(num_disj_terms);
//  std::vector < std::vector<std::vector<double> > > currRay(num_disj_terms);

  int term_ind = -1;
  std::vector<std::string> cgsName(1);
  std::vector<bool> calcAndFeasFacet(num_disj_terms);
  this->branching_lb = std::numeric_limits<double>::max();
  this->branching_ub = std::numeric_limits<double>::lowest();
  this->min_nb_obj_val = std::numeric_limits<double>::max();

  // Start with the integer-feasible solution
  if (eventHandler->isIntegerSolutionFound()) {
    term_ind++;

    // Get solution and calculate slacks
    const double* sol = eventHandler->getIntegerFeasibleSolution();
    std::vector<double> slack(solver->getNumRows(), 0.);
    const CoinPackedMatrix* mx = solver->getMatrixByRow();
    mx->times(sol, &slack[0]);
    for (int row_ind = 0; row_ind < solver->getNumRows(); row_ind++) {
      slack[row_ind] = solver->getRightHandSide()[row_ind] - slack[row_ind];
    }

    // Update the point
    CoinPackedVector point;
    double curr_nb_obj_val = 0.;
    setCompNBCoor(point, curr_nb_obj_val, params, sol, slack.data(), solver,
        probData.NBVarIndex, probData.NBReducedCost);
    const double beta = 1.0;
//        params.get(FLIP_BETA) >= 0 ? 1.0 : -1.0;
//    interPtsAndRays[term_ind].addPointOrRayToIntersectionInfo(point, beta, 0,
//        curr_nb_obj_val);
    PRLP_constraints.push_back(point);
    PRLP_rhs.push_back(beta);
#ifdef TRACE
    printf("\n## Saving integer feasible solution as term %d. ##\n", term_ind);
#endif
    calcAndFeasFacet[term_ind] = true;
    std::string tmpname = "feasSol" + std::to_string(term_ind);
    setCgsName(cgsName[0], tmpname);
    double objOffset = 0.;
    solver->getDblParam(OsiDblParam::OsiObjOffset, objOffset);
    const double objVal = dotProduct(solver->getObjCoefficients(), sol,
        solver->getNumCols()) - objOffset;
    if (objVal < this->branching_lb) {
      this->branching_lb = objVal;
    }
    if (objVal > this->branching_ub) {
      this->branching_ub = objVal;
    }
    if (curr_nb_obj_val < this->min_nb_obj_val) {
      this->min_nb_obj_val = curr_nb_obj_val;
    }
  } /* integer-feasible solution */

//  // Now we handle the normal terms
//  for (int tmp_ind = 0; tmp_ind < eventHandler->getNumNodesOnTree();
//      tmp_ind++) {
//    //    printf("\n## DEBUG: Tmp index: %d ##\n", tmp_ind);
//    const int node_id = eventHandler->getNodeIndex(tmp_ind);
//    const int orig_node_id = stats[node_id].orig_id;
//
//    PointCutsSolverInterface* tmpSolverBase =
//        dynamic_cast<PointCutsSolverInterface*>(solver->clone());
//    tmpSolverBase->disableFactorization();
//
//    // Change bounds in the solver
//    const int curr_num_changed_bounds =
//        (orig_node_id == 0) ? 0 : stats[orig_node_id].changed_var.size();
//    std::vector < std::vector<int> > commonTermIndices(curr_num_changed_bounds);
//    std::vector < std::vector<double>
//        > commonTermCoeff(curr_num_changed_bounds);
//    std::vector<double> commonTermRHS(curr_num_changed_bounds);
//    for (int i = 0; i < curr_num_changed_bounds; i++) {
//      const int col = stats[orig_node_id].changed_var[i];
//      const double coeff =
//          (stats[orig_node_id].changed_bound[i] <= 0) ? 1. : -1.;
//      const double val = stats[orig_node_id].changed_value[i];
//      commonTermIndices[i].resize(1, col);
//      commonTermCoeff[i].resize(1, coeff);
//      commonTermRHS[i] = coeff * val;
//      //      tmpSolverBase->addRow(1, &col, &coeff, rhs,
//      //          tmpSolverBase->getInfinity());
//      if (stats[orig_node_id].changed_bound[i] <= 0) {
//        tmpSolverBase->setColLower(col, val);
//      } else {
//        tmpSolverBase->setColUpper(col, val);
//      }
//    }
//
//    // Set the warm start
//    if (!(tmpSolverBase->setWarmStart(eventHandler->getBasisForNode(tmp_ind)))) {
//      error_msg(errorstring,
//          "Warm start information not accepted for node %d.\n", node_id);
//      writeErrorToLog(errorstring, params.logfile);
//      exit(1);
//    }
//
//    // Resolve
//#ifdef TRACE
//    printf("\n## Solving for parent node %d/%d. ##\n", tmp_ind + 1,
//        eventHandler->getNumNodesOnTree());
//#endif
//    tmpSolverBase->resolve();
//    if (!checkSolverOptimality(tmpSolverBase, false)) {
//      error_msg(errorstring, "Solver not proven optimal for node %d.\n",
//          node_id);
//      writeErrorToLog(errorstring, params.logfile);
//      exit(1);
//    }
//    // Sometimes we run into a few issues getting the ``right'' value
//    if (!isVal(tmpSolverBase->getObjValue(), stats[node_id].obj, 1e-3)) {
//      tmpSolverBase->resolve();
//    }
//
//#ifdef TRACE
//    double curr_nb_obj_val = tmpSolverBase->getObjValue() - probData.LPopt; //getNBObjValue(tmpSolverBase, solver, probData);
//    printf("DEBUG: Node: %d .......... Obj val:%.3f .......... NB obj val: %.3f\n", node_id,
//        tmpSolverBase->getObjValue(), curr_nb_obj_val);
//    printNodeStatistics(stats[node_id], true);
//#endif
//    if (!isVal(tmpSolverBase->getObjValue(), stats[node_id].obj, 1e-3)) {
//      std::string commonName;
//      setCgsName(commonName, curr_num_changed_bounds, commonTermIndices,
//          commonTermCoeff, commonTermRHS, false);
//#ifdef TRACE
//      printf("Bounds changed: %s.\n", commonName.c_str());
//#endif
//      double ratio = tmpSolverBase->getObjValue() / stats[node_id].obj;
//      if (ratio < 1.) {
//        ratio = 1. / ratio;
//      }
//      // Allow it ot be up to 3% off without causing an error
//      if (greaterThanVal(ratio, 1.03)) {
//        error_msg(errorstring,
//            "Objective at parent node %d is incorrect. During BB, it was %s, now it is %s.\n",
//            node_id, stringValue(stats[node_id].obj, "%1.3f").c_str(),
//            stringValue(tmpSolverBase->getObjValue(), "%1.3f").c_str());
//        writeErrorToLog(errorstring, params.logfile);
//        exit(1);
//      } else {
//        warning_msg(warnstring,
//            "Objective at parent node %d is incorrect. During BB, it was %s, now it is %s.\n",
//            node_id, stringValue(stats[node_id].obj, "%1.3f").c_str(),
//            stringValue(tmpSolverBase->getObjValue(), "%1.3f").c_str());
//      }
//    }
//
//    // Now we go to the branch or branches left for this node
//    //    const int num_ineq_per_term = curr_num_changed_bounds + 1;  // leaf nodes have one extra changed
//    const int num_ineq_per_term = 1;  // leaf nodes have one extra changed
//    const int branching_index = stats[node_id].branch_index;
//    const int branching_variable = stats[node_id].variable;
//    int branching_way = stats[node_id].way;
//    int branching_value =
//        (branching_way <= 0) ? stats[node_id].ub : stats[node_id].lb;
//
//    // Set up termIndices, termCoeff, termRHS for this term
//    // First direction
//    term_ind++;
//    termIndices[term_ind].resize(num_ineq_per_term);
//    termCoeff[term_ind].resize(num_ineq_per_term);
//    termRHS[term_ind].resize(num_ineq_per_term);
//    for (int i = 0; i < num_ineq_per_term; i++) {
//      termIndices[term_ind][i].resize(1);
//      termCoeff[term_ind][i].resize(1);
//      //      if (i == num_ineq_per_disj - 1) {
//      termIndices[term_ind][i][0] = branching_variable;
//      termCoeff[term_ind][i][0] = (branching_way <= 0) ? -1. : 1.;
//      termRHS[term_ind][i] = termCoeff[term_ind][i][0] * branching_value;
//      //      } else {
//      //        termIndices[term_ind][i][0] = stats[node_ind].changed_var[i];
//      //        termCoeff[term_ind][i][0] = (stats[node_ind].changed_way[i] <= 0) ? -1. : 1.;
//      //        termRHS[term_ind][i] = termCoeff[term_ind][i][0]
//      //            * stats[node_ind].changed_bound[i];
//      //      }
//    }
//
//    PointCutsSolverInterface* tmpSolver =
//        dynamic_cast<PointCutsSolverInterface*>(tmpSolverBase->clone());
//#ifdef TRACE
//    printf("\n## Solving first child of parent %d/%d for term %d. ##\n",
//        tmp_ind + 1, eventHandler->getNumNodesOnTree(), term_ind);
//#endif
//    calcAndFeasFacet[term_ind] =
//        (param.getTEMP() >= 100) ?
//            fixDisjTerm(tmpSolver, termIndices[term_ind], termCoeff[term_ind],
//                termRHS[term_ind]) :
//            fixDisjTermUsingBounds(tmpSolver, termIndices[term_ind],
//                termCoeff[term_ind], termRHS[term_ind]);
//
//    // Possibly strengthen
//    if (calcAndFeasFacet[term_ind]
//        && param.getParamVal(ParamIndices::STRENGTHEN_PARAM_IND) == 2) {
//      // Add Gomory cuts on disjunctive term and resolve
//      OsiCuts GMICs;
//      CglGMI GMIGen(param);
//      GMIGen.generateCuts(*tmpSolver, GMICs);
//      tmpSolver->applyCuts(GMICs);
//      tmpSolver->resolve();
//      calcAndFeasFacet[term_ind] = checkSolverOptimality(tmpSolver, true);
//    }
//
//    if (calcAndFeasFacet[term_ind]) {
//      setCgsName(cgsName[0], num_ineq_per_term, termIndices[term_ind],
//          termCoeff[term_ind], termRHS[term_ind]);
//      setCgsName(cgsName[0], curr_num_changed_bounds, commonTermIndices,
//          commonTermCoeff, commonTermRHS, true);
//      if (stats[node_id].depth + 1 < GlobalVariables::minNodeDepth) {
//        GlobalVariables::minNodeDepth = stats[node_id].depth + 1;
//      }
//      if (stats[node_id].depth + 1 > GlobalVariables::maxNodeDepth) {
//        GlobalVariables::maxNodeDepth = stats[node_id].depth + 1;
//      }
//
//      // Update strong branching info
//      if (tmpSolver->getObjValue() < this->strong_branching_lb) {
//        this->strong_branching_lb = tmpSolver->getObjValue();
//      }
//      if (tmpSolver->getObjValue() > this->strong_branching_ub) {
//        this->strong_branching_ub = tmpSolver->getObjValue();
//      }
//
//      // Get cobasis information and PR collection
//      enableFactorization(tmpSolver);
//      //      // After enabling factorization, things sometimes look different and it is worth resolving (e.g., bc1 with num_strong = 5, row price on row 458)
//      //      if ((tmpSolver->getModelPtr()->sumPrimalInfeasibilities() > param.getEPS())
//      //          || (tmpSolver->getModelPtr()->sumDualInfeasibilities()
//      //              > param.getEPS())) {
//      //        tmpSolver->initialSolve();
//      //        tmpSolver->disableFactorization();
//      //        tmpSolver->enableFactorization();
//      //      }
//      identifyCobasis(tmpVarBasicInRow[term_ind], tmpBasicStructVar[term_ind],
//          tmpBasicSlackVar[term_ind], rowOfTmpVar[term_ind],
//          rowOfOrigNBVar[term_ind], tmpNonBasicVarIndex[term_ind],
//          currBInvACol[term_ind], currRay[term_ind], tmpSolver, probData,
//          termIndices[term_ind][0][0], false, true);
//
//      vpcTimeStats.register_name(time_T1 + std::to_string(term_ind));
//      vpcTimeStats.register_name(time_T2 + std::to_string(term_ind));
//
//      if (!reachedTimeLimit(vpcTimeStats, time_T1 + "TOTAL",
//          param.getTIMELIMIT())) {
//        vpcTimeStats.start_timer(time_T1 + "TOTAL");
//        vpcTimeStats.start_timer(time_T1 + std::to_string(term_ind));
//
//        genDepth1PRCollection(interPtsAndRays[term_ind], tmpSolver, probData,
//            termIndices[term_ind], tmpVarBasicInRow[term_ind],
//            tmpBasicStructVar[term_ind], tmpBasicSlackVar[term_ind],
//            rowOfTmpVar[term_ind], rowOfOrigNBVar[term_ind],
//            tmpNonBasicVarIndex[term_ind], currRay[term_ind], term_ind, false,
//            0, max_num_cgs);
//
//        vpcTimeStats.end_timer(time_T1 + "TOTAL");
//        vpcTimeStats.end_timer(time_T1 + std::to_string(term_ind));
//
//        // Check that the points and rays we have generated satisfy the objective cut
//        // This is now done when they are generated
//        //        checkPRCollectionForFeasibility(term_ind, interPtsAndRays[term_ind],
//        //            curr_nb_obj_val, objCut);
//      }
//    }
//
//    if (tmpSolver) {
//      delete tmpSolver;
//    }
//
//    // Should we compute a second branch?
//    if (branching_index == 0) {
//      // Other direction
//      branching_way = (stats[node_id].way <= 0) ? 1 : 0;
//      branching_value =
//          (branching_way <= 0) ? branching_value - 1 : branching_value + 1;
//
//      term_ind++;
//      termIndices[term_ind].resize(num_ineq_per_term);
//      termCoeff[term_ind].resize(num_ineq_per_term);
//      termRHS[term_ind].resize(num_ineq_per_term);
//      for (int i = 0; i < num_ineq_per_term; i++) {
//        termIndices[term_ind][i].resize(1);
//        termCoeff[term_ind][i].resize(1);
//        //      if (i == num_ineq_per_disj - 1) {
//        termIndices[term_ind][i][0] = branching_variable;
//        termCoeff[term_ind][i][0] = (branching_way <= 0) ? -1. : 1.;
//        termRHS[term_ind][i] = termCoeff[term_ind][i][0] * branching_value;
//        //      } else {
//        //        termIndices[term_ind][i][0] = stats[node_ind].changed_var[i];
//        //        termCoeff[term_ind][i][0] = (stats[node_ind].changed_way[i] <= 0) ? -1. : 1.;
//        //        termRHS[term_ind][i] = termCoeff[term_ind][i][0]
//        //            * stats[node_ind].changed_bound[i];
//        //      }
//      }
//
//      // Change back to original state then add this term
//      tmpSolver =
//          dynamic_cast<PointCutsSolverInterface*>(tmpSolverBase->clone());
//#ifdef TRACE
//      printf("\n## Solving second child of parent %d/%d for term %d. ##\n",
//          tmp_ind + 1, eventHandler->getNumNodesOnTree(), term_ind);
//#endif
//      calcAndFeasFacet[term_ind] =
//          (param.getTEMP() >= 100) ?
//              fixDisjTerm(tmpSolver, termIndices[term_ind], termCoeff[term_ind],
//                  termRHS[term_ind]) :
//              fixDisjTermUsingBounds(tmpSolver, termIndices[term_ind],
//                  termCoeff[term_ind], termRHS[term_ind]);
//
//      // Possibly strengthen
//      if (calcAndFeasFacet[term_ind]
//          && param.getParamVal(ParamIndices::STRENGTHEN_PARAM_IND) == 2) {
//        // Add Gomory cuts on disjunctive term and resolve
//        OsiCuts GMICs;
//        CglGMI GMIGen(param);
//        GMIGen.generateCuts(*tmpSolver, GMICs);
//        tmpSolver->applyCuts(GMICs);
//        tmpSolver->resolve();
//        calcAndFeasFacet[term_ind] = checkSolverOptimality(tmpSolver, true);
//      }
//
//      if (calcAndFeasFacet[term_ind]) {
//        setCgsName(cgsName[0], num_ineq_per_term, termIndices[term_ind],
//            termCoeff[term_ind], termRHS[term_ind]);
//        setCgsName(cgsName[0], curr_num_changed_bounds, commonTermIndices,
//            commonTermCoeff, commonTermRHS, true);
//        if (stats[node_id].depth + 1 < GlobalVariables::minNodeDepth) {
//          GlobalVariables::minNodeDepth = stats[node_id].depth + 1;
//        }
//        if (stats[node_id].depth + 1 > GlobalVariables::maxNodeDepth) {
//          GlobalVariables::maxNodeDepth = stats[node_id].depth + 1;
//        }
//
//        // Update strong branching info
//        if (tmpSolver->getObjValue() < this->strong_branching_lb) {
//          this->strong_branching_lb = tmpSolver->getObjValue();
//        }
//        if (tmpSolver->getObjValue() > this->strong_branching_ub) {
//          this->strong_branching_ub = tmpSolver->getObjValue();
//        }
//
//        // Get cobasis information and PR collection
//        enableFactorization(tmpSolver);
//        //        // After enabling factorization, things sometimes look different and it is worth resolving (e.g., bc1 with num_strong = 5, row price on row 458)
//        //        if ((tmpSolver->getModelPtr()->sumPrimalInfeasibilities()
//        //            > param.getEPS())
//        //            || (tmpSolver->getModelPtr()->sumDualInfeasibilities()
//        //                > param.getEPS())) {
//        //          tmpSolver->initialSolve();
//        //          tmpSolver->disableFactorization();
//        //          tmpSolver->enableFactorization();
//        //        }
//        identifyCobasis(tmpVarBasicInRow[term_ind], tmpBasicStructVar[term_ind],
//            tmpBasicSlackVar[term_ind], rowOfTmpVar[term_ind],
//            rowOfOrigNBVar[term_ind], tmpNonBasicVarIndex[term_ind],
//            currBInvACol[term_ind], currRay[term_ind], tmpSolver, probData,
//            termIndices[term_ind][0][0], false, true);
//
//        vpcTimeStats.register_name(time_T1 + std::to_string(term_ind));
//        vpcTimeStats.register_name(time_T2 + std::to_string(term_ind));
//
//        if (!reachedTimeLimit(vpcTimeStats, time_T1 + "TOTAL",
//            param.getTIMELIMIT())) {
//          vpcTimeStats.start_timer(time_T1 + "TOTAL");
//          vpcTimeStats.start_timer(time_T1 + std::to_string(term_ind));
//
//          genDepth1PRCollection(interPtsAndRays[term_ind], tmpSolver, probData,
//              termIndices[term_ind], tmpVarBasicInRow[term_ind],
//              tmpBasicStructVar[term_ind], tmpBasicSlackVar[term_ind],
//              rowOfTmpVar[term_ind], rowOfOrigNBVar[term_ind],
//              tmpNonBasicVarIndex[term_ind], currRay[term_ind], term_ind, false,
//              0, max_num_cgs);
//
//          vpcTimeStats.end_timer(time_T1 + "TOTAL");
//          vpcTimeStats.end_timer(time_T1 + std::to_string(term_ind));
//
//          // Check that the points and rays we have generated satisfy the objective cut
//          // This is now done when they are generated
//          //          checkPRCollectionForFeasibility(term_ind, interPtsAndRays[term_ind],
//          //              curr_nb_obj_val, objCut);
//        }
//      }
//      if (tmpSolver) {
//        delete tmpSolver;
//      }
//    } /* end second branch computation */
//
//    if (tmpSolverBase) {
//      delete tmpSolverBase;
//    }
//  } /* end loop over nodes on tree */

#ifdef TRACE
  printf(
      "\nFinished generating cgs from partial BB. Min obj val: Structural: %1.3f, NB: %1.3f.\n",
      this->branching_lb, this->min_nb_obj_val);
#endif

  if (cbc_model)
    delete cbc_model; // deletes BBSolver
  return ExitReason::SUCCESS_EXIT;
} /* prepareDisjunction */

void CglVPC::setupPRLP(const OsiSolverInterface& si) {

} /* setupPRLP */

void CglVPC::tryObjectives(const OsiSolverInterface& PRLP) {

} /* tryObjectives */

void CglVPC::addCut(const OsiRowCut& cut, const CutType& type, OsiCuts& cuts) {
  cuts.insert(cut);
  cutType.push_back(type);
  numCutsOfType[static_cast<int>(type)]++;
  num_cuts++;
} /* addCut */

void CglVPC::initialize(const CglVPC* const source, const VPCParameters* const param) {
  if (param != NULL)
    setParams(*param);
  if (source != NULL) {
    if (param == NULL)
      this->params = source->params;
    this->exitReason = source->exitReason;
    this->timer = source->timer;
    this->cutType = source->cutType;
    this->numCutsOfType = source->numCutsOfType;
    this->numFails = source->numFails;
    this->branching_lb = source->branching_lb;
    this->branching_ub = source->branching_ub;
    this->min_nb_obj_val = source->min_nb_obj_val;
    this->ip_opt = source->ip_opt;
    this->num_cgs = source->num_cgs;
    this->num_cuts = source->num_cuts;
    this->num_disj_terms = source->num_disj_terms;
    this->num_partial_bb_nodes = source->num_partial_bb_nodes;
    this->num_pruned_nodes = source->num_pruned_nodes;
    this->num_obj_tried = source->num_obj_tried;
    this->probData = source->probData;
  }
  else {
//    if (param == NULL)
//      this->param = VPCParamsDefault;
    this->exitReason = ExitReason::UNKNOWN;
    for (int t = 0; t < static_cast<int>(VPCTimeStats::NUM_TIME_STATS); t++) {
      timer.register_name(VPCTimeStatsName[t]);
    }
    this->cutType.resize(0);
    this->numCutsOfType.resize(static_cast<int>(CutType::NUM_CUT_TYPE_STATS), 0);
    this->numFails.resize(static_cast<int>(FailureType::NUM_FAILURES), 0);
    this->branching_lb = std::numeric_limits<double>::max();
    this->branching_ub = std::numeric_limits<double>::lowest();
    this->min_nb_obj_val = std::numeric_limits<double>::max();
    this->ip_opt = std::numeric_limits<double>::max();
    this->num_cgs = 1; // 1 unless assumed otherwise
    this->num_cuts = 0;
    this->num_disj_terms = 0;
    this->num_partial_bb_nodes = 0;
    this->num_pruned_nodes = 0;
    this->num_obj_tried = 0;
  }
} /* setParams */

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
  const double time = timer.get_total_time(VPCTimeStatsName[VPCTimeStats::PRLP_SOLVE_TIME]);
  if (!reached_limit && num_obj_tried >= MANY_OBJ && time > 10.
      && time / num_obj_tried >= params.get(PRLP_TIMELIMIT)) { // checks if average PRLP solve time is too high
    reached_limit = true;
  }
//  if (reached_limit) {
//    this->exitReason = ExitReason::FAIL_LIMIT_EXIT;
//  }
  return reached_limit;
} /* reachedFailureLimit */
