/**
 * @file PRLP.cpp
 * @author A. M. Kazachkov
 * @date 2019-02-14
 */
#include "PRLP.hpp"
#include "CutHelper.hpp"
#include "nbspace.hpp"
#include "Disjunction.hpp"
#include "SolverHelper.hpp"

#include <numeric> // inner_product

#include <CglGMI.hpp>
#include <CoinTime.hpp>

#ifdef TRACE
#include "debug.hpp"
#endif

using namespace VPCParametersNamespace;

/** Default constructor */
PRLP::PRLP() {
  initialize();
} /* default constructor */

/** Owner constructor */
PRLP::PRLP(CglVPC* owner) {
  initialize();
  this->owner = owner;
} /* default constructor */

/** Copy constructor */
PRLP::PRLP(const PRLP& source) {
  initialize(&source);
} /* copy constructor */

/** Destructor */
PRLP::~PRLP() {
} /* destructor */

/** Assignment operator */
PRLP& PRLP::operator=(const PRLP& source) {
  if (this != &source) {
    initialize(&source);
  }
  return *this;
} /* assignment operator */

/** Clone */
OsiSolverInterface* PRLP::clone(bool copyData) const {
  return new PRLP(*this);
} /* clone */

/**
 * Set up solver from points and rays provided.
 * The rows are \alpha (p or r) \ge \beta = 1 (or 0 or -1)
 * @return Is prlp primal feasible?
 */
bool PRLP::setup(const double scale) {
  owner->timer.start_timer(static_cast<int>(CglVPC::VPCTimeStats::PRLP_SETUP_TIME));
  const int num_constraints = owner->prlpData.rhs.size();
  const int num_disj_terms = owner->disj()->num_terms;
  const int num_cols = owner->probData.num_cols;

  // Build the matrix as well as col lower/upper bounds
  std::vector<double> colLB(num_cols, -1 * owner->params.get(doubleConst::INF));
  std::vector<double> colUB(num_cols, owner->params.get(doubleConst::INF));
  std::vector<double> rowLB;
  rowLB.reserve(num_constraints);
  int numElementsEstimate = 0; // estimate max size of sparse coefficient matrix, for reserve

  // We will keep the matrix sorted from least to greatest orthogonality with the objective
  // The points will be first, and then the rays
  std::vector<int> rowOfPoint, rowOfRay;
  std::vector<double> objDepthPoints, objDepthRays;
  rowOfPoint.reserve(num_disj_terms);
  rowOfRay.reserve(num_disj_terms * num_cols);
  objDepthPoints.reserve(num_disj_terms);
  objDepthRays.reserve(num_disj_terms * num_cols);

  // First pass is just to get the orthogonalities and fill the vectors,
  // but we will also remove any singleton rows by setting them as bounds
  for (int r = 0; r < num_constraints; r++) {
    const int num_vec_el = owner->prlpData.constraints[r].getNumElements();
    const double vec_rhs =
        (owner->prlpData.rhs[r] != 0.) ?
            owner->prlpData.rhs[r] * scale : owner->prlpData.rhs[r];
    // Check if it is a bound
    if (num_vec_el == 1) {
      // row is coeff * var >= rhs
      const double coeff = owner->prlpData.constraints[r].getElements()[0];
      const int var = owner->prlpData.constraints[r].getIndices()[0];
      const double old_lb = colLB[var];
      const double old_ub = colUB[var];
      if (isZero(vec_rhs)) {
        if (lessThanVal(coeff, 0.)) {
          if (old_ub > 0.) {
            colUB[var] = 0.;
          }
        } else if (greaterThanVal(coeff, 0.)) {
          if (old_lb < 0.) {
            colLB[var] = 0.;
          }
        } else {
          // Do nothing, because this is a 0 >= 0 row (shouldn't happen... but who knows)
        }
      } /* if zero rhs */
      else {
        if (isVal(coeff, 0.)) {
          // Do nothing, because this is 0 >= something
          if (lessThanVal(vec_rhs, 0.)) {
            error_msg(errorstring,
                "We have an infeasible row: %.2f * x%d >= %.2f.\n", coeff,
                var, vec_rhs);
            writeErrorToLog(errorstring, owner->params.logfile);
            exit(1);
          }
        } else {
          const double new_bound = vec_rhs / coeff;
          if (lessThanVal(coeff, 0.)) {
            if (old_ub > new_bound) {
              colUB[var] = new_bound;
            }
          } else {
            if (old_lb < new_bound) {
              colLB[var] = new_bound;
            }
          }
        }
      } // else it is a non-zero rhs
    } // if num_vec_el == 1, it is a bound

    // When rhs non-zero, we input it anyway, to keep for objective function purposes
    if (num_vec_el > 1 || !isZero(vec_rhs)) {
      numElementsEstimate += num_vec_el;
      if (!isZero(vec_rhs)) {
        rowOfPoint.push_back(r);
        objDepthPoints.push_back(owner->prlpData.objViolation[r]);
      } else {
        rowOfRay.push_back(r);
        objDepthRays.push_back(owner->prlpData.objViolation[r]);
      }
    } // we add it for processing
  } // iterate over the rows in the matrix

  // Sort by orthogonality
  const int numPointsToProcess = objDepthPoints.size();
  std::vector<int> sortIndexPoints(numPointsToProcess);
  for (int i = 0; i < numPointsToProcess; i++) {
    sortIndexPoints[i] = i;
  }
  sort(sortIndexPoints.begin(), sortIndexPoints.end(),
      index_cmp_asc<const double*>(objDepthPoints.data())); // ascending order

  const int numRaysToProcess = objDepthRays.size();
  std::vector<int> sortIndexRays(numRaysToProcess);
  for (int i = 0; i < numRaysToProcess; i++) {
    sortIndexRays[i] = i;
  }
  sort(sortIndexRays.begin(), sortIndexRays.end(),
      index_cmp_asc<const double*>(objDepthRays.data())); // ascending order

  // Should we check for duplicates?
  const bool checkForDuplicates = true;
  std::vector<bool> isDuplicatePoint, isDuplicateRay;
  if (checkForDuplicates) {
    // Check points for duplicates
    isDuplicatePoint.resize(numPointsToProcess, false);
    for (int p1 = 1; p1 < numPointsToProcess; p1++) {
      const int ind1 = sortIndexPoints[p1];
      const int r1 = rowOfPoint[ind1];
      const double ortho1 = objDepthPoints[ind1];
      if (!isVal(objDepthPoints[sortIndexPoints[p1-1]], ortho1)) {
        continue;
      }
      int howDissimilar = 2;
      const CoinShallowPackedVector vec1 = owner->prlpData.constraints[r1];
      const double rhs1 = owner->prlpData.rhs[r1];
      // Iterate through previous points until bottom or the orthogonality is different
      for (int p2 = p1 - 1; p2 >= 0; p2--) {
        if (isDuplicatePoint[p2]) {
          continue;
        }
        const int ind2 = sortIndexPoints[p2];
        const double ortho2 = objDepthPoints[ind2];
        if (!isVal(ortho1, ortho2)) {
          break;
        }
        const int r2 = rowOfPoint[ind2];
        const CoinShallowPackedVector vec2 = owner->prlpData.constraints[r2];
        const double rhs2 = owner->prlpData.rhs[r2];
        // Ensure constants line up
        if (!((greaterThanVal(rhs1, 0.) && greaterThanVal(rhs2, 0.))
            || (lessThanVal(rhs1, 0.) && lessThanVal(rhs2, 0.)))) {
          continue;
        }
        howDissimilar = isRowDifferent(vec1, rhs1, vec2, rhs2, owner->params.get(doubleConst::DIFFEPS));
        if (howDissimilar == -1) {
          // This means the current row is better somehow, and we should replace r2
          isDuplicatePoint[p2] = true;
          numElementsEstimate -= vec2.getNumElements();
        }
        if (howDissimilar != 2) {
          isDuplicatePoint[p1] = true;
          numElementsEstimate -= vec1.getNumElements();
          break;
        }
      } // iterate over previously added points that have same objDepth
    } // iterate over points to find duplicates

    // Check rays for duplicates
    isDuplicateRay.resize(numRaysToProcess, false);
    for (int p1 = 1; p1 < numRaysToProcess; p1++) {
      const int ind1 = sortIndexRays[p1];
      const int r1 = rowOfRay[ind1];
      const CoinShallowPackedVector vec1 = owner->prlpData.constraints[r1];
      const double rhs1 = owner->prlpData.rhs[r1];
      const double ortho1 = objDepthRays[ind1];
      int howDissimilar = 2;
      // Iterate through previous rays until bottom or the orthogonality is different
      for (int p2 = p1 - 1; p2 >= 0; p2--) {
        if (isDuplicateRay[p2]) {
          continue;
        }
        const int ind2 = sortIndexRays[p2];
        const double ortho2 = objDepthRays[ind2];
        if (!isVal(ortho1, ortho2)) {
          break;
        }
        const int r2 = rowOfRay[ind2];
        const CoinShallowPackedVector vec2 = owner->prlpData.constraints[r2];
        const double rhs2 = owner->prlpData.rhs[r2];
        howDissimilar = isRowDifferent(vec1, rhs1, vec2, rhs2, owner->params.get(doubleConst::DIFFEPS));
        if (howDissimilar == -1) {
          // This means the current row is better somehow, and we should replace r2
          isDuplicateRay[p2] = true;
          numElementsEstimate -= vec2.getNumElements();
//          rowToDelete.push_back(disjTermStart[t2] + r2);
        }
        if (howDissimilar != 2) {
          isDuplicateRay[p1] = true;
          numElementsEstimate -= vec1.getNumElements();
          break;
        }
      } // iterate over previously added rays that have same objDepth
    } // iterate over rays to find duplicates
  } // check for duplicates

  // Now input into the matrix
  CoinPackedMatrix mx; // create a new col-ordered matrix
  mx.reverseOrdering(); // make it row-ordered
  mx.setDimensions(0, num_cols);
  mx.reserve(numPointsToProcess + numRaysToProcess, numElementsEstimate, false);

  // First the points
  ortho.clear();
  ortho.reserve(numPointsToProcess + numRaysToProcess);
  for (int p = 0; p < numPointsToProcess; p++) {
    if (checkForDuplicates && isDuplicatePoint[p]) {
      continue;
    }
    const int ind = sortIndexPoints[p];
    const int r = rowOfPoint[ind];
    const double rhs = owner->prlpData.rhs[r] * scale;
    mx.appendRow(owner->prlpData.constraints[r].getNumElements(),
        owner->prlpData.constraints[r].getIndices(),
        owner->prlpData.constraints[r].getElements());
    rowLB.push_back(rhs);
    numPoints++;
    ortho.push_back(objDepthPoints[ind]);
  } // input points

  // Then the rays
  for (int p = 0; p < numRaysToProcess; p++) {
    if (checkForDuplicates && isDuplicateRay[p]) {
      continue;
    }
    const int ind = sortIndexRays[p];
    const int r = rowOfRay[ind];
    const double rhs = owner->prlpData.rhs[r];
    mx.appendRow(owner->prlpData.constraints[r].getNumElements(),
        owner->prlpData.constraints[r].getIndices(),
        owner->prlpData.constraints[r].getElements());
    rowLB.push_back(rhs);
    numRays++;
    ortho.push_back(objDepthRays[ind]);
  } // input rays

  // Clean matrix
  mx.cleanMatrix();

  // Filter out columns that are empty
  const bool shouldDelete = true;
  if (shouldDelete) {
    std::vector<int> deleteColIndex;
    const int* length = mx.countOrthoLength();
    for (int col_ind = 0; col_ind < mx.getNumCols(); col_ind++) {
      if (length[col_ind] == 0) {
        const bool ignore_lb = isNegInfinity(colLB[col_ind]) || isZero(colLB[col_ind]);
        const bool ignore_ub = isInfinity(colUB[col_ind]) || isZero(colUB[col_ind]);
        if (!ignore_lb || !ignore_ub) {
          nonZeroColIndex.push_back(col_ind);
        } else {
          deleteColIndex.push_back(col_ind);
        }
      } else {
        nonZeroColIndex.push_back(col_ind);
      }
    }
    if (length) {
      delete[] length;
      length = NULL;
    }
    // Delete the columns from mx, as well as from colLB and colUB
    mx.deleteCols(deleteColIndex.size(), deleteColIndex.data());
    for (int pos = deleteColIndex.size() - 1; pos >= 0; pos--) {
      colLB.erase(colLB.begin() + deleteColIndex[pos]);
      colUB.erase(colUB.begin() + deleteColIndex[pos]);
    }
  }

  // Load problem into prlp
  // Defaults (in COIN-OR):
  // colLB: 0 **** This should be -inf (except when changed above)
  // colUB: inf
  // obj coeff: 0
  // rowSense: >=
  // rhs: 0 **** We want to change this to beta for the points
  // rowRange: 0
  this->loadProblem(mx, colLB.data(), colUB.data(), NULL, NULL,
      rowLB.data(), NULL);
  this->disableFactorization();
  //this->getModelPtr()->cleanMatrix();

#ifdef TRACE
  printf(
      "\n## Check feasibility using all zeroes objective. ##\n");
#endif

  owner->timer.end_timer(static_cast<int>(CglVPC::VPCTimeStats::PRLP_SETUP_TIME));

  if (numPoints < 0 || numRays < 0) {
    return true;
  }
  this->density =
    (double) this->getMatrixByCol()->getNumElements()
    / (this->getNumRows() * this->getNumCols());

  const double max_time = CoinMax(60., 10 * owner->params.get(PRLP_TIMELIMIT));
  setLPSolverParameters(this, owner->params.get(VERBOSITY), max_time);
#ifdef USE_CLP
  this->getModelPtr()->setMaximumIterations(static_cast<int>(std::numeric_limits<int>::max()));
  this->getModelPtr()->setNumberIterations(0);
#endif
  this->setHintParam(OsiDoPresolveInInitial, owner->params.get(intConst::PRLP_PRESOLVE) >= 1);
  this->setHintParam(OsiDoPresolveInResolve, owner->params.get(intConst::PRLP_PRESOLVE) >= 2);

  // Check that the prlp is feasible for the zero objective function
  owner->timer.start_timer(CglVPC::ObjectiveTypeName[static_cast<int>(CglVPC::ObjectiveType::DUMMY_OBJ)] + "_TIME");

  const double start = CoinCpuTime();
  this->initialSolve();
  const double end = CoinCpuTime();

  // Set PRLP time limit based on this initial solve, if the time is not yet set and it is not permanently set to infinity
  if (isVal(owner->params.get(PRLP_TIMELIMIT), -1.)) {
    owner->params.set(PRLP_TIMELIMIT,
        CoinMax(owner->params.get(doubleConst::MIN_PRLP_TIMELIMIT), (end - start)));
    setTimeLimit(this, owner->params.get(PRLP_TIMELIMIT));
  }

  owner->timer.end_timer(CglVPC::ObjectiveTypeName[static_cast<int>(CglVPC::ObjectiveType::DUMMY_OBJ)] + "_TIME");

  // Check that prlp is not primal infeasible
  // This can happen if 0 \in cone(\rayset), for instance, in the nonbasic space.
  // We might also get that the cut LP is "bad" and does not solve quickly
  // in which case we abort with numerical issues as the error.
  bool retval;
  if (this->isProvenOptimal()) {
    retval = true;
  } else if (this->isProvenPrimalInfeasible() && !this->isAbandoned()) {
    warning_msg(errorstring,
        "PRLP is primal infeasible; giving up on this point-ray collection.\n");
    owner->numFails[static_cast<int>(CglVPC::FailureType::PRIMAL_INFEASIBLE_NO_OBJ)]++;
    retval = false;
  } else if (this->isIterationLimitReached()
      || this->getModelPtr()->hitMaximumIterations()) {
    warning_msg(errorstring,
        "PRLP is having numerical issues; giving up on this point-ray collection.\n");
    owner->numFails[static_cast<int>(CglVPC::FailureType::NUMERICAL_ISSUES_NO_OBJ)]++;
    retval = false;
  } else if (this->isProvenDualInfeasible() && !this->isAbandoned()) {
    // This should never happen, but we are going to chock it up to numerical problems
    warning_msg(errorstring,
        "PRLP is having numerical issues; giving up on this point-ray collection.\n");
    owner->numFails[static_cast<int>(CglVPC::FailureType::NUMERICAL_ISSUES_NO_OBJ)]++;
    retval = false;
  } else {
    // Something else happened; still bad...
    warning_msg(errorstring,
        "PRLP is failing, probably due to numerical issues; giving up on this point-ray collection.\n");
    owner->numFails[static_cast<int>(CglVPC::FailureType::UNKNOWN)]++;
    retval = false;
  }
  return retval;
} /* setup */

/****************** PROTECTED **********************/

void PRLP::initialize(const PRLP* const source) {
  if (source != NULL) {
    this->owner = source->owner;
    this->nonZeroColIndex = source->nonZeroColIndex;
    this->ortho = source->ortho;
    this->numPoints = source->numPoints;
    this->numRays = source->numRays;
    this->density = source->density;
    this->num_obj_tried = source->num_obj_tried;
    this->num_cuts = source->num_cuts;
    this->num_failures = source->num_failures;
  }
  else {
    this->owner = NULL;
    this->nonZeroColIndex.resize(0);
    this->ortho.resize(0);
    this->numPoints = 0;
    this->numRays = 0;
    this->density = 0.;
    this->num_obj_tried = 0;
    this->num_cuts = 0;
    this->num_failures = 0;
  }
} /* initialize */

void PRLP::setObjectiveFromStructuralPoint(const double* const pointVals,
    const double* const pointSlack, const OsiSolverInterface* const origSolver,
    const bool inNBSpace) {
  CoinPackedVector vec;
  if (inNBSpace) {
    double obj_viol;
    setCompNBCoor(vec, obj_viol, owner->params, pointVals, pointSlack,
        origSolver, owner->probData.NBVarIndex, owner->probData.NBReducedCost);
    addToObjectiveFromPackedVector(this, &vec, true, 1.,
        (nonZeroColIndex.size() > 0) ? &nonZeroColIndex : NULL);
  } else {
    if (nonZeroColIndex.size() > 0) {
      for (int i : nonZeroColIndex) {
        this->setObjCoeff(i, pointVals[nonZeroColIndex[i]]);
      }
    } else {
      this->setObjective(pointVals);
    }
  }
} /* setObjectiveFromStructuralPoint */

int PRLP::tryOneObjective(std::vector<int>& numTimesTightRow,
    std::vector<int>& numTimesTightColLB, std::vector<int>& numTimesTightColUB,
    OsiCuts& cuts, const OsiSolverInterface* const origSolver,
    const double beta, const OsiCuts* const structSICs, const bool inNBSpace,
    const CglVPC::ObjectiveType cutHeur, const bool tryExtraHard) {
  int return_code = genCut(cuts, origSolver, beta, structSICs,
      inNBSpace, cutHeur, tryExtraHard);
  if (this->isProvenOptimal()) {
    updateStepForTargetedCutGeneration(numTimesTightRow, numTimesTightColLB,
        numTimesTightColUB);
  }
  return return_code;
} /* tryOneObjective */

int PRLP::resolvePRLP(const bool tryExtraHard) {
  this->getModelPtr()->setNumberIterations(0);
  const double timeLimit = (tryExtraHard) ?  -1. : owner->params.get(PRLP_TIMELIMIT); // maybe go unlimited time this round
  setTimeLimit(this, timeLimit);
  this->resolve();

//  // If the objective value seems weird, we may want to try an initial solve
//  if (this->isProvenOptimal()) {
//    const double objValue = this->getObjValue();
//    if (isZero(objValue, owner->params.get(EPS))
//        && !isZero(objValue, owner->params.get(doubleConst::DIFFEPS))) {
//      this->getModelPtr()->setNumberIterations(0);
//      setTimeLimit(this, timeLimit);
//      this->initialSolve();
//    }
//  }

  checkSolverOptimality(this, false, timeLimit);

  // If "nan" objective, then something went terribly wrong; try again
  // This can happen when repeatedly solving an OsiClpSolverInterface LP
  if (!this->isProvenPrimalInfeasible() && std::isnan(this->getObjValue())) {
#ifdef TRACE
    warning_msg(warnstr,
        "Encountered NaN objective value.\n");
#endif
    this->getModelPtr()->setNumberIterations(0);
    setTimeLimit(this, timeLimit);
    this->initialSolve();
    checkSolverOptimality(this, false, timeLimit);
  }

  if (this->isProvenOptimal()
      && isNegInfinity(this->getObjValue(), -owner->params.get(doubleConst::INF))) {
    // We essentially have a ray, though it is not being reported as such
    // Try resolving because this fixes the issue sometimes
    setTimeLimit(this, timeLimit);
    this->initialSolve();
    checkSolverOptimality(this, false, timeLimit);
  }

  if (this->isProvenOptimal()) {
    return 0;
  }
  else if (this->isProvenDualInfeasible()) {
    return -1 * (static_cast<int>(CglVPC::FailureType::DUAL_INFEASIBLE) + 1);
  }
  else if (this->isIterationLimitReached()) {
    return -1 * (static_cast<int>(CglVPC::FailureType::ITERATION_LIMIT) + 1);
  }
  else if (this->getModelPtr()->hitMaximumIterations()) {
    return -1 * (static_cast<int>(CglVPC::FailureType::TIME_LIMIT) + 1);
  }
  else if (this->isAbandoned()) {
    return -1 * (static_cast<int>(CglVPC::FailureType::ABANDONED) + 1);
  }
  return -1 * (static_cast<int>(CglVPC::FailureType::UNKNOWN) + 1);
} /* resolvePRLP */

/**
 * @brief Generate a cut, if possible, from the solver
 * Rhs is given as beta
 *
 * @return -1 * (fail index + 1) if it is something that means to discontinue this solver
 */
int PRLP::genCut(OsiCuts& cuts, const OsiSolverInterface* const origSolver,
    const double beta, const OsiCuts* const structSICs, const bool inNBSpace,
    const CglVPC::ObjectiveType cutHeur, const bool tryExtraHard) {
  owner->timer.start_timer(CglVPC::VPCTimeStatsName[static_cast<int>(CglVPC::VPCTimeStats::PRLP_SOLVE_TIME)]);
  this->num_obj_tried++;
  owner->num_obj_tried++;
  owner->numObjFromHeur[static_cast<int>(cutHeur)]++;
  if (owner->reachedCutLimit()) {
    return exitGenCut(-1 * (static_cast<int>(CglVPC::FailureType::CUT_LIMIT) + 1), cutHeur);
  }

  // First we save the old solution; this is just to have a quick way to check that the same solution was found
  std::vector<double> old_solution;
  if (this->num_obj_tried > 0 && this->isProvenOptimal()) {
    old_solution.assign(this->getColSolution(), this->getColSolution() + this->getNumCols());
  }

  // Resolve
  int return_code = resolvePRLP(tryExtraHard);
  if (return_code < 0 && return_code != -1 * (static_cast<int>(CglVPC::FailureType::DUAL_INFEASIBLE) + 1)) {
    // Homogeneous cuts never cut the LP optimum in NB space, but we leave that decision to later
    return exitGenCut(return_code, cutHeur);
  }

  // If it is proven optimal, we can simply obtain the solution
  // If it is dual infeasible, then it is unbounded (or infeasible, but it would say so)
  // in which case we need to get the ray and use that as the solution
  if (this->isProvenOptimal()) {
    // Try to generate the cut (but first, test against previous solution)
    if (!old_solution.empty()) {
      double maxDiff = 0.;
      for (int ind = 0; ind < this->getNumCols(); ind++) {
        const double diff = std::abs(this->getColSolution()[ind] - old_solution[ind]);
        if (diff > maxDiff) {
          maxDiff = diff;
        }
      }
      if (isZero(maxDiff, owner->probData.EPS)) {
        return exitGenCut(-1 * (static_cast<int>(CglVPC::FailureType::DUPLICATE_VPC) + 1), cutHeur);
      }
    } // check whether it is the same solution as before

    return_code = genCutHelper(cuts, origSolver, beta, structSICs,
        inNBSpace, cutHeur);
  } // check if proven optimal
  else if (this->isProvenDualInfeasible()) {
    // For the time being this seems to not be working well
    // For bell3a for instance, when cutting 64 rays, there is a cut that appears
    // both as a primal ray and as a feasible solution, so it has both >= 0 and >= -1
    // constant sides, and it causes infeasibility of the original problem when added
    return_code = -1 * (static_cast<int>(CglVPC::FailureType::DUAL_INFEASIBLE) + 1);
  } // generate rays? (currently non-functional)

  if (!this->isProvenOptimal() && !this->isProvenDualInfeasible()) {
    warning_msg(errorstring,
        "PRLP is not proven optimal or dual infeasible. Primal infeasible? %d.\n",
        this->isProvenPrimalInfeasible());
  }
  return exitGenCut(return_code, cutHeur);
} /* genCut */

int PRLP::genCutHelper(OsiCuts& cuts,
    const OsiSolverInterface* const origSolver, const double beta,
    const OsiCuts* const structSICs, const bool inNBSpace,
    const CglVPC::ObjectiveType cutHeur) {
  // We might get here without this being optimal
  if (!this->isProvenOptimal()) {
    return 0; // errors are tabulated in exitGenCutsFromCutSolver
  } // may have run into issues
  if (owner->reachedCutLimit()) {
    return -1 * (static_cast<int>(CglVPC::FailureType::CUT_LIMIT) + 1);
  }

  int num_cuts_generated = 0;
  OsiRowCut currCut;
  if (inNBSpace) {
    setCutFromNBCoefficients(&currCut, this->nonZeroColIndex,
        this->getColSolution(), beta, origSolver,
        owner->probData.NBVarIndex, true, owner->probData.EPS);
  } else {
    setOsiRowCut(&currCut, nonZeroColIndex, this->getNumCols(),
        this->getColSolution(), beta, owner->probData.EPS);
  }

  // Clean
  const int clean_error = cleanCut(&currCut, origSolver, owner->params,
      owner->probData.minAbsCoeff, owner->probData.maxAbsCoeff,
      owner->probData.EPS, beta,
      !inNBSpace);
  if (clean_error != 0) {
    return clean_error;
  }

  // Ensure the cut is not a duplicate, and some other checks
  const bool generateAnywayIfDuplicateSIC = true; // in principle, when this is true, we do not need to check ortho with SICs; in our tests, we might anyway if we want to get the relevant stats
  const bool replaceOldWorseCut = owner->canReplaceGivenCuts;
  const bool checkIfDuplicateSIC = false && (structSICs != NULL);
  const double currCutNorm = currCut.row().twoNorm();
  const double violation = currCut.violated(origSolver->getColSolution()); // beta normalization fixes the violation in nb-space
  int duplicateSICIndex = -1, duplicateGICIndex = -1;
  int minOrthogonalityIndex = -1;
  double orthogonalityWithSICs = 1., orthogonalityWithGICs = 1.;
  if (checkIfDuplicateSIC) {
    howDuplicate(*structSICs, currCut, 0, duplicateSICIndex, minOrthogonalityIndex,
        orthogonalityWithSICs, owner->params.get(MIN_ORTHOGONALITY),
        owner->probData.EPS);
  }
  const bool duplicateSICFlag = (duplicateSICIndex >= 0);
  const bool orthogonalitySICFailFlag = !duplicateSICFlag
      && (orthogonalityWithSICs < owner->params.get(MIN_ORTHOGONALITY));
  const int startIndex = replaceOldWorseCut ? 0 : owner->init_num_cuts;
  if (generateAnywayIfDuplicateSIC || (!duplicateSICFlag && !orthogonalitySICFailFlag)) {
    // Check again all of the cuts because we don't want a cut duplicate to an old one,
    // even if we are not replacing the old ones
    howDuplicate(cuts, currCut, 0, duplicateGICIndex,
        minOrthogonalityIndex, orthogonalityWithGICs,
        owner->params.get(MIN_ORTHOGONALITY), owner->probData.EPS);
  }
  const bool duplicateGICFlag = (duplicateGICIndex >= 0);
  bool orthogonalityGICFailFlag = !duplicateGICFlag
      && (orthogonalityWithGICs < owner->params.get(MIN_ORTHOGONALITY));
  if (orthogonalityGICFailFlag && (minOrthogonalityIndex >= startIndex)) {
    // If rhs is better or sparser than existing cut, replace the cut
    OsiRowCut* oldCut = cuts.rowCutPtr(minOrthogonalityIndex);
    const int cut_num_el = oldCut->row().getNumElements();
    const bool isSparser = currCut.row().getNumElements() < cut_num_el;
    const double old_effectiveness = oldCut->effectiveness();
    const bool isMoreEffective = (violation / currCutNorm) > old_effectiveness;
    if (isSparser || isMoreEffective) {
      // Adjust statistics
      owner->numCutsFromHeur[static_cast<int>(owner->objType[minOrthogonalityIndex - startIndex])]--;
      owner->numCutsFromHeur[static_cast<int>(cutHeur)]++;
      owner->numCutsOfType[static_cast<int>(owner->cutType[minOrthogonalityIndex - startIndex])]--;
      owner->numCutsOfType[static_cast<int>(CglVPC::CutType::VPC)]++;

      // Replace old cut
      owner->objType[minOrthogonalityIndex - startIndex] = cutHeur;
      owner->cutType[minOrthogonalityIndex - startIndex] = CglVPC::CutType::VPC;
      currCut.setEffectiveness(violation / currCutNorm);
      *oldCut = currCut; // TODO double check this works

      orthogonalityGICFailFlag = true;
    }
  } // orthogonalityGICFail (replace cut possibly)
  else {
    const bool toAdd = (generateAnywayIfDuplicateSIC || (!duplicateSICFlag && !orthogonalitySICFailFlag))
        && !duplicateGICFlag && !orthogonalityGICFailFlag;
    if (toAdd) {
      currCut.setEffectiveness(violation / currCutNorm);
      owner->addCut(currCut, cuts, CglVPC::CutType::VPC, cutHeur);
      this->num_cuts++;
      num_cuts_generated++;
      owner->numFails[static_cast<int>(CglVPC::FailureType::DUPLICATE_SIC)] += duplicateSICFlag;
      owner->numFails[static_cast<int>(CglVPC::FailureType::ORTHOGONALITY_SIC)] += orthogonalitySICFailFlag;
      return num_cuts_generated;
    } else {
      // Note that these are all mutually exclusive IF we are checking for duplicates
      if (generateAnywayIfDuplicateSIC) {
        owner->numFails[static_cast<int>(CglVPC::FailureType::DUPLICATE_SIC)] += duplicateSICFlag;
        owner->numFails[static_cast<int>(CglVPC::FailureType::ORTHOGONALITY_SIC)] += orthogonalitySICFailFlag;
      } else if (duplicateSICFlag) {
        return -1 * (static_cast<int>(CglVPC::FailureType::DUPLICATE_SIC) + 1);
      } else if (orthogonalitySICFailFlag) {
        return -1 * (static_cast<int>(CglVPC::FailureType::ORTHOGONALITY_SIC) + 1);
      }
      if (duplicateGICFlag) {
        return -1 * (static_cast<int>(CglVPC::FailureType::DUPLICATE_VPC) + 1);
      }
      if (orthogonalityGICFailFlag) {
        return -1 * (static_cast<int>(CglVPC::FailureType::ORTHOGONALITY_VPC) + 1);
      }
    }
  }
  return num_cuts_generated;
} /* genCutHelper */

/**
 * Returns error status or num_cuts_generated
 */
int PRLP::exitGenCut(const int num_cuts_generated, const CglVPC::ObjectiveType cutHeur) {
  owner->timer.end_timer(CglVPC::VPCTimeStatsName[static_cast<int>(CglVPC::VPCTimeStats::PRLP_SOLVE_TIME)]);
  int return_status = num_cuts_generated;
  if (return_status >= 0 && !this->isProvenOptimal()) {
    // Sometimes it simply takes too long to solve optimally
    if (this->isIterationLimitReached()) {
      return_status = -1 * (static_cast<int>(CglVPC::FailureType::ITERATION_LIMIT)+ 1);
    }
    // The next will be true if either # iterations or time limit reached
    else if (hitTimeLimit(this)) {
      return_status = -1 * (static_cast<int>(CglVPC::FailureType::TIME_LIMIT)+ 1);
    }
    // Numerical difficulties may have been encountered
    else if (this->isAbandoned()) {
      return_status = -1 * (static_cast<int>(CglVPC::FailureType::ABANDONED)+ 1);
    }
    // The problem may be unbounded
    else if (this->isProvenDualInfeasible()) {
      return_status = -1 * (static_cast<int>(CglVPC::FailureType::DUAL_INFEASIBLE)+ 1);
    }
    // We may get an infeasible problem
    else if (this->isProvenPrimalInfeasible()) {
      return_status = -1 * (static_cast<int>(CglVPC::FailureType::PRIMAL_INFEASIBLE)+ 1);
    }
    // Something super strange may have happened
    else {
      return_status = -1 * (static_cast<int>(CglVPC::FailureType::UNKNOWN)+ 1);
    }

    // Make sure we did not generate cuts
    if (num_cuts_generated  > 0) {
      error_msg(errorstring,
          "Somehow we generated %d cuts from a non-optimal PRLP; not possible. Error: %d.\n",
          num_cuts_generated, -1 * (return_status + 1));
      writeErrorToLog(errorstring, owner->params.logfile);
      exit(1);
    }
  } // may have run into issues

  if (return_status < 0) {
    const int error_index = -1 * (return_status + 1);
    owner->numFails[error_index]++;
  }

  if (return_status <= 0){
    this->num_failures++;
    owner->num_failures++;
    owner->numFailsFromHeur[static_cast<int>(cutHeur)]++;
  }

  // Reset timer and number iterations
  this->getModelPtr()->setNumberIterations(0);
  setTimeLimit(this, owner->params.get(PRLP_TIMELIMIT));
  return return_status;
} /* exitGenCutsFromCutSolver */

/**
 * Collect all the points and rays and sort them by the objective
 */
void PRLP::setupForTargetedCutGeneration(std::vector<rowAndActivity>& pointIndex,
    std::vector<rowAndActivity>& rayIndex) {
  pointIndex.clear();
  rayIndex.clear();
  for (int row_ind = 0; row_ind < this->getNumRows(); row_ind++) {
    rowAndActivity tmp;
    tmp.row = row_ind;
    tmp.activity = ortho[row_ind];
    if (isZero(this->getRightHandSide()[row_ind])) {
      rayIndex.push_back(tmp);
    } else {
      pointIndex.push_back(tmp);
    }
  }
} /* setupForTargetedCutGeneration */

/**
 * After a cut has been generated, mark all the tight rows
// * For the rays, combine them into one and add it to compareRays
 * In addition, update the set of rays to be considered
 */
int PRLP::updateStepForTargetedCutGeneration(std::vector<int>& numTimesTightRow,
    std::vector<int>& numTimesTightColLB, std::vector<int>& numTimesTightColUB) const {
  if (!this->isProvenOptimal()) {
    return 0;
  }

  int num_changed = 0;
  // Check each row
  for (int row_ind = 0; row_ind < this->getNumRows(); row_ind++) {
    const double activity = this->getRowActivity()[row_ind];
    const double rhs = this->getRightHandSide()[row_ind];
    if (isVal(activity, rhs)) {
      if (numTimesTightRow[row_ind] == 0) {
        num_changed++;
      }
      if (numTimesTightRow[row_ind] >= 0) {
        numTimesTightRow[row_ind]++;
      } else {
        numTimesTightRow[row_ind]--;
      }
    }
  } // check rows

  // Also process the axis directions
  for (int col_ind = 0; col_ind < this->getNumCols(); col_ind++) {
    const double val = this->getColSolution()[col_ind];
    const double lb = this->getColLower()[col_ind];
    const double ub = this->getColUpper()[col_ind];
    if (!isNegInfinity(lb) && isVal(val, lb)) {
      if (numTimesTightColLB[col_ind] == 0) {
        num_changed++;
      }
      if (numTimesTightColLB[col_ind] >= 0) {
        numTimesTightColLB[col_ind]++;
      } else {
        numTimesTightColLB[col_ind]--;
      }
    }
    if (!isInfinity(ub) && isVal(val, ub)) {
      if (numTimesTightColUB[col_ind] == 0) {
        num_changed++;
      }
      if (numTimesTightColUB[col_ind] >= 0) {
        numTimesTightColUB[col_ind]++;
      } else {
        numTimesTightColUB[col_ind]--;
      }
    }
  } // check columns

  return num_changed;
} /* updateStepForTargetedCutGeneration */

/**
 * After a new cut, we need to update the sorted list of rays that we are going to use
 * We need to remove all the rays that are now tight and we need to update the orthogonalities
 * with respect to the new ray (we calculate with respect to the cut norm actually)
 */
void PRLP::updateRaySetForTargetedCutGeneration(std::set<compareRayInfo>& sortedRays,
    const std::vector<rowAndActivity>& rayIndex, const int init_num_cuts, int& num_old_cuts,
    const OsiCuts& cuts) {
  const int num_cuts = cuts.sizeCuts();

  std::set<compareRayInfo> newSortedRays;
  for (compareRayInfo curr_ray_info : sortedRays) {
    const int ray_ind = curr_ray_info.index;
    const int row_ind = rayIndex[ray_ind].row;

    // If this ray has a negative row, we skip it
    if (row_ind < 0) {
      continue;  // not a candidate; already has been tight
    }

    double this_total_ortho = 0.;
    const CoinShallowPackedVector vec = this->getMatrixByRow()->getVector(row_ind);
    for (int cut_ind = num_old_cuts; cut_ind < num_cuts; cut_ind++) {
      const CoinPackedVector newOrthogonalRay = cuts.rowCutPtr(cut_ind)->row();

      // The norm to the cut is orthogonal to the rays we would actually like to be orthogonal to
      // In other words, we would like to be parallel to the cut norm
      const double this_ortho = std::abs(getParallelism(vec, newOrthogonalRay));

      if (this_ortho < curr_ray_info.min_ortho) {
        curr_ray_info.min_ortho = this_ortho;
      }
      this_total_ortho += this_ortho;
    }
    // We redo the average
    // (We use num_old_cuts + 1 not num_old_cuts because we initialized with ortho wrt objective)
    curr_ray_info.avg_ortho =
        (curr_ray_info.avg_ortho
        * (num_old_cuts - init_num_cuts + 1) + this_total_ortho)
        / (num_cuts - init_num_cuts + 1);
    newSortedRays.insert(curr_ray_info);
  } // update the sorted rays set
  sortedRays.swap(newSortedRays);
  num_old_cuts = num_cuts;
} /* updateRaySetForTargetedCutGeneration */

/**
 * Try getting cuts tight on the branching lb
 * Returns 0 if error, 1 if no error
 */
int PRLP::findCutsTightOnPoint(std::vector<int>& numTimesTightRow,
    std::vector<int>& numTimesTightColLB, std::vector<int>& numTimesTightColUB,
    const int point_row_ind, const CglVPC::ObjectiveType& cutHeur,
    const double beta, OsiCuts& cuts,
    const OsiSolverInterface* const origSolver, const OsiCuts* const structSICs,
    const std::string& timeName, const bool inNBSpace,
    const int MAX_NUM_OBJ_PER_POINT) {
  const int num_rows = this->getNumRows();
  const int num_cols = this->getNumCols();
  int num_cuts_per_point = 0;
  int check_ind = -1; //strong_lb_row_ind; // row that we last _subtracted_ from the objective
  const int BAD_RETURN_CODE = -1 * (static_cast<int>(CglVPC::FailureType::PRIMAL_INFEASIBLE) + 1);
  int ret_val = 1;

#ifdef TRACE_CUT_BOUND
  // Add a debug solver in order to track the bound after adding the cuts we have generated
  OsiSolverInterface* tmpSolver = origSolver->clone();
  OsiCuts tmpCuts = cuts;
  applyCutsCustom(tmpSolver, tmpCuts);
#endif

  const CoinPackedMatrix* mat = this->getMatrixByRow();
  const CoinShallowPackedVector orig_point = mat->getVector(point_row_ind);
  CoinPackedVector* vec = NULL; // ray that will be added/subtracted
  addToObjectiveFromPackedVector(this, &orig_point, false);

  // Try orig_point as an objective
  int tmp_return_code = tryOneObjective(numTimesTightRow, numTimesTightColLB,
      numTimesTightColUB, cuts, origSolver, beta, structSICs,
      inNBSpace, cutHeur, true);
  if (tmp_return_code > 0) {
    num_cuts_per_point += tmp_return_code;
#ifdef TRACE_CUT_BOUND
    for (int add_ind = 0; add_ind < tmp_return_code; add_ind++) {
      tmpCuts.insert(tmpCuts.rowCut(tmpCuts.sizeCuts() + add_ind));
    }
    applyCutsCustom(tmpSolver, tmpCuts, tmpCuts.sizeCuts() - tmp_return_code,
        tmp_return_code, tmpSolver->getObjValue());
#endif
  }
  setConstantObjectiveFromPackedVector(this, 0., orig_point.getNumElements(), orig_point.getIndices());

  if (tmp_return_code == BAD_RETURN_CODE
      || owner->reachedCutLimit()
      || owner->reachedTimeLimit(timeName, owner->params.get(TIMELIMIT))
      || owner->reachedFailureLimit(num_cuts, num_failures)) {
#ifdef TRACE_CUT_BOUND
    if (tmpSolver) {
      delete tmpSolver;
    }
#endif
    return tmp_return_code != BAD_RETURN_CODE;
  }

  // There better exist a cut that is tight on the point, or we give up
  if (!this->isProvenOptimal()
      || !isVal(this->getRowActivity()[point_row_ind], this->getRightHandSide()[point_row_ind])
      || MAX_NUM_OBJ_PER_POINT <= 1) {
#ifdef TRACE_CUT_BOUND
    if (tmpSolver) {
      delete tmpSolver;
    }
#endif
    return 1;
  }

  // Setup; there are different options for how to choose objectives here
  const int mode_param = owner->params.get(intConst::MODE_OBJ_PER_POINT);
  // 0 = all of the rows that are not tight, and subtract as they get tight
  // 1 = one point/ray at time
  // 2 = keep trying even if the point/ray has become tight in the process?
  // 0 is more expensive at each step and not clear that it would be better
  const int mode_ones = mode_param % 10;
  // 0x: only rays
  // 1x: points+rays
  // 2x: points+rays+variables
  const int mode_param_tens = (mode_param % 100 - (mode_param % 10)) / 10;
  // 0xx: small to large angle with obj (ascending)
  // 1xx: large to small angle with obj (descending)
  // 2xx: small to large slack (ascending)
  // 3xx: large to small slack (descending)
  const int mode_param_hundreds = (mode_param % 1000 - (mode_param % 100)) / 100;
  const bool usePoints = mode_param_tens > 0; // use points as objectives
  const bool useVariables = mode_param_tens > 1; // use variables as objectives
  const int shouldSort = (mode_ones > 0) * mode_param_hundreds; // sort by row activity (no need for mode 0)

  // Trying to find (more) cuts tight on the given point, so set the row to equality temporarily
  const double orig_lb = this->getRowLower()[point_row_ind];
  const double orig_ub = this->getRowUpper()[point_row_ind];
  this->setRowUpper(point_row_ind, orig_lb);

  // Set up vectors to be used in setting the objective
  std::vector<int> indexToCheck; // index order is first rows, then col lb, then col ub
  std::vector<double> sortCriterion;
  std::vector<int> savedNumTimesTightRow(numTimesTightRow);
  std::vector<int> savedNumTimesTightColLB(numTimesTightColLB);
  std::vector<int> savedNumTimesTightColUB(numTimesTightColUB);

  this->enableFactorization(); // note that this changes things slightly sometimes
  std::vector<int> varBasicInRow(num_rows);
  this->getBasics(varBasicInRow.data());

  // Set up slack and rhs vectors (how to sort which points/rays to try)
  for (int row_ind = 0; row_ind < num_rows; row_ind++) {
    const int var = varBasicInRow[row_ind];
    double curr_activity = 0., curr_rhs = 0.;
    if (var < num_cols) {
      if (!useVariables) {
        continue; // skip columns if desired
      }
      const double lb = this->getColLower()[var];
      const double ub = this->getColUpper()[var];
      curr_activity = this->getColSolution()[var];
      if (!isNegInfinity(lb) && numTimesTightColLB[var] >= 0) {
        curr_rhs = lb;
        if (!isZero(curr_activity - curr_rhs)) {
          indexToCheck.push_back(var + num_rows);
          if (shouldSort == 0 || shouldSort == 1) {
            sortCriterion.push_back(
                inNBSpace ?
                    owner->probData.NBReducedCost[var] :
                    std::abs(origSolver->getReducedCost()[var]));
          } else if (shouldSort == 2 || shouldSort == 3) {
            sortCriterion.push_back(curr_activity - curr_rhs);
          }
        }
      }
      if (!isInfinity(ub) && numTimesTightColUB[var] >= 0) {
        curr_rhs = ub;
        if (!isZero(curr_activity - curr_rhs)) {
          indexToCheck.push_back(var + num_rows + num_cols);
          if (shouldSort == 0 || shouldSort == 1) {
            sortCriterion.push_back(
                inNBSpace ?
                    owner->probData.NBReducedCost[var] :
                    std::abs(origSolver->getReducedCost()[var]));
          } else if (shouldSort == 2 || shouldSort == 3) {
            sortCriterion.push_back(curr_rhs - curr_activity); // stay positive
          }
        }
      }
    } else {
      const int row = var - num_cols; // should actually be same as row_ind, but just in case...
      curr_rhs = this->getRightHandSide()[row];
      if (!usePoints && !isZero(curr_rhs) && numTimesTightRow[row] >= 0) {
        continue; // skip points if desired
      }
      curr_activity = this->getRowActivity()[row];
      if (!isZero(curr_activity - curr_rhs)) {
        indexToCheck.push_back(row);
        if (shouldSort == 0 || shouldSort == 1) {
          if (isZero(curr_rhs)) {
            sortCriterion.push_back(ortho[row]); // already normalized
          } else {
            // ortho[row] has c * p2 (current point)
            // We wish to add c * (p2 - p1) / norm(p2 - p1)
            if (vec) {
              delete vec;
            }
            vec = new CoinPackedVector(mat->getVectorSize(row),
                mat->getIndices() + mat->getVectorFirst(row),
                mat->getElements() + mat->getVectorFirst(row), false);
            CoinPackedVector tmpVec;
            packedSortedVectorSum(tmpVec, 1., (*vec), -1., orig_point, 1e-20);
            sortCriterion.push_back(
                (ortho[row] - ortho[point_row_ind]) / tmpVec.twoNorm());
          }
        } else if (shouldSort == 2 || shouldSort == 3) {
          sortCriterion.push_back((curr_activity - curr_rhs) / getRowTwoNorm(row, mat));
        }
      }
    }
  } // set up slack and rhs vectors
  this->disableFactorization();

  // The real # to check
  const int numToCheck = indexToCheck.size();

  // Sort by descending/ascending activity, if requested
  std::vector<int> sortIndex(numToCheck);
  for (int i = 0; i < numToCheck; i++) {
    sortIndex[i] = i;
  }
  if (shouldSort == 0 || shouldSort == 2) {
    std::sort(sortIndex.begin(), sortIndex.end(),
        index_cmp_asc<const std::vector<double>&>(sortCriterion)); // ascending
  } else if (shouldSort == 1 || shouldSort == 3) {
    std::sort(sortIndex.begin(), sortIndex.end(),
        index_cmp_dsc<const std::vector<double>&>(sortCriterion)); // descending
  }

  // Prepare first new objective
  std::vector<int> indexAddedToObjective;
  while (check_ind < numToCheck - 1) {
    check_ind++;
    const int curr_index = indexToCheck[sortIndex[check_ind]];
    if (vec) {
      delete vec;
      vec = NULL;
    }
    if (curr_index < num_rows) {
      vec = new CoinPackedVector(mat->getVectorSize(curr_index),
          mat->getIndices() + mat->getVectorFirst(curr_index),
          mat->getElements() + mat->getVectorFirst(curr_index), false);
    } else if (curr_index < num_rows + num_cols) {
      const int col = curr_index - num_rows;
      const int index[1] = { col };
      const double value[1] = { 1. };
      vec = new CoinPackedVector(1, index, value);
    } else {
      const int col = curr_index - num_rows - num_cols;
      const int index[1] = { col };
      const double value[1] = { -1. };
      vec = new CoinPackedVector(1, index, value);
    }
    addToObjectiveFromPackedVector(this, vec, false, 1.);
    if (mode_ones > 0) {
      break;
    } else {
      indexAddedToObjective.push_back(check_ind);
    }
  } // prepare first new objective

  if ((mode_ones > 0 && check_ind >= numToCheck - 1)
      || (mode_ones == 0 && indexAddedToObjective.empty())) {
    setConstantObjectiveFromPackedVector(this, 0.);
    this->setRowUpper(point_row_ind, orig_ub);
    if (vec) {
      delete vec;
    }
#ifdef TRACE_CUT_BOUND
    if (tmpSolver) {
      delete tmpSolver;
    }
#endif
    return 1;
  }

  for (int try_ind = 1; try_ind < MAX_NUM_OBJ_PER_POINT; try_ind++) {
    if (tmp_return_code == BAD_RETURN_CODE
        || owner->reachedCutLimit()
        || owner->reachedTimeLimit(timeName, owner->params.get(TIMELIMIT))
        || owner->reachedFailureLimit(num_cuts, num_failures)) {
      ret_val = tmp_return_code != BAD_RETURN_CODE;
      break;
    }
    if (mode_ones > 0 && check_ind >= numToCheck - 1) {
      break;
    }
    tmp_return_code = tryOneObjective(numTimesTightRow, numTimesTightColLB,
        numTimesTightColUB, cuts, origSolver, beta, structSICs, inNBSpace,
        cutHeur);
    if (tmp_return_code > 0) {
      num_cuts_per_point += tmp_return_code;
#ifdef TRACE_CUT_BOUND
      for (int add_ind = 0; add_ind < tmp_return_code; add_ind++) {
        tmpCuts.insert(tmpCuts.rowCut(tmpCuts.sizeCuts() + add_ind));
      }
      applyCutsCustom(tmpSolver, tmpCuts, tmpCuts.sizeCuts() - tmp_return_code,
          tmp_return_code, tmpSolver->getObjValue());
#endif
    }
    if (try_ind + 1 < MAX_NUM_OBJ_PER_POINT) { //&& this->isProvenOptimal()) {
      if (mode_ones > 0) {
        addToObjectiveFromPackedVector(this, vec, false, -1.);
        while (check_ind < numToCheck - 1) {
          check_ind++;
          const int curr_index = indexToCheck[sortIndex[check_ind]];
          if (vec) {
            delete vec;
            vec = NULL;
          }
          if (curr_index < num_rows) {
            if ((mode_ones == 1)
                && (numTimesTightRow[curr_index]
                    != savedNumTimesTightRow[curr_index])) {
              continue;
            }
            vec = new CoinPackedVector(mat->getVectorSize(curr_index),
                mat->getIndices() + mat->getVectorFirst(curr_index),
                mat->getElements() + mat->getVectorFirst(curr_index), false);
          } else if (curr_index < num_rows + num_cols) {
            const int col = curr_index - num_rows;
            if ((mode_ones == 1)
                && (numTimesTightColLB[col] != savedNumTimesTightColLB[col])) {
              continue;
            }
            const int index[1] = { col };
            const double value[1] = { 1. };
            vec = new CoinPackedVector(1, index, value);
          } else {
            const int col = curr_index - num_rows - num_cols;
            if ((mode_ones == 1)
                && (numTimesTightColUB[col] != savedNumTimesTightColUB[col])) {
              continue;
            }
            const int index[1] = { col };
            const double value[1] = { -1. };
            vec = new CoinPackedVector(1, index, value);
          }
          addToObjectiveFromPackedVector(this, vec, false, 1.);
          break;
        }
      } // mode = 0
      else {
        // Remove all rows that are tight but were not before
        int num_removed = 0, num_to_remove = 0;
        for (int tmp_ind = 0; tmp_ind < (int) indexAddedToObjective.size(); tmp_ind++) {
          check_ind = indexAddedToObjective[tmp_ind];
          if (check_ind < 0) {
            continue;
          }
          num_to_remove++;
          const int curr_index = indexToCheck[sortIndex[check_ind]];
          double curr_activity = 0., curr_rhs = 0.;
          if (curr_index < num_rows) {
            curr_activity = this->getRowActivity()[curr_index];
            curr_rhs = this->getRightHandSide()[curr_index];
          } else if (curr_index < num_rows + num_cols) {
            const int var = curr_index - num_rows;
            curr_activity = this->getColSolution()[var];
            curr_rhs = this->getColLower()[var];
          } else {
            const int var = curr_index - num_rows - num_cols;
            curr_activity = this->getColSolution()[var];
            curr_rhs = this->getColUpper()[var];
          }
          if (isVal(curr_activity, curr_rhs)) {
            num_removed++;
            indexAddedToObjective[tmp_ind] += 1;
            indexAddedToObjective[tmp_ind] *= -1;
            if (vec) {
              delete vec;
              vec = NULL;
            }
            if (curr_index < num_rows) {
              vec = new CoinPackedVector(mat->getVectorSize(curr_index),
                  mat->getIndices() + mat->getVectorFirst(curr_index),
                  mat->getElements() + mat->getVectorFirst(curr_index), false);
            } else if (curr_index < num_rows + num_cols) {
              const int col = curr_index - num_rows;
              const int index[1] = { col };
              const double value[1] = { 1. };
              vec = new CoinPackedVector(1, index, value);
            } else {
              const int col = curr_index - num_rows - num_cols;
              const int index[1] = { col };
              const double value[1] = { -1. };
              vec = new CoinPackedVector(1, index, value);
            }
            addToObjectiveFromPackedVector(this, vec, false, -1.);
          }
        }
        if (num_removed * num_to_remove == 0 || (num_to_remove - num_removed == 0)) {
          break;
        }
      } // mode = 1
    } // try_ind = 1 ... MAX_NUM_TRIES
  } // try max_depth times
  setConstantObjectiveFromPackedVector(this, 0.);
  this->setRowUpper(point_row_ind, orig_ub);
  if (vec) {
    delete vec;
  }
#ifdef TRACE_CUT_BOUND
  if (tmpSolver) {
    delete tmpSolver;
  }
#endif
  return ret_val;
} /* findCutsTightOnPoint */

/**
 * @brief Try the (1,...,1) obj, all points, then select among rays ``intelligently''
 * @return Number generated cuts
 */
int PRLP::targetStrongAndDifferentCuts(const double beta, OsiCuts& cuts,
    const OsiSolverInterface* const origSolver, const OsiCuts* const structSICs,
    const std::string& timeName) {
  const bool inNBSpace = owner->params.get(intConst::NB_SPACE);
  if (owner->reachedCutLimit()
      || owner->reachedTimeLimit(timeName, owner->params.get(TIMELIMIT))) {
    return 0;
  }

  // Setup step
  std::vector<rowAndActivity> pointIndex, rayIndex;
  setupForTargetedCutGeneration(pointIndex, rayIndex);

  // numTimesTight will keep how frequently that row (or column) was optimal for some objective
  // If the value is < 0, then we should not ever re-try that objective
  std::vector<int> numTimesTightRow(this->getNumRows(), 0);
  std::vector<int> numTimesTightColLB(this->getNumCols(), 0);
  std::vector<int> numTimesTightColUB(this->getNumCols(), 0);

  int num_points_tried = 0, num_rays_tried = 0;
  int return_code = 0;
  const int BAD_RETURN_CODE = -1 * (static_cast<int>(CglVPC::FailureType::PRIMAL_INFEASIBLE) + 1);
  const int MAX_NUM_POINTS_TO_TRY = owner->params.get(USE_TIGHT_POINTS);
  const int MAX_NUM_OBJ_PER_POINT =
      owner->params.get(intConst::NUM_OBJ_PER_POINT) < 0 ?
          -1 * owner->params.get(intConst::NUM_OBJ_PER_POINT) * owner->getCutLimit() :
          owner->params.get(intConst::NUM_OBJ_PER_POINT);
  const int MAX_NUM_RAYS_TO_TRY =
      owner->params.get(USE_TIGHT_RAYS) < 0 ?
          -1 * owner->params.get(USE_TIGHT_RAYS) * std::ceil(std::sqrt(rayIndex.size())) :
          owner->params.get(USE_TIGHT_RAYS);
  const int MAX_NUM_UNIT_VECTORS_TO_TRY =
      owner->params.get(USE_UNIT_VECTORS) < 0 ?
          -1 * owner->params.get(USE_UNIT_VECTORS) * std::ceil(std::sqrt(owner->probData.num_cols)) :
          CoinMin(static_cast<int>(this->nonZeroColIndex.size()), owner->params.get(USE_UNIT_VECTORS));
  const std::string PRLP_SOLVE_TIME = CglVPC::VPCTimeStatsName[static_cast<int>(CglVPC::VPCTimeStats::PRLP_SOLVE_TIME)];

  // First, the all ones objective
  if (owner->params.get(USE_ALL_ONES) == 1) {
    const CglVPC::ObjectiveType cutHeur = CglVPC::ObjectiveType::ALL_ONES;
    const std::string currTimeName = CglVPC::ObjectiveTypeName[static_cast<int>(cutHeur)] + "_TIME";
    owner->timer.start_timer(currTimeName);
#ifdef TRACE
    printf("\n## Try all ones objective; cuts so far: %d. ##\n",
        cuts.sizeCuts());
#endif
    setConstantObjectiveFromPackedVector(this, 1.);
    return_code = tryOneObjective(numTimesTightRow, numTimesTightColLB,
        numTimesTightColUB, cuts, origSolver, beta, structSICs, inNBSpace,
        cutHeur);
    setConstantObjectiveFromPackedVector(this, 0.);
    owner->timer.end_timer(currTimeName);

    if (return_code == BAD_RETURN_CODE
        || owner->reachedCutLimit()
        || owner->reachedTimeLimit(timeName, owner->params.get(TIMELIMIT))
        || owner->reachedFailureLimit(num_cuts, num_failures)) {
      return num_cuts; // abandon
    }
  } // all ones heuristic

  // Cut away the post-Gomory cut optimum
  if (owner->params.get(USE_ITER_BILINEAR) >= 1) {
    // Try the bilinear optimization; this usually does not add many cuts
    const CglVPC::ObjectiveType cutHeur = CglVPC::ObjectiveType::ITER_BILINEAR;
    const std::string currTimeName = CglVPC::ObjectiveTypeName[static_cast<int>(cutHeur)] + "_TIME";
    owner->timer.start_timer(currTimeName);
#ifdef TRACE
    printf(
        "\n## Using iterative procedure to get deepest cut after adding Gomory cuts; cuts so far: %d. ##\n",
        cuts.sizeCuts());
#endif
//    const int init_num_iter_bil_obj = num_obj_tried;
//    const int init_num_iter_bil_cuts = cuts.sizeCuts();
    iterateDeepestCutPostGomory(cuts, origSolver, beta, structSICs, inNBSpace);
    setConstantObjectiveFromPackedVector(this, 0.);
//    const int final_num_iter_bil_obj = num_obj_tried - init_num_iter_bil_obj;
//    const int final_num_iter_bil_cuts = cuts.sizeCuts() - init_num_iter_bil_cuts;
    owner->timer.end_timer(currTimeName);

    if (return_code == BAD_RETURN_CODE
        || owner->reachedCutLimit()
        || owner->reachedTimeLimit(timeName, owner->params.get(TIMELIMIT))
        || owner->reachedFailureLimit(num_cuts, num_failures)) {
      return num_cuts; // abandon
    }
  } // separate post-Gomory cut optimum

  // Reset things with unlimited time
#ifdef TRACE
  printf(
      "\n## Resetting solver before using subsequent heuristics; cuts so far: %d. ##\n",
      cuts.sizeCuts());
#endif
  owner->timer.start_timer(PRLP_SOLVE_TIME);
  this->getModelPtr()->setNumberIterations(0);
  setTimeLimit(this, -1.);
  this->initialSolve();
  owner->timer.end_timer(PRLP_SOLVE_TIME);

  // Let us try to target the strong branching lower bound
  // This means we find something tight on the point with the lowest objective value
  // Since that point cannot be cut away, that is the best bound we can hope for
  const int disj_lb_row_ind =
      (pointIndex[0].row < 0) ?
          (-1 * (pointIndex[0].row + 1)) : pointIndex[0].row;
  if (owner->params.get(USE_DISJ_LB) > 0) {
    const CglVPC::ObjectiveType cutHeur = CglVPC::ObjectiveType::DISJ_LB;
    const std::string currTimeName = CglVPC::ObjectiveTypeName[static_cast<int>(cutHeur)] + "_TIME";
    owner->timer.start_timer(currTimeName);
#ifdef TRACE
    printf(
        "\n## Try finding a cut tight on the branching lb point; cuts so far: %d. ##\n",
        cuts.sizeCuts());
#endif
    const int goodReturn = findCutsTightOnPoint(numTimesTightRow,
        numTimesTightColLB, numTimesTightColUB, disj_lb_row_ind, cutHeur, beta,
        cuts, origSolver, structSICs, timeName, inNBSpace,
        MAX_NUM_OBJ_PER_POINT);
    owner->timer.end_timer(currTimeName);

    if (!goodReturn || return_code == BAD_RETURN_CODE
        || owner->reachedCutLimit()
        || owner->reachedTimeLimit(timeName, owner->params.get(TIMELIMIT))
        || owner->reachedFailureLimit(num_cuts, num_failures)) {
      return num_cuts; // abandon
    }
  } // branching lb heuristic

  // For each of the points that has not been tight so far, try it
  if (MAX_NUM_POINTS_TO_TRY > 0) {
    bool goodReturn = true;
    const CglVPC::ObjectiveType cutHeur = CglVPC::ObjectiveType::TIGHT_POINTS;
    const std::string currTimeName = CglVPC::ObjectiveTypeName[static_cast<int>(cutHeur)] + "_TIME";
    owner->timer.start_timer(currTimeName);
#ifdef TRACE
    printf(
        "\n## Try finding a cut tight on each intersection point; cuts so far: %d. %d total points to try, %d objectives left. ##\n",
        cuts.sizeCuts(), (int) pointIndex.size(), MAX_NUM_POINTS_TO_TRY);
#endif
    for (int ind = 1; ind < (int) pointIndex.size() - 1; ind++) { // skip the first point (same as the strong lb point)
      if (!goodReturn || return_code == BAD_RETURN_CODE
          || owner->reachedCutLimit()
          || owner->reachedTimeLimit(timeName, owner->params.get(TIMELIMIT))
          || owner->reachedFailureLimit(num_cuts, num_failures)) {
        owner->timer.end_timer(currTimeName);
        return num_cuts; // abandon
      }
      if (num_points_tried >= MAX_NUM_POINTS_TO_TRY) {
        break;
      }
      const int row_ind = pointIndex[ind].row;
      if (row_ind < 0 || numTimesTightRow[row_ind] != 0) {
        continue;
      }
      num_points_tried++;
      pointIndex[ind].row = -1 * (row_ind + 1);
      goodReturn = findCutsTightOnPoint(numTimesTightRow, numTimesTightColLB,
          numTimesTightColUB, disj_lb_row_ind, cutHeur, beta, cuts, origSolver,
          structSICs, timeName, inNBSpace, MAX_NUM_OBJ_PER_POINT);
    } // tight on points
    owner->timer.end_timer(currTimeName);
  } // MAX_NUM_POINTS_TO_TRY > 0

  if (return_code == BAD_RETURN_CODE
      || owner->reachedCutLimit()
      || owner->reachedTimeLimit(timeName, owner->params.get(TIMELIMIT))
      || owner->reachedFailureLimit(num_cuts, num_failures)) {
    return num_cuts; // abandon
  }

  // Next, we should try the nonbasic directions individually, because these have worked well before
  if (MAX_NUM_UNIT_VECTORS_TO_TRY > 0) {
    const CglVPC::ObjectiveType cutHeur = CglVPC::ObjectiveType::UNIT_VECTORS;
    const std::string currTimeName = CglVPC::ObjectiveTypeName[static_cast<int>(cutHeur)] + "_TIME";
    owner->timer.start_timer(currTimeName);
#ifdef TRACE
    printf(
        "\n## Try finding deepest cut with respect to the axis directions; cuts so far: %d. %d dirns, %d to try. ##\n",
        cuts.sizeCuts(), this->getNumCols(), MAX_NUM_UNIT_VECTORS_TO_TRY);
#endif
    const int numNonZero = nonZeroColIndex.size();
    std::vector<int> sortIndex(numNonZero);
    std::vector<double> sortCriterion(numNonZero);
    for (int i = 0; i < numNonZero; i++) {
      sortIndex[i] = i;
      sortCriterion[i] = inNBSpace ? owner->probData.NBReducedCost[nonZeroColIndex[i]] : origSolver->getReducedCost()[nonZeroColIndex[i]];
    }
    std::sort(sortIndex.begin(), sortIndex.end(),
        index_cmp_dsc<const std::vector<double>&>(sortCriterion)); // descending
    int num_unit_vectors_tried = 0;
    for (int tmp_ind = 0; tmp_ind < (int) sortIndex.size(); tmp_ind++) {
      if (num_unit_vectors_tried >= MAX_NUM_UNIT_VECTORS_TO_TRY) {
        break;
      }
      if (return_code == BAD_RETURN_CODE
          || owner->reachedCutLimit()
          || owner->reachedTimeLimit(timeName, owner->params.get(TIMELIMIT))
          || owner->reachedFailureLimit(num_cuts, num_failures)) {
        owner->timer.end_timer(currTimeName);
        return num_cuts; // abandon
      }
      const int nb_dir = sortIndex[tmp_ind];
      if (numTimesTightColLB[nb_dir] != 0) {
        continue;
      }
      this->setObjCoeff(nb_dir, 1.);
      num_unit_vectors_tried++;
      return_code = tryOneObjective(numTimesTightRow, numTimesTightColLB,
          numTimesTightColUB, cuts, origSolver, beta, structSICs, inNBSpace,
          cutHeur);
      if (this->isProvenOptimal()) {
        if (greaterThanVal(this->getColSolution()[nb_dir],
            this->getColLower()[nb_dir])) {
          this->setColLower(nb_dir, this->getObjValue());
        }
      }
      this->setObjCoeff(nb_dir, 0.);
    } // nb directions
    owner->timer.end_timer(currTimeName);
  } // try axis directions in the PRLP

  if (return_code == BAD_RETURN_CODE
      || owner->reachedCutLimit()
      || owner->reachedTimeLimit(timeName, owner->params.get(TIMELIMIT))
      || owner->reachedFailureLimit(num_cuts, num_failures)) {
    return num_cuts; // abandon
  }

  const CoinPackedMatrix* mat = this->getMatrixByRow();
  if (MAX_NUM_RAYS_TO_TRY > 0) {
    // Continue with rays
    const CglVPC::ObjectiveType cutHeur = CglVPC::ObjectiveType::TIGHT_RAYS;
    const std::string currTimeName = CglVPC::ObjectiveTypeName[static_cast<int>(cutHeur)] + "_TIME";
    owner->timer.start_timer(currTimeName);
  #ifdef TRACE
    printf(
        "\n## Try finding a cut tight on the intersection rays; cuts so far: %d: %d total rays to try, %d objectives left. ##\n",
        cuts.sizeCuts(), (int) rayIndex.size(), MAX_NUM_RAYS_TO_TRY);
  #endif
    for (int ind = 0; ind < (int) rayIndex.size(); ind++) {
      if (num_rays_tried >= MAX_NUM_RAYS_TO_TRY) {
        break;
      }
      if (return_code == BAD_RETURN_CODE
          || owner->reachedCutLimit()
          || owner->reachedTimeLimit(timeName, owner->params.get(TIMELIMIT))
          || owner->reachedFailureLimit(num_cuts, num_failures)) {
        owner->timer.end_timer(currTimeName);
        return num_cuts; // abandon
      }
      const int curr_index = ind; //curr_ray_info.index
      const int row_ind = rayIndex[curr_index].row;
      if (row_ind < 0 || numTimesTightRow[row_ind] != 0) {
        continue;
      }
      rayIndex[curr_index].row = -1 * (row_ind + 1);
      const CoinShallowPackedVector vec = mat->getVector(row_ind);
      addToObjectiveFromPackedVector(this, &vec, false);
      num_rays_tried++;
      return_code = tryOneObjective(numTimesTightRow, numTimesTightColLB,
          numTimesTightColUB, cuts, origSolver, beta, structSICs, inNBSpace,
          cutHeur);
      if (return_code < 0) {
        numTimesTightRow[row_ind]++;
        numTimesTightRow[row_ind] *= -1;
      }
      setConstantObjectiveFromPackedVector(this, 0., vec.getNumElements(),
          vec.getIndices());
    } // tight on rays
    owner->timer.end_timer(currTimeName);
  } // MAX_NUM_RAYS_TO_TRY > 0

  if (return_code == BAD_RETURN_CODE
      || owner->reachedCutLimit()
      || owner->reachedTimeLimit(timeName, owner->params.get(TIMELIMIT))
      || owner->reachedFailureLimit(num_cuts, num_failures)) {
    return num_cuts; // abandon
  }

  if (num_points_tried < MAX_NUM_POINTS_TO_TRY) {
    // We may still not hit the cut limit but have exhausted all points and rays
    // In this case, we will try the points then the rays in the order they are in
    // (sorted by objective value)
    // and hope that we get new cuts along the way?
    bool goodReturn = true;
    const CglVPC::ObjectiveType cutHeur = CglVPC::ObjectiveType::TIGHT_POINTS2;
    const std::string currTimeName = CglVPC::ObjectiveTypeName[static_cast<int>(cutHeur)] + "_TIME";
    owner->timer.start_timer(currTimeName);
#ifdef TRACE
    printf(
        "\n## Try being tight on the intersection points _again_; cuts so far: %d: %d total points to try, %d objectives left. ##\n",
        cuts.sizeCuts(), (int) pointIndex.size(), MAX_NUM_POINTS_TO_TRY - num_points_tried);
#endif
    for (int ind = 0; ind < (int) pointIndex.size(); ind++) {
      if (!goodReturn || return_code == BAD_RETURN_CODE
          || owner->reachedCutLimit()
          || owner->reachedTimeLimit(timeName, owner->params.get(TIMELIMIT))
          || owner->reachedFailureLimit(num_cuts, num_failures)) {
        owner->timer.end_timer(currTimeName);
        return num_cuts; // abandon
      }
      if (num_points_tried >= MAX_NUM_POINTS_TO_TRY) {
        break;
      }
      const int row_ind =
        (pointIndex[ind].row < 0) ?
        -1 * (pointIndex[ind].row + 1) : pointIndex[ind].row;
      if (numTimesTightRow[row_ind] < 0) {
        continue;
      }
      if (row_ind == disj_lb_row_ind) {
        continue;
      }
      num_points_tried++;
      goodReturn = findCutsTightOnPoint(numTimesTightRow, numTimesTightColLB,
          numTimesTightColUB, disj_lb_row_ind, cutHeur, beta, cuts, origSolver,
          structSICs, timeName, inNBSpace, MAX_NUM_OBJ_PER_POINT);
    } // tight on points
    owner->timer.end_timer(currTimeName);
  } // any objectives left for points?

  if (return_code == BAD_RETURN_CODE
      || owner->reachedCutLimit()
      || owner->reachedTimeLimit(timeName, owner->params.get(TIMELIMIT))
      || owner->reachedFailureLimit(num_cuts, num_failures)) {
    return num_cuts; // abandon
  }

  if (num_rays_tried < MAX_NUM_RAYS_TO_TRY) {
    // Continue with rays
    const CglVPC::ObjectiveType cutHeur = CglVPC::ObjectiveType::TIGHT_RAYS2;
    const std::string currTimeName = CglVPC::ObjectiveTypeName[static_cast<int>(cutHeur)] + "_TIME";
    owner->timer.start_timer(currTimeName);
#ifdef TRACE
    printf(
        "\n## Try the intersection rays _again_; cuts so far: %d: %d total rays to try, %d objectives left. ##\n",
        cuts.sizeCuts(), (int) rayIndex.size(), MAX_NUM_RAYS_TO_TRY - num_rays_tried);
#endif
    for (int ind = 0; ind < (int) rayIndex.size(); ind++) {
      if (return_code == BAD_RETURN_CODE
          || owner->reachedCutLimit()
          || owner->reachedTimeLimit(timeName, owner->params.get(TIMELIMIT))
          || owner->reachedFailureLimit(num_cuts, num_failures)) {
        owner->timer.end_timer(currTimeName);
        return num_cuts; // abandon
      }
      if (num_rays_tried >= MAX_NUM_RAYS_TO_TRY) {
        break;
      }
      const int row_ind =
        (rayIndex[ind].row < 0) ?
        -1 * (rayIndex[ind].row + 1) : rayIndex[ind].row;
      if (numTimesTightRow[row_ind] < 0) {
        continue;
      }
      const CoinShallowPackedVector vec = mat->getVector(row_ind);
      addToObjectiveFromPackedVector(this, &vec, false);
      num_rays_tried++;
      return_code = tryOneObjective(numTimesTightRow, numTimesTightColLB,
          numTimesTightColUB, cuts, origSolver, beta, structSICs, inNBSpace,
          cutHeur);
      setConstantObjectiveFromPackedVector(this, 0., vec.getNumElements(), vec.getIndices());
    } // tight on rays again
    owner->timer.end_timer(currTimeName);
  } // any objectives left for rays?

  return num_cuts;
} /* targetStrongAndDifferentCuts */

/**
 * We want to solve the bilinear program
 * min <\alpha, x>
 * s.t.
 * <\alpha, p> \ge 1 \forall p \in \pointset
 * <\alpha, r> \ge 0 \forall r \in \rayset
 * x \in Pbar (P + all SICs)
 */
int PRLP::iterateDeepestCutPostGomory(OsiCuts & cuts,
    const OsiSolverInterface* const origSolver, const double beta,
    const OsiCuts* structSICs, const bool inNBSpace) {
  if (owner->reachedCutLimit()
      || (owner->params.get(USE_ITER_BILINEAR) <= 0)) {
    return 0;
  }

  int num_cuts_gen = 0;
  const int MAX_NUM_ITER_BIL_CUTS = owner->params.get(intParam::USE_ITER_BILINEAR);
  const int ITER_PER_CUT = 1;
  const CglVPC::ObjectiveType cutHeur = CglVPC::ObjectiveType::ITER_BILINEAR;
  bool stationary_point_flag = false;
  const double initSolveTime =
      owner->timer.get_time(
          CglVPC::VPCTimeStatsName[static_cast<int>(CglVPC::VPCTimeStats::INIT_SOLVE_TIME)]);
  const double timeLimit = owner->params.get(PRLP_TIMELIMIT);

  SolverInterface* PostGomorySolver = dynamic_cast<SolverInterface*>(origSolver->clone());
  setLPSolverParameters(PostGomorySolver, owner->params.get(VERBOSITY), 2 * initSolveTime);

  // Add GMICs
  OsiCuts GMICs;
  if (structSICs == NULL) {
    PostGomorySolver->initialSolve();
    CglGMI GMIGen;
    GMIGen.generateCuts(*PostGomorySolver, GMICs);
  } else {
    GMICs = *structSICs;
  }
  applyCutsCustom(PostGomorySolver, GMICs);
  if (hitTimeLimit(PostGomorySolver) || PostGomorySolver->isIterationLimitReached()) {
    return exitFromIterateDeepestCutPostGomory(num_cuts_gen, PostGomorySolver);
  }
  if (!PostGomorySolver->isProvenOptimal()) {
    error_msg(error_str, "Adding Gomory cuts causes LP (with original objective) to not be optimal. This cannot happen!\n"); // NB: could actually happen in subspace
    exit(1);
  }
  setTimeLimit(PostGomorySolver, timeLimit);

  // Save this post-Gomory point as the objective for the PRLP
  setObjectiveFromStructuralPoint(PostGomorySolver->getColSolution(), NULL, origSolver, inNBSpace);
  num_cuts_gen = genCut(cuts, origSolver, beta, structSICs, inNBSpace, cutHeur);
  if (num_cuts_gen < 0) {
    if (PostGomorySolver) {
      delete PostGomorySolver;
    }
    return 0;
  }
  if (MAX_NUM_ITER_BIL_CUTS == 1) {
    return exitFromIterateDeepestCutPostGomory(num_cuts_gen, PostGomorySolver);
  }

  // We may need to get primal rays
  // In case of dual infeasibility, there should be at least one
  // During this procedure, since we are using dual rays, do not use presolve
  this->setHintParam(OsiDoPresolveInInitial, false);
  this->setHintParam(OsiDoPresolveInResolve, false);
  PostGomorySolver->setHintParam(OsiDoPresolveInInitial, false);
  PostGomorySolver->setHintParam(OsiDoPresolveInResolve, false);

  std::vector<double> prevSolution(PostGomorySolver->getColSolution(),
      PostGomorySolver->getColSolution() + PostGomorySolver->getNumCols());
  OsiRowCut structSpaceCut;
  int i, j;
  for (j = 1; j < MAX_NUM_ITER_BIL_CUTS; j++) {
    for (i = 0; i < ITER_PER_CUT; i++) {
      if (owner->reachedCutLimit()
          || owner->reachedFailureLimit(num_cuts, num_failures)) {
        break;
      }

      /***********
       *  STEP 1
       ***********/
      { // First set the objective for PostGomorySolver using the PRLP solution
        std::vector<double*> PRLPPrimalRay;
        const double* sol = NULL;
        if (this->isProvenOptimal()) {
          sol = this->getColSolution();
        } else if (this->isProvenDualInfeasible()) {
          PRLPPrimalRay = this->getPrimalRays(1); // there may be more than one, which all lead to cuts

          // It has happened (e.g., instance rlp2)
          if (PRLPPrimalRay[0] == NULL) {
            this->getModelPtr()->setNumberIterations(0);
            setTimeLimit(this, timeLimit);
            this->initialSolve();
            checkSolverOptimality(this, false, timeLimit);
            if (this->isProvenDualInfeasible()) {
              // Free memory
              if (PRLPPrimalRay.size() > 0 && PRLPPrimalRay[0]) {
                delete[] PRLPPrimalRay[0];
                PRLPPrimalRay.clear();
              }
              PRLPPrimalRay = this->getPrimalRays(1);
            }
            if (!this->isProvenDualInfeasible() || !PRLPPrimalRay[0]) {
              warning_msg(warnstring,
                  "Cut %d, iter %d: "
                  "Issues obtaining ray when PRLP is unbounded when trying to find deepest post-Gomory cut. "
                  "Abandoning this function.\n",
                  j, i);
              // Free memory
              if (PRLPPrimalRay.size() > 0 && PRLPPrimalRay[0]) {
                delete[] PRLPPrimalRay[0];
                PRLPPrimalRay.clear();
              }
              return exitFromIterateDeepestCutPostGomory(num_cuts_gen, PostGomorySolver);
            }
          }

  //        scalePrimalRay(this->getNumCols(), &PRLPPrimalRay[0]); // amk 2019-02-26: should reimplement
          sol = PRLPPrimalRay[0];
        } // check if PRLP is bounded

        // Recall that PostGomorySolver is in the structural space,
        // while the PRLP is potentially in the NB space
        if (sol) {
          if (inNBSpace) {
            setCutFromNBCoefficients(&structSpaceCut, this->nonZeroColIndex,
                sol, beta, origSolver, owner->probData.NBVarIndex, true,
                owner->probData.EPS);
            addToObjectiveFromPackedVector(PostGomorySolver,
                &structSpaceCut.row(), true);
          } else {
            if ((int) nonZeroColIndex.size() > 0) {
              for (int i = 0; i < (int) nonZeroColIndex.size(); i++) {
                PostGomorySolver->setObjCoeff(nonZeroColIndex[i],
                    this->getColSolution()[i]);
              }
            } else {
              PostGomorySolver->setObjective(this->getColSolution());
            }
          }
        } // check that sol exists

        // Free memory
        if (PRLPPrimalRay.size() > 0 && PRLPPrimalRay[0]) {
          delete[] PRLPPrimalRay[0];
          PRLPPrimalRay.clear();
        }
      } // set PostGomorySolver objective from the PRLP solution

      /***********
       *  STEP 2
       ***********/
      { // Resolve PostGomorySolver so that we can extract the new solution from it
        setTimeLimit(PostGomorySolver, timeLimit);
        PostGomorySolver->resolve();
        checkSolverOptimality(PostGomorySolver, false, timeLimit); // with pk1, it has happened that the PostGomorySolver is dual infeasible but does not say it is proven so
        if (!PostGomorySolver->isProvenOptimal()
            && !PostGomorySolver->isProvenDualInfeasible()) {
          error_msg(errorstring,
              "Cut %d, iter %d: Iterative solver not proven optimal or dual infeasible.\n",
              j, i);
          writeErrorToLog(errorstring, owner->params.logfile);
          exit(1);
        }

        if (isNegInfinity(PostGomorySolver->getObjValue(), -owner->params.get(doubleConst::INF))) {
          // We essentially have a ray, though it is not being reported as such
          // Try resolving because this fixes the issue sometimes
          setTimeLimit(PostGomorySolver, timeLimit);
          PostGomorySolver->initialSolve();
          checkSolverOptimality(PostGomorySolver, false, timeLimit);

          if (!PostGomorySolver->isProvenDualInfeasible()) {
            warning_msg(warnstring,
                "Cut %d, iter %d: Scaling issues when trying to find deepest post-Gomory cut. "
                "Abandoning this function.\n",
                j, i);
            return exitFromIterateDeepestCutPostGomory(num_cuts_gen, PostGomorySolver);
          }
        }
      } // get new PostGomorySolver solution

      /***********
       *  STEP 3
       ***********/
      { // Set PRLP objective from PostGomorySolution
        const double* sol = NULL;
        std::vector<double*> PostGomorySolverPrimalRay;
        if (PostGomorySolver->isProvenOptimal()) {
          sol = PostGomorySolver->getColSolution();
        } else if (PostGomorySolver->isProvenDualInfeasible()) {
          // Obtain the dual ray
          PostGomorySolverPrimalRay = PostGomorySolver->getPrimalRays(1);

          // It has happened (e.g., instance rlp2)
          if (PostGomorySolverPrimalRay[0] == NULL) {
            setTimeLimit(PostGomorySolver, timeLimit);
            PostGomorySolver->initialSolve();
            checkSolverOptimality(PostGomorySolver, false, timeLimit);
            if (PostGomorySolver->isProvenDualInfeasible()) {
              // Free memory
              if (PostGomorySolverPrimalRay.size() > 0 && PostGomorySolverPrimalRay[0]) {
                delete[] PostGomorySolverPrimalRay[0];
                PostGomorySolverPrimalRay.clear();
              }
              PostGomorySolverPrimalRay = PostGomorySolver->getPrimalRays(1);
            }
            if (!PostGomorySolver->isProvenDualInfeasible()
                || !PostGomorySolverPrimalRay[0]) {
              warning_msg(warnstring,
                  "Cut %d, iter %d: "
                  "Issues obtaining ray when PostGomorySolver is unbounded while getting deepest post-Gomory cut. "
                  "Abandoning this function.\n",
                  j, i);
              // Free memory
              if (PostGomorySolverPrimalRay.size() > 0 && PostGomorySolverPrimalRay[0]) {
                delete[] PostGomorySolverPrimalRay[0];
                PostGomorySolverPrimalRay.clear();
              }
              return exitFromIterateDeepestCutPostGomory(num_cuts_gen, PostGomorySolver);
            }
          } // exit early if issues arise getting the primal ray
  //        scalePrimalRay(PostGomorySolver->getNumCols(), &SICSolverPrimalRay[0]); // NB: again, should reimplement
          sol = PostGomorySolverPrimalRay[0];
        } // check whether PostGomorySolver is optimal or unbdd

        // Check if we have reached a stationary point
        const double parallelism = getParallelism(
            PostGomorySolver->getNumCols(), prevSolution.data(), sol);

        // Reset prevSolution to be current solution, so we can set objective
        prevSolution.assign(sol, sol + PostGomorySolver->getNumCols());

        // Free memory
        if (PostGomorySolverPrimalRay.size() > 0 && PostGomorySolverPrimalRay[0]) {
          delete[] PostGomorySolverPrimalRay[0];
          PostGomorySolverPrimalRay.clear();
        }

        // Break if we reached the stationary point (after freeing memory)
        if (isVal(parallelism, 1.)) {
          stationary_point_flag = true;
          break;
        }
      } // set PRLP objective from PostGomorySolver solution

      /***********
       *  STEP 4
       ***********/
      // Resolve PRLP without generating a cut
      setObjectiveFromStructuralPoint(prevSolution.data(), NULL, origSolver, inNBSpace);
      if (resolvePRLP() < 0) {
        warning_msg(warnstring,
            "Cut %d, iter %d: Issues resolving PRLP when finding deepest post-Gomory cut. "
            "Abandoning this function.\n",
            j, i);
        return exitFromIterateDeepestCutPostGomory(num_cuts_gen,
            PostGomorySolver);
      }
    } // possibly do several iterations before trying to get a cut again

    if (i == 0) break; // no iterations done
    const int return_code = genCut(cuts, origSolver, beta, structSICs, inNBSpace, cutHeur);
    if (return_code >= 0) {
      num_cuts_gen += return_code;
    } else {
      return exitFromIterateDeepestCutPostGomory(num_cuts_gen, PostGomorySolver);
    }
    if (stationary_point_flag)
      break;
  } // loop at most MAX_NUM_ITER_BIL_CUTS times

  return exitFromIterateDeepestCutPostGomory(num_cuts_gen, PostGomorySolver);
} /* iterateDeepestCutPostMSIC */

int PRLP::exitFromIterateDeepestCutPostGomory(const int num_cuts_gen,
    OsiSolverInterface* const PostGomorySolver) {
  // Delete solvers
  if (PostGomorySolver) {
    delete PostGomorySolver;
  }

  this->getModelPtr()->setNumberIterations(0);
  this->getModelPtr()->setMaximumSeconds(owner->params.get(PRLP_TIMELIMIT));

  this->setHintParam(OsiDoPresolveInInitial, owner->params.get(intConst::PRLP_PRESOLVE) >= 1);
  this->setHintParam(OsiDoPresolveInResolve, owner->params.get(intConst::PRLP_PRESOLVE) >= 2);

  return num_cuts_gen;
} /* exitFromIterateDeepestCutPostGomory */
