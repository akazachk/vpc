/**
 * @file CglVPC.cpp
 * @author A. M. Kazachkov
 * @date 2018-12-24
 */
#include "CglVPC.hpp"

#include <cmath> // abs, floor, ceil
#include <limits> // numeric_limits
#include <algorithm> // std::min_element, std::max_element

// COIN-OR files
#include <CbcModel.hpp>
#include <CglGMI.hpp>

// Project files
#include "CutHelper.hpp" // badViolation
#include "Disjunction.hpp"
#include "PRLP.hpp"
#include "SolverHelper.hpp"
#include "SolverInterface.hpp"
#include "nbspace.hpp"
#include "utility.hpp"
#include "VPCEventHandler.hpp"
using namespace VPCParametersNamespace;

// Various pre-built disjunction options
#include "PartialBBDisjunction.hpp"
#include "SplitDisjunction.hpp"

#ifdef TRACE
#include "vpc_debug.hpp"
#endif

/// Match DisjExitReason status to ExitReason status
CglVPC::ExitReason matchStatus(const DisjExitReason status) {
  if (status == DisjExitReason::SUCCESS_EXIT)
    return CglVPC::ExitReason::SUCCESS_EXIT;
  else if (status == DisjExitReason::OPTIMAL_SOLUTION_FOUND_EXIT)
    return CglVPC::ExitReason::OPTIMAL_SOLUTION_FOUND_EXIT;
  else if (status == DisjExitReason::TOO_FEW_TERMS_EXIT)
    return CglVPC::ExitReason::TOO_FEW_TERMS_EXIT;
  else if (status == DisjExitReason::NO_DISJUNCTION_EXIT)
    return CglVPC::ExitReason::NO_DISJUNCTION_EXIT;
  else
    return CglVPC::ExitReason::UNKNOWN;
} /* matchStatus */

const std::vector<std::string> CglVPC::ExitReasonName {
  "SUCCESS",
  "CUT_LIMIT",
  "FAIL_LIMIT",
  "OPTIMAL_SOLUTION_FOUND",
  "PRLP_INFEASIBLE",
  "PRLP_NUMERICAL_ISSUES",
  "PRLP_TIME_LIMIT",
  "TIME_LIMIT",
  "TOO_FEW_TERMS",
  "NO_CUTS_LIKELY",
  "NO_DISJUNCTION",
  "UNKNOWN"
}; /* ExitReasonName */
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
const std::vector<std::string> CglVPC::ObjectiveTypeName {
  "DUMMY_OBJ",
  "ALL_ONES",
  "CUT_VERTICES",
  "ITER_BILINEAR",
  "UNIT_VECTORS",
  "DISJ_LB",
  "TIGHT_POINTS",
  "TIGHT_RAYS",
  "TIGHT_POINTS2",
  "TIGHT_RAYS2",
  "USER",
  "OBJ_CUT",
  "ONE_SIDED"
}; /* ObjectiveTypeName */
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
  "DLB_EQUALS_DUB_NO_OBJ",
  "DLB_EQUALS_LPOPT_NO_OBJ",
  "PRIMAL_INFEASIBLE_NO_OBJ",
  "NUMERICAL_ISSUES_NO_OBJ",
  "UNKNOWN"
}; /* FailureTypeName */
const std::string CglVPC::time_T1 = "TIME_TYPE1_";
const std::string CglVPC::time_T2 = "TIME_TYPE2_";

/**
 * @details This allows us to change between restricting number of cuts per split and total number of cuts easily
 * The cut limit is either across all cut-generating sets
 * or it is divided among the cut-generating sets (either as the limit / numFracVars, or as a fixed number per cgs)
 * If CUTLIMIT = 0 => no cuts
 * If CUTLIMIT > 0 => absolute cut limit
 * If CUTLIMIT < 0 => cut limit per cgs
 */
int CglVPC::getCutLimit(const int CUTLIMIT, const int numFracVar) {
  if (CUTLIMIT == 0) {
    return 0;
    //return std::numeric_limits<int>::max();
  } else if (CUTLIMIT > 0) {
    return CUTLIMIT;
  } else {
    return (numFracVar <= 0) ? 0 : std::ceil((-1. * CUTLIMIT * numFracVar));
  }
} /* getCutLimit */

int CglVPC::getCutLimit() const {
  return params.get(CUTLIMIT);
} /* getCutLimit */

/// Calls initialize()
CglVPC::CglVPC() {
  initialize();
} /* default constructor */

/// Calls \link initialize() initialize(NULL, &param) \endlink
CglVPC::CglVPC(const VPCParameters& param, const int round_ind) {
  initialize(NULL, &param);
  this->num_rounds = round_ind;
} /* param constructor */

/// Calls CglCutGenerator(source) and \link initialize() initialize(&source) \endlink
CglVPC::CglVPC(const CglVPC& source) : CglCutGenerator(source) {
  initialize(&source);
} /* copy constructor */

/// Deletes #disjunction and #prlp
CglVPC::~CglVPC() {
  if (disjunction && ownsDisjunction) { delete disjunction; }
  if (prlp) delete prlp;
} /* destructor */

/// Calls \link initialize() initialize(&source) \endlink
CglVPC& CglVPC::operator=(const CglVPC& source) {
  if (this != &source) {
    initialize(&source);
  }
  return *this;
} /* assignment operator */

/// Creates new CglVPC instance via copy constructor \link CglVPC(const CglVPC&) CglVPC(*this) \endlink
CglCutGenerator* CglVPC::clone() const {
  return new CglVPC(*this);
} /* clone */

void CglVPC::setParams(const VPCParameters& param) {
  this->params = param;
} /* setParams */

/**
 * @details #disjunction will be deleted if it is currently not NULL and #ownsDisjunction is true
 */
void CglVPC::setDisjunction(
    /// [in] Disjunction to be used; will be cloned via Disjunction::clone() if \p ownIt is true
    Disjunction* const sourceDisj,
    /// [in] -1, 0, or 1; -1 means use #ownsDisjunction value from class
    int ownIt) {
  if (ownIt < 0) {
    ownIt = this->ownsDisjunction;
  }
  if (this->disjunction && this->ownsDisjunction)
    delete disjunction;
  this->ownsDisjunction = ownIt;
  if (ownIt)
    this->disjunction = sourceDisj->clone();
  else
    this->disjunction = sourceDisj;
} /* setDisjunction */

void CglVPC::setUserObjectives(
    /// [in] set of user objectives in structural space, saved in #user_objectives
    const std::vector<std::vector<double> >& obj) {
  this->user_objectives.reserve(obj.size());
  for (const auto& v : obj) {
    this->user_objectives.push_back(v);
  }
} /* setUserObjectives */

void CglVPC::setUserTightPoints(
    /// [in] set of points in structural space, for which the VPCs should be as tight as possible
    const std::vector<std::vector<double> >& tight_points) {
  this->user_tight_points.reserve(tight_points.size());
  for (const auto& v : tight_points) {
    this->user_tight_points.push_back(v);
  }
} /* setUserTightPoints */

void CglVPC::generateCuts(const OsiSolverInterface& si, OsiCuts& cuts, const CglTreeInfo info) {
  // Time starts here, and will end when finish is called
  timer.start_timer(VPCTimeStatsName[static_cast<int>(VPCTimeStats::TOTAL_TIME)]);

  CglVPC::ExitReason status = CglVPC::ExitReason::UNKNOWN;
  if (params.get(intParam::DISJ_TERMS) == 0) {
    status = CglVPC::ExitReason::NO_DISJUNCTION_EXIT;
    finish(status);
    return;
  }
  if (reachedTimeLimit(VPCTimeStats::TOTAL_TIME, params.get(TIMELIMIT))) {
    status = CglVPC::ExitReason::TIME_LIMIT_EXIT;
    finish(status);
    return;
  }

  // Reset things in preparation for round of cuts, in case a previous round was done using this generator
  setupAsNew();
  init_num_cuts = cuts.sizeCuts();
  if (init_num_cuts == 0) {
    this->cutType.resize(0);
    this->objType.resize(0);
  }
  else if (this->canReplaceGivenCuts) {
    // If we are going to be able to replace given cuts,
    // then it must be the case that the current cutType and cutHeurVec
    // should correspond to those old cuts
    const int num_old_cut_type = cutType.size();
    if (num_old_cut_type != init_num_cuts) {
      error_msg(errorstring,
          "# given cuts: %d. # old cuts: %d. "
          "Either set canReplaceGivenCuts to false, or ensure that old cuts are "
          "accounted for in the cutType and cutHeurVec members.\n",
          init_num_cuts, num_old_cut_type);
      writeErrorToLog(errorstring, params.logfile);
      exit(1);
    }
  }

  // Set cut limit
  params.set(CUTLIMIT,
      CglVPC::getCutLimit(params.get(CUTLIMIT),
          si.getFractionalIndices().size()));
  if (reachedCutLimit() && !use_temp_option(std::abs(params.get(intParam::TEMP)), TempOptions::GEN_TIKZ_STRING)) {
    if (si.getFractionalIndices().size() > 0) {
      status = CglVPC::ExitReason::CUT_LIMIT_EXIT;
    } else {
      status = CglVPC::ExitReason::NO_DISJUNCTION_EXIT;
    }
    finish(status);
    return;
  }

  if (!disjunction) {
    mode = static_cast<VPCMode>(params.get(MODE));
    if (mode == VPCMode::CUSTOM) {
      error_msg(errorstring,
          "Mode chosen is CUSTOM but no disjunction is set.\n");
      writeErrorToLog(errorstring, params.logfile);
      exit(1);
    }
    ownsDisjunction = true;
    if (mode == VPCMode::PARTIAL_BB) {
      if (params.get(intParam::DISJ_TERMS) < 2) {
        status = CglVPC::ExitReason::NO_DISJUNCTION_EXIT;
        finish(status);
        return;
      }
      printf(
          "\n## Starting VPC generation from partial branch-and-bound tree with up to %d disjunctive terms. ##\n",
          params.get(intParam::DISJ_TERMS));
      disjunction = new PartialBBDisjunction(this->params);
      dynamic_cast<PartialBBDisjunction*>(disjunction)->num_rounds = this->num_rounds;
      dynamic_cast<PartialBBDisjunction*>(disjunction)->timer = &timer;
    } else if (mode == VPCMode::SPLITS) {
      printf("\n## Starting VPC generation from one split. ##\n");
      disjunction = new SplitDisjunction(this->params);
      dynamic_cast<SplitDisjunction*>(disjunction)->timer = &timer;
    } else {
      error_msg(errorstring,
          "Mode that is chosen has not yet been implemented for VPC generation: %s.\n",
          VPCModeName[static_cast<int>(mode)].c_str());
      writeErrorToLog(errorstring, params.logfile);
      exit(1);
    }
  } else {
    mode = VPCMode::CUSTOM;
  }

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

  // Make sure we are doing a minimization problem; this is just to make later
  // comparisons simpler (i.e., a higher LP obj after adding the cut is better).
  ensureMinimizationObjective(solver);

  if (mode != VPCMode::CUSTOM) {
    // Get disjunctive terms and obtain their optimal bases
    // (If mode is custom, i.e., disjunction is given to us,
    // then we assume that the disjunction is already prepared)
    DisjExitReason disjstatus = disjunction->prepareDisjunction(solver);
    status = matchStatus(disjstatus);
    if (status == CglVPC::ExitReason::OPTIMAL_SOLUTION_FOUND_EXIT) {
      warning_msg(warnstr,
          "An integer (optimal) solution with value %.6g was found prior while getting disjunction.\n",
          disjunction->integer_obj);
          //" We will generate between n and 2n cuts, restricting the value of each variable.\n");
      /*const double* solution = disjunction->integer_sol.data();
      if (solution) {
        for (int col = 0; col < solver->getNumCols(); col++) {
          const double val = solution[col];

          // Check which of the bounds needs to be fixed
          for (int b = 0; b < 2; b++) {
            if ((b == 0 && greaterThanVal(val, solver->getColLower()[col]))
                || (b == 1 && lessThanVal(val, solver->getColUpper()[col]))) {
              const double mult = (b == 0) ? 1. : -1.;
              const double el = mult * 1.;

              OsiRowCut currCut;
              currCut.setLb(mult * val);
              currCut.setRow(1, &cdol, &el, false);
              addCut(currCut, cuts, CutType::OPTIMALITY_CUT,
                  ObjectiveType::ONE_SIDED);
            }
          }
        } // iterate over columns and add optimality cut if needed
      }*/

      // Add objective cut
      OsiRowCut objCut;
      std::vector<int> indices;
      std::vector<double> vals;
      indices.reserve(solver->getNumCols());
      vals.reserve(solver->getNumCols());
      for (int i = 0; i < solver->getNumCols(); i++) {
        if (!isVal(solver->getObjCoefficients()[i], probData.EPS)) {
          indices.push_back(i);
          vals.push_back(solver->getObjCoefficients()[i]);
        }
      }
      objCut.setRow(indices.size(), indices.data(), vals.data(), false);
      double offset = 0.;
      solver->getDblParam(OsiDblParam::OsiObjOffset, offset);
      objCut.setLb(disjunction->integer_obj + offset);
      addCut(objCut, cuts, CutType::OPTIMALITY_CUT,
          ObjectiveType::OBJ_CUT);
    } // exit out early if integer-optimal solution found
    if (status != CglVPC::ExitReason::SUCCESS_EXIT) {
      finish(status);
      return;
    }
  }

  // Make a copy of the solver to allow for fixing variables and changed bounds at root
  SolverInterface* vpcsolver = dynamic_cast<SolverInterface*>(si.clone());
  
  // Save the V-polyhedral relaxations of each optimal basis in the terms of the PRLP
  status = setupConstraints(vpcsolver, cuts);
  if (status != CglVPC::ExitReason::SUCCESS_EXIT) {
    if (vpcsolver) { delete vpcsolver; }
    finish(status);
    return;
  }

  // Generate cuts from the PRLP
  status = tryObjectives(cuts, vpcsolver, NULL);
  if (vpcsolver) { delete vpcsolver; }
  finish(status);
} /* generateCuts */

/**
 * \return Whether cut was added
 */
bool CglVPC::addCut(
    /// [in] Cut being added
    const OsiRowCut& cut,
    /// [in,out] Current set of cuts
    OsiCuts& cuts,
    /// [in] Track the type of cut, in cases we are taking a multimodal approach, e.g., different types of disjunctions
    const CutType& type,
    /// [in] Track the cut heuristic generating this cut
    const ObjectiveType& cutHeur,
    /// [in] Solver pointer, the solution to which will be used to check violation of the cut
    const OsiSolverInterface* const origSolver,
    /// [in] Whether to check violation (when adding one-sided cuts, we check, whereas the other VPCs have already been vetted)
    const bool check_violation) {
  if (check_violation) {
    if (badViolation(&cut, origSolver,
          params.get(VPCParametersNamespace::doubleConst::MIN_VIOL_ABS),
          params.get(VPCParametersNamespace::doubleConst::MIN_VIOL_REL))) {
      return false;
    }
  }
  cuts.insert(cut);
  cutType.push_back(type);
  numCutsOfType[static_cast<int>(type)]++;
  objType.push_back(cutHeur);
  numCutsFromHeur[static_cast<int>(cutHeur)]++;
  num_cuts++;
  return true;
} /* addCut */

//<!------------------ PROTECTED ---------------------->/

/**
 * @details Reset _some_ things (those corresponding to a previous run of this generator).
 * E.g., we do not reset timing, the cutType vector, or the cutHeurVec.
 * The latter two should not be changed and need to correspond to the cuts passed into generateCuts
 * (in order to enable replacing the cuts in PRLP).
 */
void CglVPC::setupAsNew() {
  this->exitReason = CglVPC::ExitReason::UNKNOWN;
  this->numCutsOfType.clear();
  this->numCutsOfType.resize(static_cast<int>(CutType::NUM_CUT_TYPES), 0);
  this->numCutsFromHeur.clear();
  this->numCutsFromHeur.resize(static_cast<int>(ObjectiveType::NUM_OBJECTIVE_TYPES), 0);
  this->numObjFromHeur.clear();
  this->numObjFromHeur.resize(static_cast<int>(ObjectiveType::NUM_OBJECTIVE_TYPES), 0);
  this->numFailsFromHeur.clear();
  this->numFailsFromHeur.resize(static_cast<int>(ObjectiveType::NUM_OBJECTIVE_TYPES), 0);
  this->numFails.clear();
  this->numFails.resize(static_cast<int>(FailureType::NUM_FAILURE_TYPES), 0);
  this->init_num_cuts = 0;
  this->num_cuts = 0;
  this->num_obj_tried = 0;
  this->num_failures = 0;
  this->probData.EPS = this->params.get(EPS);
  if (!this->isSetupForRepeatedUse) {
    this->canReplaceGivenCuts = false;
    this->cutType.resize(0);
    this->objType.resize(0);
    this->num_rounds = 0;
  }
  this->user_objectives.resize(0);
  this->user_tight_points.resize(0);
} /* setupAsNew */

/**
 * @details
 * With neither argument provided, sets up default values for class members through calling #setupAsNew, setting #mode = Partial_BB (in #VPCMode), and initializing #timer.
 * Calls setParams(const VPCParametersNamespace::VPCParameters&) if \p param != NULL, copies over information from \p source if provided.
 */
void CglVPC::initialize(
  /// [in] Everything will be copied from here if provided, except #param when \p param != NULL
  const CglVPC* const source,
  /// [in] Set of parameters to use for the run
  const VPCParameters* const param) {
  if (param != NULL)
    setParams(*param);
  if (source != NULL) {
    if (param == NULL)
      setParams(source->params);
    setDisjunction(source->disjunction, source->ownsDisjunction);
    this->canReplaceGivenCuts = source->canReplaceGivenCuts;
    this->isSetupForRepeatedUse = source->isSetupForRepeatedUse;
    this->mode = source->mode;
    this->exitReason = source->exitReason;
    this->timer = source->timer;
    this->prlp = source->prlp;
    this->cutType = source->cutType;
    this->objType = source->objType;
    this->numCutsOfType = source->numCutsOfType;
    this->numCutsFromHeur = source->numCutsFromHeur;
    this->numObjFromHeur = source->numObjFromHeur;
    this->numFailsFromHeur = source->numFailsFromHeur;
    this->numFails = source->numFails;
    this->init_num_cuts = source->init_num_cuts;
    this->num_rounds = source->num_rounds;
    this->num_cuts = source->num_cuts;
    this->num_obj_tried = source->num_obj_tried;
    this->num_failures = source->num_failures;
    this->probData = source->probData;
    this->prlpData = source->prlpData;
    this->user_objectives = source->user_objectives;
    this->user_tight_points = source->user_tight_points;
  }
  else {
    this->mode = VPCMode::PARTIAL_BB;
    if (prlp) delete prlp;
    this->prlp = NULL;
    this->ownsDisjunction = false;
    this->isSetupForRepeatedUse = false;
    this->disjunction = NULL;
    for (int t = 0; t < static_cast<int>(VPCTimeStats::NUM_TIME_STATS); t++) {
      timer.register_name(VPCTimeStatsName[t]);
    }
    for (int t = 0; t < static_cast<int>(ObjectiveType::NUM_OBJECTIVE_TYPES); t++) {
      timer.register_name(ObjectiveTypeName[t] + "_TIME");
    }
    setupAsNew();
  }
} /* initialize */

/**
 * Get problem data such as min/max coeff, problem-specific epsilon,
 * nonbasic variables, row in which each variable is basic, etc.
 */
void CglVPC::getProblemData(
    /// [in,out] Solver being used to determine the nonbasic space; note that the basis and/or solution may change due to enableFactorization
    OsiSolverInterface* const solver,
    /// [out] Where to save the data
    ProblemData& probData, 
    /// [in] If this is a subproblem, then we may want to pass the original problem data (to save the locations of the original nonbasic variables)
    const ProblemData* const origProbData,
    /// [in] Whether to enable factorization (can change solution slightly)
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

  // May also need to save rays of C1, where the coefficients of inv(B) * A
  // are sometimes negated because we have
  // \bar x = inv(B) * b - inv(B) * A * x_N
  const int numNB = probData.NBVarIndex.size();
  int tempIndex = 0;
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
            "However, for nb var %d (real var %d), it is %e. Small enough error that we only send warning.\n",
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

CglVPC::ExitReason CglVPC::setupConstraints(OsiSolverInterface* const vpcsolver, OsiCuts& cuts) {
  /***********************************************************************************
   * Change initial bounds
   ***********************************************************************************/
  const int num_changed_bounds = this->disjunction->common_changed_var.size();
  int num_fixed = 0;
  for (int i = 0; i < num_changed_bounds; i++) {
    if (reachedCutLimit()) {
      printf("CglVPC::setupConstraints: Reached cut limit from one-sided cuts due to bounds fixed at the root.\n");
      return CglVPC::ExitReason::CUT_LIMIT_EXIT;
    }
    const int col = this->disjunction->common_changed_var[i];
    if (this->disjunction->common_changed_bound[i] <= 0) {
      vpcsolver->setColLower(col, this->disjunction->common_changed_value[i]);
    } else {
      vpcsolver->setColUpper(col, this->disjunction->common_changed_value[i]);
    }

    const double mult = (this->disjunction->common_changed_bound[i] <= 0) ? 1. : -1.;
    const double el = mult * 1.;
    const double val = this->disjunction->common_changed_value[i];
    OsiRowCut currCut;
    currCut.setLb(mult * val);
    currCut.setRow(1, &col, &el, false);
    // TODO set effectiveness?
    addCut(currCut, cuts, CutType::ONE_SIDED_CUT, ObjectiveType::ONE_SIDED, vpcsolver, true);

    if (isVal(vpcsolver->getColLower()[col], vpcsolver->getColUpper()[col])) {
      num_fixed++;
    }
  } // account for bounds changed at root

  const int num_added_ineqs = this->disjunction->common_ineqs.size();
  for (int i = 0; i < num_added_ineqs; i++) {
    if (reachedCutLimit()) {
      printf("CglVPC::setupConstraints: Reached cut limit from one-sided cuts added at the root.\n");
      return CglVPC::ExitReason::CUT_LIMIT_EXIT;
    }
    OsiRowCut* currCut = &(this->disjunction->common_ineqs[i]);
    vpcsolver->applyRowCuts(1, currCut); // hopefully this works
    addCut(*currCut, cuts, CutType::ONE_SIDED_CUT, ObjectiveType::ONE_SIDED, vpcsolver, true);
  }

#ifdef TRACE
//  printf(
//      "\n## Total number changed bounds: %d. Number fixed: %d. Number added inequalities: %d. ##\n",
//      num_changed_bounds, num_fixed, num_added_ineqs);
#endif
  if (num_changed_bounds + num_fixed + num_added_ineqs > 0) {
    vpcsolver->resolve();
    if (!checkSolverOptimality(vpcsolver, false)) {
      error_msg(errorstring,
          "CglVPC::setupConstraints: Solver not proven optimal after updating bounds from root node.\n");
      writeErrorToLog(errorstring, params.logfile);
      exit(1);
    }

    // For debugging purposes, verify that the stored root objective matches
    if (num_added_ineqs == 0) {
      const double newval = vpcsolver->getObjValue();
      const double savedval = this->disjunction->root_obj;
      const double epsilon = params.get(doubleConst::DIFFEPS);
      if (!isVal(newval, savedval, epsilon)) {
        vpcsolver->resolve();
      }
      if (!isVal(newval, savedval, epsilon)) {
        double ratio = 1.;
        if (isZero(newval, epsilon) && isZero(savedval, epsilon)) {
          // nothing to do, keep ratio = 1.
          ratio = 1.;
        }
        else if (isZero(newval, epsilon) || isZero(savedval, epsilon)) {
          // ratio is 1 + abs(diff between values, since one of these values is zero)
          ratio = 1. + std::abs(newval - savedval);
        }
        else {
          ratio = newval / savedval;
          if (ratio < 1.) {
            ratio = 1. / ratio;
          }
        }

        // Allow it to be up to 3% off without causing an error
        if (greaterThanVal(ratio, 1.03)) {
          error_msg(errorstring,
              "CglVPC::setupConstraints: Root objective calculated after bound changes (%f) does not match stored value in disjunction (%f), with a difference of %e.\n",
            vpcsolver->getObjValue(), this->disjunction->root_obj, std::abs(this->disjunction->root_obj - vpcsolver->getObjValue()));
          writeErrorToLog(errorstring, params.logfile);
          exit(1);
        } else {
          warning_msg(warnstring,
              "CglVPC::setupConstraints: Root objective calculated after bound changes (%f) does not match stored value in disjunction (%f), with a difference of %e.\n",
            vpcsolver->getObjValue(), this->disjunction->root_obj, std::abs(this->disjunction->root_obj - vpcsolver->getObjValue()));
        }
      } // check root objective value matches
    } // check if any inequalities have been added to the root node through the disjunction
  } // compute root LP basis after changing bounds at root and adding globally valid inequalities

  // Save problem data for nonbasic space usage and decide on problem-specific epsilon
#ifdef TRACE
  printf("\n## CglVPC::setupConstraints: Saving root node optimal solution for nonbasic space usage ##\n");
#endif
  getProblemData(vpcsolver, this->probData);
#ifdef TRACE
  printf("CglVPC::setupConstraints: Problem-specific epsilon set to %.10e\n", this->probData.EPS);
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
   * Get bases and generate VPCs
   ***********************************************************************************/
  //const int num_disj_terms = this->disjunction->num_terms; // may be different from terms.size() due to integer-feasible terms?
  const int dim = vpcsolver->getNumCols(); // NB: treating fixed vars as at bound

  int terms_added = -1;

  // Start with the integer-feasible solution
  if (!(this->disjunction->integer_sol.empty())) {
    terms_added++;

    // Get solution and calculate slacks
    const double* sol = this->disjunction->integer_sol.data();
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
    printf("\n## CglVPC::setupConstraints: Saving integer feasible solution as term %d. ##\n", terms_added);
#endif
    //term->is_feasible = true;

    // Update disjunctive bound info
    // This is here rather than in the Disjunction class,
    // because it is unclear whether, in that class,
    // the user computes with the variables changed at the root
    double objOffset = 0.;
    vpcsolver->getDblParam(OsiDblParam::OsiObjOffset, objOffset);
    const double objVal =
        dotProduct(vpcsolver->getObjCoefficients(), sol, dim) - objOffset;
    this->disjunction->updateObjValue(objVal);
//    this->disjunction->updateNBObjValue(curr_nb_obj_val);
  } // integer-feasible solution

  // Now we handle the normal terms
  const int num_normal_terms = this->disjunction->terms.size();
  for (int tmp_ind = 0; tmp_ind < num_normal_terms; tmp_ind++) {
    DisjunctiveTerm* term = &(this->disjunction->terms[tmp_ind]);
    terms_added++;

    if (!term->is_feasible && !params.get(intParam::RECYCLED_DISJUNCTION)) {
      continue;
    }

    SolverInterface* termSolver = dynamic_cast<SolverInterface*>(vpcsolver->clone());
    termSolver->disableFactorization();

    // Change bounds in the solver
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
        termSolver->setColLower(col, val);
      } else {
        termSolver->setColUpper(col, val);
      }
    }

    const int curr_num_added_ineqs = term->ineqs.size();
    for (int i = 0; i < curr_num_added_ineqs; i++) {
      OsiRowCut* currCut = &term->ineqs[i];
      termSolver->applyRowCuts(1, currCut); // hopefully this works
    }

    // Set the warm start
    if (term->basis && !(termSolver->setWarmStart(term->basis))) {
      error_msg(errorstring,
          "CglVPC::setupConstraints: Warm start information not accepted for term %d/%d.\n",
          tmp_ind + 1, num_normal_terms);
      writeErrorToLog(errorstring, params.logfile);
      exit(1);
    }

    // Resolve and check the objective matches
    // taking this off trace only because I need it for tracking number of terms in current release
    printf("\n## CglVPC::setupConstraints: Solving for term %d/%d. ##\n",
        tmp_ind + 1, num_normal_terms);
    termSolver->resolve();
    term->is_feasible = checkSolverOptimality(termSolver, true);
    //enableFactorization(termSolver, params.get(doubleParam::EPS)); // this may change the solution slightly

    // Check objective value against what is stored (if it is finite)
    // this value isn't going to make sense if we recycle the disjunction
    if (term->is_feasible && !isInfinity(term->obj) && !params.get(intParam::RECYCLED_DISJUNCTION)) {
      // Sometimes we run into a few issues getting the ``right'' value
      double newval = termSolver->getObjValue();
      const double savedval = term->obj;
      const double epsilon = params.get(doubleConst::DIFFEPS);
      // double check to make sure CLP didn't trip on its own shoelaces
      if (!isVal(newval, savedval, epsilon)) {
        termSolver->resolve();
        newval = termSolver->getObjValue();
      }
      if (!isVal(newval, savedval, epsilon)) {
        double ratio = 1.;
        if (isZero(newval, epsilon) && isZero(savedval, epsilon)) {
          // nothing to do, keep ratio = 1.
          ratio = 1.;
        }
        else if (isZero(newval, epsilon) || isZero(savedval, epsilon)) {
          // ratio is 1 + abs(diff between values, since one of these values is zero)
          ratio = 1. + std::abs(newval - savedval);
        }
        else {
          ratio = newval / savedval;
          if (ratio < 1.) {
            ratio = 1. / ratio;
          }
        }

        // Allow it to be up to 3% off without causing an error
        if (greaterThanVal(ratio, 1.03)) {
          error_msg(errorstring,
              "CglVPC::setupConstraints: Objective at disjunctive term %d/%d is incorrect. Before, it was %s, now it is %s.\n",
              tmp_ind, num_normal_terms, stringValue(term->obj, "%1.3f").c_str(),
              stringValue(termSolver->getObjValue(), "%1.3f").c_str());
          writeErrorToLog(errorstring, params.logfile);
          exit(1);
        } else {
          warning_msg(warnstring,
              "CglVPC::setupConstraints: Objective at disjunctive term %d/%d is incorrect. Before, it was %s, now it is %s.\n",
              tmp_ind, num_normal_terms, stringValue(term->obj, "%1.3f").c_str(),
              stringValue(termSolver->getObjValue(), "%1.3f").c_str());
        }
#ifdef TRACE
        std::string commonName;
        Disjunction::setCgsName(commonName, curr_num_changed_bounds, commonTermIndices,
            commonTermCoeff, commonTermRHS, false);
        printf("CglVPC::setupConstraints: Bounds changed: %s.\n", commonName.c_str());
#endif
      } // check that objective value matches
    }

    // Possibly strengthen the term
    if (term->is_feasible && params.get(intParam::STRENGTHEN) == 2) {
    // Add Gomory cuts on disjunctive term and resolve
      OsiCuts GMICs;
      CglGMI GMIGen;
      GMIGen.generateCuts(*termSolver, GMICs);
      termSolver->applyCuts(GMICs);
      termSolver->resolve();
      term->is_feasible = checkSolverOptimality(termSolver, true);
    }
    if (term->is_feasible) {
      // Update disjunctive bound info
      // This is here rather than in the Disjunction class,
      // because it is unclear whether, in that class,
      // the user computes with the variables changed at the root
      this->disjunction->updateObjValue(termSolver->getObjValue());

      timer.register_name(CglVPC::time_T1 + std::to_string(terms_added));
      timer.register_name(CglVPC::time_T2 + std::to_string(terms_added));

      if (!reachedTimeLimit(CglVPC::time_T1 + "TOTAL", params.get(TIMELIMIT))) {
        timer.start_timer(CglVPC::time_T1 + "TOTAL");
        timer.start_timer(CglVPC::time_T1 + std::to_string(terms_added));

        // Get cobasis information and PR collection
        enableFactorization(termSolver, probData.EPS);
        ProblemData tmpData;
        getProblemData(termSolver, tmpData, &probData, false);
        genDepth1PRCollection(vpcsolver, termSolver,
            probData, tmpData, terms_added);

        timer.end_timer(CglVPC::time_T1 + "TOTAL");
        timer.end_timer(CglVPC::time_T1 + std::to_string(terms_added));
      }
    } // compute PR collection if the disj term is feasible
    if (termSolver) {
      delete termSolver;
    }
  } // end loop over nodes on tree

//  {
//    // We already checked that the objective cut is valid 
//    // Now we can check that the tilted objective cut is valid
//    // We need the variable branched at the root node, the value of the variable at the root, and a bound for each side
//    // Let boundD be the "down" bound, and boundU be the "up" bound
//    // Suppose boundD <= boundU (recall this is a minimization problem)
//    // Suppose also that x_k \in [\ell_k, u_k]
//    // To ensure validity, the tilted objective cut is:
//    //   c \dot x >= boundD + (boundU - boundD) * (x_k - \floor{\bar{x}_k}) / (u_k - \floor{\bar{x}_k})
//    // A similar cut can be formulated when boundU < boundD:
//    //   c \dot x >= boundU + (boundD - boundU) * (\ceil{\bar{x}_k} - \ell_k) / (\ceil{\bar{x}_k} - \ell_k)
//    // Both cuts can be formulated as
//    //   c \dot x >= smallerBound + (biggerBound - smallerBound) * (x_k - roundedVal) / (globalBound - roundedVal)
//    int var = -1;
//    double val = 0.;
//    double boundD = std::numeric_limits<double>::lowest();
//    double boundU = std::numeric_limits<double>::lowest();
//    try {
//      RootTerm root = dynamic_cast<PartialBBDisjunction*>(this->disjunction)->root;
//      var = root.var;
//      val = root.val;
//      boundD = root.boundD;
//      boundU = root.boundU;
//    } catch (std::exception& e) {
//    }
//    if (var >= 0) {
//      if (!isInfinity(std::abs(boundD)) && !isInfinity(std::abs(boundU))) {
//        // We get x_k = \bar{x}_k - \sum_{j \in \NB} \bar{a}_{kj} x_j, 
//        // where \bar{a}_{kj} is coefficient (k,j) after multiplying A by the basis inverse
//        // The original objective (c \dot x) in the nonbasic space is read from probData.NBReducedCost
//        std::vector<double> basisRow(vpcsolver->getNumCols());
//        enableFactorization(vpcsolver, probData.EPS);
//        vpcsolver->getBInvARow(var, &(basisRow[0]));
//        // Set up things based on which bound is smaller
//        const double smallerBound = (boundD <= boundU) ? boundD : boundU;
//        const double higherBound = (boundD <= boundU) ? boundU : boundD;
//        const double globalBound = (boundD <= boundU) ? vpcsolver->getColUpper()[var] : vpcsolver->getColLower()[var];
//        const double roundedVal = (boundD <= boundU) ? std::floor(val) : std::ceil(val);
//        const double Delta = higherBound - smallerBound;
//        const double coeff = Delta / (globalBound - roundedVal);
//        const double rhs = smallerBound - Delta * roundedVal / (globalBound - roundedVal) - probData.lp_opt + coeff * val;
//
//
//      } else if (isInfinity(std::abs(boundD))) {
//      } else if (isInfinity(std::abs(boundU))) {
//      } else {
//        // Both sides infeasible, which means the problem is infeasible; error
//        error_msg(errorstring, "Both sides of root split variable (index %d, value %g) are infeasible, which would imply the instance is infeasible.\n", var, val);
//        writeErrorToLog(errorstring, params.logfile);
//        exit(1);
//      }
//    }
//  } // check tilted objective cut validity

#ifdef TRACE
  printf(
      "\nCglVPC::setupConstraints: Finished setting up constraints. Min obj val: Structural: %1.3e, NB: %1.3e.\n",
      this->disjunction->best_obj, this->disjunction->best_obj - this->probData.lp_opt);
#endif
  return CglVPC::ExitReason::SUCCESS_EXIT;
} /* setupConstraints */

/**
 * IN NON-BASIC SPACE
 * Get Point and Rays from corner polyhedron defined by current optimum at solver
 * Assumed to be optimal already
 */
void CglVPC::genDepth1PRCollection(const OsiSolverInterface* const vpcsolver,
    const OsiSolverInterface* const tmpSolver, const ProblemData& origProbData,
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
  const double beta = params.get(PRLP_FLIP_BETA) >= 0 ? 1.0 : -1.0;
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
            "Term %d (point) does not satisfy the objective cut. Activity: %.8f. RHS: %.8f. Small enough violation that we only send a warning.\n",
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
    if (isNonBasicFixedVar(tmpSolver, NBVar)) {
      continue; // fixed variable means that we can work with the lower-dimensional cone
    }

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
            "Term %d ray %d does not satisfy the objective cut. Activity: %.8f. RHS: %.8f. Small enough violation that we only send a warning.\n",
            term_ind, ray_ind, calcRedCost, 0.);
      }
    }
    prlpData.addConstraint(currRayNB, 0.0, term_ind, curr_nb_obj_val);
    // double objDepth = curr_nb_obj_val / currRayNB.twoNorm();
  } // iterate over nonbasic vars
} /* genDepth1PRCollection */

CglVPC::ExitReason CglVPC::tryObjectives(OsiCuts& cuts,
    const OsiSolverInterface* const origSolver, const OsiCuts* const structSICs) {
  if (reachedTimeLimit(VPCTimeStats::TOTAL_TIME, params.get(TIMELIMIT))) {
    return CglVPC::ExitReason::TIME_LIMIT_EXIT;
  }
  if (reachedCutLimit()) {
    return CglVPC::ExitReason::CUT_LIMIT_EXIT;
  }

  CglVPC::ExitReason status = CglVPC::ExitReason::UNKNOWN;

  // We can scale the rhs for points by min_nb_obj_val
  const double min_nb_obj_val = this->disjunction->best_obj - this->probData.lp_opt;
  const bool useScale = true && !isInfinity(std::abs(min_nb_obj_val));
  const double scale = (!useScale || lessThanVal(min_nb_obj_val, 1.)) ? 1. : min_nb_obj_val;
  double beta = params.get(intParam::PRLP_FLIP_BETA) >= 0 ? scale : -1. * scale;

  // Check that we can actually use this disjunction
  // We reject it if the best term and worst term have the same obj value as the LP opt
  const bool LP_OPT_IS_NOT_CUT = !greaterThanVal(min_nb_obj_val, 0.0);
  const bool DLB_EQUALS_DUB = !greaterThanVal(this->disjunction->worst_obj, this->disjunction->best_obj);
  const bool flipBeta = params.get(PRLP_FLIP_BETA) > 0; // if the PRLP is infeasible, we will try flipping the beta
  if (LP_OPT_IS_NOT_CUT) {
    this->numFails[static_cast<int>(FailureType::DLB_EQUALS_LPOPT_NO_OBJ)]++;
  }
  if (DLB_EQUALS_DUB) {
    this->numFails[static_cast<int>(FailureType::DLB_EQUALS_DUB_NO_OBJ)]++;
  }
  if (LP_OPT_IS_NOT_CUT && DLB_EQUALS_DUB && !flipBeta) {
    return CglVPC::ExitReason::NO_CUTS_LIKELY_EXIT;
  }

  printf("\n## CglVPC: Finished setting up constraints. Trying objectives. ##\n");
  const int init_num_obj = this->num_obj_tried;
  const int init_num_cuts = this->num_cuts;
  const int init_num_failures = this->num_failures;

  // Check conditions for using PRLP: 
  // - the lp optimum is cut away
  // - the disjunctive lower and upper bounds do not coincide
  if (!LP_OPT_IS_NOT_CUT || !DLB_EQUALS_DUB) {
    prlp = new PRLP(this);
    setLPSolverParameters(prlp, params.get(intParam::VERBOSITY), params.get(doubleParam::PRLP_TIMELIMIT));
    status = prlp->setup(scale);
  //  printf("# rows: %d\t # cols: %d\n", prlp->getNumRows(), prlp->getNumCols());
  //  printf("# points: %d\t # rays: %d\n", prlp->numPoints, prlp->numRays);

    // Try user objectives
    if (status == CglVPC::ExitReason::SUCCESS_EXIT) {
      for (const auto& v : this->user_objectives) {
        prlp->solve(cuts, v, origSolver, beta, structSICs);
        if (reachedTimeLimit(VPCTimeStats::TOTAL_TIME, params.get(TIMELIMIT))) {
          status = CglVPC::ExitReason::TIME_LIMIT_EXIT;
          break;
        } else if (reachedCutLimit()) {
          status = CglVPC::ExitReason::CUT_LIMIT_EXIT;
          break;
        } else if (reachedFailureLimit(num_cuts - init_num_cuts, num_failures - init_num_failures)) {
          status = CglVPC::ExitReason::FAIL_LIMIT_EXIT;
          break;
        } else {
          status = CglVPC::ExitReason::SUCCESS_EXIT;
        }
      }
    } // try user objectives

    // Try user tight points
    if (status == CglVPC::ExitReason::SUCCESS_EXIT) {
      // TODO
      /*for (const auto& v : this->user_tight_points) {
        //prlp->solve(cuts, v, origSolver, beta, structSICs);
        if (reachedTimeLimit(VPCTimeStats::TOTAL_TIME, params.get(TIMELIMIT))) {
          status = CglVPC::ExitReason::TIME_LIMIT_EXIT;
          break;
        } else if (reachedCutLimit()) {
          status = CglVPC::ExitReason::CUT_LIMIT_EXIT;
          break;
        } else if (reachedFailureLimit(num_cuts - init_num_cuts, num_failures - init_num_failures)) {
          status = CglVPC::ExitReason::FAIL_LIMIT_EXIT;
          break;
        } else {
          status = CglVPC::ExitReason::SUCCESS_EXIT;
        }
      }*/
    } // try user tight points

    // Try cut generation targeting points and rays
    if (status == CglVPC::ExitReason::SUCCESS_EXIT) {
      timer.start_timer(VPCTimeStatsName[static_cast<int>(VPCTimeStats::GEN_CUTS_TIME)]);
      prlp->targetStrongAndDifferentCuts(beta, cuts, origSolver, structSICs,
          VPCTimeStatsName[static_cast<int>(VPCTimeStats::TOTAL_TIME)]);
      timer.end_timer(VPCTimeStatsName[static_cast<int>(VPCTimeStats::GEN_CUTS_TIME)]);
      if (reachedTimeLimit(VPCTimeStats::TOTAL_TIME, params.get(TIMELIMIT))) {
        status = CglVPC::ExitReason::TIME_LIMIT_EXIT;
      } else if (reachedCutLimit()) {
        status = CglVPC::ExitReason::CUT_LIMIT_EXIT;
      } else if (reachedFailureLimit(num_cuts - init_num_cuts, num_failures - init_num_failures)) {
        status = CglVPC::ExitReason::FAIL_LIMIT_EXIT;
      } else {
        status = CglVPC::ExitReason::SUCCESS_EXIT;
      }
    }
  } // check conditions for using PRLP

  if (flipBeta) {
    error_msg(errorstring,
        "Currently, flipping beta is not implemented. Need to check that within setup, we can correctly switch the rhs.\n");
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }

#ifdef TRACE
  printf("\n## CglVPC: Finished trying %d objectives. Generated %d cuts. Total num cuts: %d. ##\n",
      this->num_obj_tried - init_num_obj, this->num_cuts - init_num_cuts, cuts.sizeCuts());
#endif

  if ((this->num_obj_tried - init_num_obj)
      != (this->num_cuts + this->num_failures)
          - (init_num_cuts + init_num_failures)) {
    error_msg(errorstring,
        "num_obj_tried (%d) not equal to num_cuts (%d) + num_failures (%d)\n",
        num_obj_tried - init_num_obj, num_cuts - init_num_cuts, num_failures - init_num_failures);
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }

  return status;
} /* tryObjectives */

/**
 * @details There are four types of checks. The first three have non-decreasing success requirements.
 * 1. Few cuts have been generated:
 *      default is FEW_CUTS = 1 cut, and the success threshold is at least 1 cut every 20 obj (fail ratio = .95).
 * 2. Many cuts have been generated:
 *      default is MANY_CUTS = .25 * CUT_LIMIT, and the success threshold is at least 1 cut every 10 obj (fail ratio = .90).
 * 3. Many obj have been tried, i.e., we have been trying for a long time so we better be successful super often:
 *      MANY_OBJ = max(FEW_CUTS / (1-few_cuts_fail_threshold), MANY_CUTS / (1-many_cuts_fail_threshold));
 *      default = max(20, 2.5 * CUT_LIMIT) and the success threshold is at least 1 cut every 5 obj (fail ratio = .80).
 * 4. Time is too long and we are not too successful:
 *       \# obj tried >= MANY_OBJ && time >= 10 && average time / obj >= CUTSOLVER_TIMELIMIT + 1
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
bool CglVPC::reachedFailureLimit(
    /// Number of cuts currently generated
    const int num_cuts, 
    /// Number of objectives tried that did not lead to cuts
    const int num_fails,
    /// Active if FEW_CUTS (= 1 by default) generated; by default at least 1 cut every 20 obj (fail ratio = .95)
    const double few_cuts_fail_threshold,
    /// Active if MANY_CUTS (.25 * CUT_LIMIT by default) generated; by default at least 1 cut every 10 obj (fail ratio = .90)
    const double many_cuts_fail_threshold,
    /// Active if many objectives have been tried (MANY_OBJ = max(FEW_CUTS / (1 - \p few_cuts_fail_threshold), MANY_CUTS / (1 - \p many_cuts_fail_threshold)), and by default corresponds to at least 1 cut every 5 obj (fail ratio = .80)
    const double many_obj_fail_threshold,
    /// CURRENTLY NOT USED: Time is long enough and we are not successful
    const double time_fail_threshold) const {
  const int num_obj_tried = num_cuts + num_fails;
  if (num_obj_tried == 0) {
    return false;
  }
  const int CUT_LIMIT = getCutLimit();
  const int FEW_CUTS = 1;
  const int NO_CUTS_OBJ_LIMIT = std::ceil(FEW_CUTS / (1. - few_cuts_fail_threshold));
  const int MANY_CUTS = std::ceil(.25 * CUT_LIMIT);
  int MANY_OBJ = CoinMax(
      std::ceil(FEW_CUTS / (1. - few_cuts_fail_threshold)),
      std::ceil(MANY_CUTS / (1. - many_cuts_fail_threshold)));
  if (MANY_OBJ < 0) // in case of overflow
    MANY_OBJ = std::numeric_limits<int>::max();
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
  if (reached_limit) {
//    this->exitReason = CglVPC::ExitReason::FAIL_LIMIT_EXIT;
    printf("CglVPC: Reached failure limit with %d cuts and %d fails.\n", num_cuts, num_fails);
  }
  return reached_limit;
} /* reachedFailureLimit */

void CglVPC::finish(CglVPC::ExitReason exitReason) {
  const std::string timeName = VPCTimeStatsName[static_cast<int>(VPCTimeStats::TOTAL_TIME)];
  if (exitReason == CglVPC::ExitReason::TIME_LIMIT_EXIT) {
    printf("CglVPC: Reached %s time limit %f < current time %f.\n",
        timeName.c_str(), params.get(TIMELIMIT), timer.get_total_time(timeName));
  } else if (exitReason == CglVPC::ExitReason::CUT_LIMIT_EXIT) {
    printf("CglVPC: Reached cut limit %d.\n", getCutLimit());
  } else if (exitReason == CglVPC::ExitReason::FAIL_LIMIT_EXIT) {
    // printed when reachedFailureLimit() is called
  }
  this->exitReason = exitReason;
  this->timer.end_timer(timeName);
  printf("CglVPC: Finishing with exit reason: %s. Number cuts: %d.\n", CglVPC::ExitReasonName[static_cast<int>(exitReason)].c_str(), num_cuts);
} /* finish */
