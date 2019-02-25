// Name:     DisjunctionHelper.cpp
// Author:   A. M. Kazachkov
// Date:     2018-02-25
//-----------------------------------------------------------------------------
#include "DisjunctionHelper.hpp"

// Project files
#include "SolverHelper.hpp"
#include "TimeStats.hpp"
#include "utility.hpp"
#include "VPCParameters.hpp"

// Disjunctions
#include "PartialBBDisjunction.hpp"
#include "SplitDisjunction.hpp"

//void generateRoundOfCuts(OsiCuts& vpcs, OsiSolverInterface* solver, TimeStats& timer, VPCParameters& params) {
//  timer.start_timer(OverallTimeStats::GEN_VPC_TIME);
//  CglVPC gen(params);
//  gen.generateCuts(*solver, vpcs); // solution may change slightly due to enable factorization called in getProblemData...
//  timer.end_timer(OverallTimeStats::GEN_VPC_TIME);
//
//  printf(
//      "\n## Finished VPC generation (exit reason: %s) using partial tree with # disjunctive terms = %d, # cuts generated = %d. Now solving for new objective value. ##\n",
//      ExitReasonName[static_cast<int>(gen.exitReason)].c_str(),
//      gen.disj->num_terms, vpcs.sizeCuts());
//  double init_obj_value = solver->getObjValue();
//  solver->applyCuts(vpcs);
//  solver->resolve();
//  checkSolverOptimality(solver, false);
//  const double new_obj_value = solver->getObjValue();
//
//  printf(
//      "\n## Initial obj value: %1.6f. New obj value: %1.6f. Disj lb: %1.6f. ##\n",
//      init_obj_value, new_obj_value, gen.disj->best_obj);
//} /* generateRoundOfCuts */

ExitReason setDisjunctions(std::vector<Disjunction*>& disjVec,
    const OsiSolverInterface* const si, const VPCParameters& params,
    CglVPC::VPCMode mode) {
//  CglVPC::VPCMode mode = static_cast<CglVPC::VPCMode>(params.get(MODE));
  if (mode == CglVPC::VPCMode::PARTIAL_BB) {
    if (params.get(intParam::DISJ_TERMS) < 2) {
      return ExitReason::NO_DISJUNCTION_EXIT;
    }
    PartialBBDisjunction* disj = new PartialBBDisjunction(params);
    ExitReason status = disj->prepareDisjunction(si);
    disjVec.push_back(disj);
    return status;
//      if (status == ExitReason::PARTIAL_BB_OPTIMAL_SOLUTION_FOUND_EXIT) {
//        warning_msg(warnstr,
//            "An integer (optimal) solution was found prior while getting disjunction. " "We will generate between n and 2n cuts, restricting the value of each variable.\n");
//        const double* solution = disj->integer_sol.data();
//        for (int col = 0; col < static_cast<int>(disj->integer_sol.size());
//            col++) {
//          const double val = solution[col];
//
//          // Check which of the bounds needs to be fixed
//          for (int b = 0; b < 2; b++) {
//            if ((b == 0 && greaterThanVal(val, solver->getColLower()[col]))
//                || (b == 1 && lessThanVal(val, solver->getColUpper()[col]))) {
//              const double mult = (b == 0) ? 1. : -1.;
//              const double el = mult * 1.;
//
//              OsiRowCut currCut;
//              currCut.setLb(mult * val);
//              currCut.setRow(1, &col, &el, false);
//              addCut(currCut, CutType::OPTIMALITY_CUT, cuts);
//            }
//          }
//        } // iterate over columns and add optimality cut if needed
//      }
  } // PARTIAL_BB
  else if (mode == CglVPC::VPCMode::SPLITS) {
    if (generateSplitDisjunctions(disjVec, si, params)) {
      return ExitReason::SUCCESS_EXIT;
    }
  } else {
    error_msg(errorstring,
        "Mode that is chosen has not yet been implemented for VPC generation: %s.\n",
        CglVPC::VPCModeName[static_cast<int>(mode)].c_str());
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }
  return ExitReason::UNKNOWN;
} /* setDisjunctions */

int generateSplitDisjunctions(std::vector<Disjunction*>& disjVec, const OsiSolverInterface* const si, const VPCParameters& params) {
  OsiVectorInt fracCore = si->getFractionalIndices(params.get(doubleConst::AWAY));
  if (fracCore.size() == 0)
    return 0;

  int num_splits = 0;
  std::vector<int> fracCoreSelected;
  fracCoreSelected.reserve(fracCore.size());

  // To save changed variable bounds at root node (bound <= 0 is LB, bound = 1 is UB)
  std::vector<int> common_changed_var;
  std::vector<int> common_changed_bound;
  std::vector<double> common_changed_value;

  // Set up solver for hot start
  SolverInterface* solver;
  try {
    solver = dynamic_cast<SolverInterface*>(si->clone());
  } catch (std::exception& e) {
    error_msg(errorstring,
        "Unable to clone solver into desired SolverInterface.\n");
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }
#ifdef USE_CLP
  try {
    setupClpForStrongBranching(dynamic_cast<OsiClpSolverInterface*>(solver));
  } catch (std::exception& e) {
    // It's okay, we can continue
  }
#endif
  solver->enableFactorization();
  solver->markHotStart();

  for (int var : fracCore) {
    if (num_splits >= params.get(DISJ_TERMS))
      break;

    const double val = solver->getColSolution()[var];
    const double floorxk = std::floor(val);
    const double ceilxk = std::ceil(val);

    if (!si->isInteger(var)) {
      error_msg(errorstring, "Chosen variable %d is not an integer variable.\n", var);
      writeErrorToLog(errorstring, params.logfile);
      exit(1);
    }
    if (isVal(val, floorxk, params.get(doubleConst::AWAY))
          || isVal(val, ceilxk, params.get(doubleConst::AWAY))) {
      error_msg(errorstring, "Chosen variable %d is not fractional (value: %1.6e).\n", var, val);
      writeErrorToLog(errorstring, params.logfile);
      exit(1);
    }

    const double origLB = solver->getColLower()[var];
    const double origUB = solver->getColUpper()[var];
    bool downBranchFeasible = true, upBranchFeasible = true;

    // Check down branch
    solver->setColUpper(var, floorxk);
    solveFromHotStart(solver, var, true, origUB, floorxk);
    if (solver->isProvenOptimal()) {
//      addTerm(var, 1, floorxk, solver);
    } else if (solver->isProvenPrimalInfeasible()) {
      downBranchFeasible = false;
    } else {
      // Something strange happened
      error_msg(errorstring,
          "Down branch is neither optimal nor primal infeasible on variable %d (value %e).\n",
          var, val);
      writeErrorToLog(errorstring, params.logfile);
      exit(1);
    }
    solver->setColUpper(var, origUB);

    // Return to previous state
    solver->solveFromHotStart();

    // Check up branch
    solver->setColLower(var, ceilxk);
    solveFromHotStart(solver, var, false, origLB, ceilxk);
    if (solver->isProvenOptimal()) {
//      addTerm(var, 0, ceilxk, solver);
    } else if (solver->isProvenPrimalInfeasible()) {
      upBranchFeasible = false;
    } else {
      // Something strange happened
      error_msg(errorstring,
          "Up branch is neither optimal nor primal infeasible on variable %d (value %e).\n",
          var, val);
      writeErrorToLog(errorstring, params.logfile);
      exit(1);
    }
    solver->setColLower(var, origLB);

    // Return to original state
    solver->solveFromHotStart();

    // Check if some side of the split is infeasible
    if (!downBranchFeasible || !upBranchFeasible) {
      if (!downBranchFeasible && !upBranchFeasible) {
        // Infeasible problem
        error_msg(errorstring,
            "Infeasible problem due to integer variable %d (value %e).\n",
            var, val);
        writeErrorToLog(errorstring, params.logfile);
        exit(1);
      }
      // If one side infeasible, can delete last term added and fix variables instead
      common_changed_var.push_back(var);
      if (!downBranchFeasible) { // set lb to ceilxk
        common_changed_bound.push_back(0);
        common_changed_value.push_back(ceilxk);
      }
      if (!upBranchFeasible) { // set ub to floorxk
        common_changed_bound.push_back(1);
        common_changed_value.push_back(floorxk);
      }
    } // check infeasibility
    else {
      num_splits++;
      fracCoreSelected.push_back(var);
    }
  } // loop through fractional core
  solver->unmarkHotStart();
  solver->disableFactorization();
  if (solver)
    delete solver;

  for (int var : fracCoreSelected) {
    SplitDisjunction* disj = new SplitDisjunction;
    disj->var = var;
    disj->prepareDisjunction(si);
    disjVec.push_back(disj);
  }
  { // DEBUG
    for (int i = 0; i < num_splits; ++i) {
      printf("Var: %d\n", dynamic_cast<SplitDisjunction*>(disjVec[i])->var);
    }
  }
  return num_splits;
} /* generateCutsFromSplitDisjunctions */
