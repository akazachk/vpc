// Name:     preprocess.cpp
// Author:   A. M. Kazachkov
// Date:     2019-03-16
//-----------------------------------------------------------------------------
#include "preprocess.hpp"

// Project files
#include "BBHelper.hpp"
#include "GurobiHelper.hpp"
#include "SolverHelper.hpp"
#include "VPCParameters.hpp"

// COIN-OR
#include <OsiSolverInterface.hpp>

/**
 * Perform preprocessing and get statistics
 */
void performCleaning(const VPCParameters& orig_params,
    OsiSolverInterface* solver, std::string& filename, const double ip_obj,
    const int CLEANING_MODE_OPTION) {
  const bool DO_BRANCHING_WITH_SICS = false;
  const bool DO_STRONG_BRANCHING = false;
  const bool DO_DEFAULT = true;
  const bool DO_CUTSON = false;
  const bool DO_CUTSOFF = false;
  const int numOrigRows = solver->getNumRows();
  const int numOrigCols = solver->getNumCols();
  int numBoundsChanged = 0;
  int numSBFixed = 0;
  const double origLPOpt = solver->getObjValue();
  VPCParameters params = orig_params;

#ifdef TRACE
  printf("Collecting information about original instance.\n");
#endif

  // Get number of integer and binary variables
  // As well as primal and dual degeneracy
  int numNonZero = solver->getMatrixByRow()->getNumElements();
  int numInteger = 0, numBinary = 0;
  int origPrimalDegen = 0, origDualDegen = 0;
  for (int col = 0; col < solver->getNumCols(); col++) {
    if (solver->isInteger(col)) {
      numInteger++;
      if (solver->isBinary(col)) {
        numBinary++;
      }
    }

    const bool isBasic = isBasicCol(solver, col);
    const double val = solver->getColSolution()[col];
    const double lb = solver->getColLower()[col];
    const double ub = solver->getColUpper()[col];
    if (isBasic && (isVal(val, lb) || isVal(val, ub))) {
      origPrimalDegen++;
    }
    if (!isBasic && isZero(solver->getReducedCost()[col])) {
      origDualDegen++;
    }
  } /* get col info */
  for (int row = 0; row < solver->getNumRows(); row++) {
    const bool isBasic = isBasicSlack(solver, row);
    const double val = std::abs(solver->getRowActivity()[row] - solver->getRightHandSide()[row]);
    if (isBasic && isZero(val)) {
      origPrimalDegen++;
    }
    if (!isBasic && isZero(solver->getRowPrice()[row])) {
      origDualDegen++;
    }
  } /* get row info */

  // Save the original branch-and-bound strategy
  const int strategy = params.get(BB_STRATEGY);
  const int default_yesusercuts_strategy = strategy;
  const int cutson_yesusercuts_strategy =
      enable_bb_option(enable_bb_option(
          disable_bb_option(strategy, BB_Strategy_Options::all_cuts_off),
          BB_Strategy_Options::all_cuts_on),
          BB_Strategy_Options::user_cuts);
//      ((strategy & ~(BB_Strategy_Options::all_cuts_off)) | BB_Strategy_Options::all_cuts_on) | BB_Strategy_Options::user_cuts;
  const int cutsoff_yesusercuts_strategy =
      enable_bb_option(enable_bb_option(
          disable_bb_option(
              disable_bb_option(strategy, BB_Strategy_Options::all_cuts_on),
              BB_Strategy_Options::gmics_on),
          BB_Strategy_Options::all_cuts_off), BB_Strategy_Options::user_cuts);
//      ((strategy & ~(BB_Strategy_Options::all_cuts_on) & ~(BB_Strategy_Options::gmics_on)) | BB_Strategy_Options::all_cuts_off) | BB_Strategy_Options::user_cuts;

  // Cbc
  BBInfo origBBInfoCbcDefault, origBBInfoCbcCutsOn, origBBInfoCbcCutsOff;
  if (use_bb_option(strategy, BB_Strategy_Options::cbc)) {
    if (DO_DEFAULT) {
      // Default
  #ifdef TRACE
    printf("Performing branch-and-bound with Cbc (default) on original instance.\n");
  #endif
      params.set(BB_STRATEGY, default_yesusercuts_strategy);
      doBranchAndBoundNoCuts(params, solver, origBBInfoCbcDefault);
    }

    if (DO_CUTSON) {
      // Cuts on
#ifdef TRACE
  printf("Performing branch-and-bound with Cbc (cuts on) on original instance.\n");
#endif
      params.set(BB_STRATEGY,
          cutson_yesusercuts_strategy);
      doBranchAndBoundNoCuts(params, solver, origBBInfoCbcCutsOn);
    }

    if (DO_CUTSOFF) {
      // Cuts off
#ifdef TRACE
  printf("Performing branch-and-bound with Cbc (cuts off) on original instance.\n");
#endif
      params.set(BB_STRATEGY,
          cutsoff_yesusercuts_strategy);
      doBranchAndBoundNoCuts(params, solver, origBBInfoCbcCutsOff);
    }
  } /* cbc */

  // Cplex
  BBInfo origBBInfoCplexDefault, origBBInfoCplexCutsOn, origBBInfoCplexCutsOff;
#ifdef VPC_USE_CPLEX
  if (strategy & BB_Strategy_Options::cplex) {
    if (DO_DEFAULT) {
      // Default
  #ifdef TRACE
      printf("Performing branch-and-bound with Cplex (default) on original instance.\n");
  #endif
      params.set(BB_STRATEGY, default_yesusercuts_strategy);
      doBranchAndBoundWithCplexCallable(filename.c_str(), origBBInfoCplexDefault);
    }

    if (DO_CUTSON) {
      // Cuts on
#ifdef TRACE
      printf("Performing branch-and-bound with Cplex (cuts on) on original instance.\n");
#endif
      params.set(BB_STRATEGY, cutson_yesusercuts_strategy);
      doBranchAndBoundWithCplexCallable(filename.c_str(), origBBInfoCplexCutsOn);
    }

    if (DO_CUTSOFF) {
      // Cuts off
#ifdef TRACE
      printf("Performing branch-and-bound with Cplex (cuts off) on original instance.\n");
#endif
      params.set(BB_STRATEGY, cutsoff_yesusercuts_strategy);
      doBranchAndBoundWithCplexCallable(filename.c_str(), origBBInfoCplexCutsOff);
    }
  } /* cplex */
#endif

  // Gurobi
  BBInfo origBBInfoGurobiDefault, origBBInfoGurobiCutsOn, origBBInfoGurobiCutsOff;
#ifdef USE_GUROBI
  if (use_bb_option(strategy, BB_Strategy_Options::gurobi)) {
    if (DO_DEFAULT) {
      // Default
  #ifdef TRACE
      printf("Performing branch-and-bound with Gurobi (default) on original instance.\n");
  #endif
      doBranchAndBoundWithGurobi(params, default_yesusercuts_strategy, filename.c_str(), origBBInfoGurobiDefault);
    }

    if (DO_CUTSON) {
      // Cuts on
#ifdef TRACE
      printf("Performing branch-and-bound with Gurobi (cuts on) on original instance.\n");
#endif
      doBranchAndBoundWithGurobi(params, cutson_yesusercuts_strategy, filename.c_str(), origBBInfoGurobiCutsOn);
    }

    if (DO_CUTSOFF) {
      // Cuts off
#ifdef TRACE
      printf("Performing branch-and-bound with Gurobi (cuts off) on original instance.\n");
#endif
      doBranchAndBoundWithGurobi(params, cutsoff_yesusercuts_strategy, filename.c_str(), origBBInfoGurobiCutsOff);
    }
  } /* gurobi */
#endif

  /********** Now we do the cleaning **********/
  std::string presolved_name_stub =
      (CLEANING_MODE_OPTION <= 1) ? filename + "_presolved" : "";
  const std::string cleaned_name =
      (CLEANING_MODE_OPTION <= 1) ?
          presolved_name_stub : filename + "_cleaned.mps";

  // First get presolved opt using CPLEX / Gurobi
  double presolvedLPOptCplex = 0., presolvedLPOptGurobi = 0.;
#ifdef VPC_USE_CPLEX
  if (strategy & BB_Strategy_Options::cplex) {
#ifdef TRACE
    printf("Presolve model with Cplex.\n");
#endif
    presolveModelWithCplexCallable(
        filename.c_str(), presolvedLPOptCplex,
        presolved_name_stub); // returns mps.gz
    solver->readMps(presolved_name_stub.c_str());

    // Make sure we are doing a minimization problem; this is just to make later
    // comparisons simpler (i.e., a higher LP obj after adding the cut is better).
    if (solver->getObjSense() < param.getEPS()) {
      printf(
          "\n## Detected maximization problem. Negating objective function to make it minimization. ##\n");
      solver->setObjSense(1.0);
      const double* obj = solver->getObjCoefficients();
      for (int col = 0; col < solver->getNumCols(); col++) {
        solver->setObjCoeff(col, -1. * obj[col]);
      }
      double objOffset = 0.;
      solver->getDblParam(OsiDblParam::OsiObjOffset, objOffset);
      if (objOffset != 0.) {
        solver->setDblParam(OsiDblParam::OsiObjOffset, -1. * objOffset);
      }
    }

    // Perform initial solve
    solver->initialSolve();
    if (!checkSolverOptimality(solver, false)) {
      error_msg(errorstring, "After initial solve, solver is not optimal.\n");
      writeErrorToLog(errorstring, params.logfile);
      exit(1);
    }

    if (CLEANING_MODE_OPTION > 1) {
      remove(presolved_name_stub.c_str()); // remove the temporary file
    }
  }
#endif // use_cplex
#ifdef USE_GUROBI
  if (use_bb_option(strategy, BB_Strategy_Options::gurobi)) {
#ifdef TRACE
    printf("Presolve model with Gurobi.\n");
#endif
    presolveModelWithGurobi(params, strategy, filename.c_str(),
        presolvedLPOptGurobi, presolved_name_stub, ip_obj); // returns mps.gz
    solver->readMps(presolved_name_stub.c_str());

    // Make sure we are doing a minimization problem; this is just to make later
    // comparisons simpler (i.e., a higher LP obj after adding the cut is better).
    if (solver->getObjSense() < 1e-3) {
      printf(
          "\n## Detected maximization problem. Negating objective function to make it minimization. ##\n");
      solver->setObjSense(1.0);
      const double* obj = solver->getObjCoefficients();
      for (int col = 0; col < solver->getNumCols(); col++) {
        solver->setObjCoeff(col, -1. * obj[col]);
      }
      double objOffset = 0.;
      solver->getDblParam(OsiDblParam::OsiObjOffset, objOffset);
      if (objOffset != 0.) {
        solver->setDblParam(OsiDblParam::OsiObjOffset, -1. * objOffset);
      }
    }

    // Perform initial solve
    solver->initialSolve();
    if (!checkSolverOptimality(solver, false)) {
      error_msg(errorstring, "After initial solve, solver is not optimal.\n");
      writeErrorToLog(errorstring, params.logfile);
      exit(1);
    }

    if (CLEANING_MODE_OPTION > 1) {
      remove(presolved_name_stub.c_str()); // remove the temporary file
    }
  }
#endif // use_gurobi

  // Now clean using our own methods
  bool is_clean = (CLEANING_MODE_OPTION <= 1);
  while (!is_clean) {
#ifdef TRACE
    printf("Clean model with custom method (remove one-sided split disjunctions and tighten bounds, iteratively).\n");
#endif
    is_clean = true; //cleanProblem(solver, numBoundsChanged, numSBFixed);
  }

  // Get new solver info, including primal and dual degeneracy values
  const int cleanedNumCols = solver->getNumCols();
  const int cleanedNumRows = solver->getNumRows();
  const int cleanedNumNonZero = solver->getMatrixByRow()->getNumElements();
  int cleanedNumInteger = 0, cleanedNumBinary = 0;
  int cleanedPrimalDegen = 0, cleanedDualDegen = 0;
  int cleanedNumSICs = 0, cleanedNumSICsRd2 = 0, cleanedNumSICsStr = 0, cleanedNumSICsStrRd2 = 0;
  double cleanedSICOpt = 0., cleanedSICOptRd2 = 0.;
  double cleanedSICStrOpt = 0., cleanedSICStrOptRd2 = 0.;

  // Resolve
  if (cleanedNumNonZero > 0) {
    solver->resolve();
    if (!checkSolverOptimality(solver, false)) {
      error_msg(errorstring, "After cleaning, solver is not optimal.\n");
      writeErrorToLog(errorstring, params.logfile);
      exit(1);
    }
  }
  const double cleanedLPOpt = (cleanedNumNonZero > 0) ? solver->getObjValue() : ip_obj;

  if (CLEANING_MODE_OPTION > 1) {
    // Save cleaned LP to in directory
    solver->writeMps(cleaned_name.c_str(), "", solver->getObjSense());
  }

  // Set up new BBInfos
  BBInfo cleanedBBInfoCbcDefault, cleanedSICStrBBInfoCbcDefault;
  BBInfo cleanedBBInfoCbcCutsOn, cleanedSICStrBBInfoCbcCutsOn;
  BBInfo cleanedBBInfoCbcCutsOff, cleanedSICStrBBInfoCbcCutsOff;

  BBInfo cleanedBBInfoCplexDefault, cleanedSICStrBBInfoCplexDefault;
  BBInfo cleanedBBInfoCplexCutsOn, cleanedSICStrBBInfoCplexCutsOn;
  BBInfo cleanedBBInfoCplexCutsOff, cleanedSICStrBBInfoCplexCutsOff;

  BBInfo cleanedBBInfoGurobiDefault, cleanedSICStrBBInfoGurobiDefault;
  BBInfo cleanedBBInfoGurobiCutsOn, cleanedSICStrBBInfoGurobiCutsOn;
  BBInfo cleanedBBInfoGurobiCutsOff, cleanedSICStrBBInfoGurobiCutsOff;

  const bool was_cleaned = (numSBFixed > 0) || (numBoundsChanged > 0)
      || (cleanedNumCols != numOrigCols)
      || (cleanedNumRows != numOrigRows)
      || (cleanedNumNonZero != numNonZero);
  if (cleanedNumNonZero > 0 && was_cleaned) {
#ifdef TRACE
  printf("Collecting information about cleaned instance.\n");
#endif
    for (int col = 0; col < solver->getNumCols(); col++) {
      if (solver->isInteger(col)) {
        cleanedNumInteger++;
        if (solver->isBinary(col)) {
          cleanedNumBinary++;
        }
      }

      const bool isBasic = isBasicCol(solver, col);
      const double val = solver->getColSolution()[col];
      const double lb = solver->getColLower()[col];
      const double ub = solver->getColUpper()[col];
      if (isBasic && (isVal(val, lb) || isVal(val, ub))) {
        cleanedPrimalDegen++;
      }
      if (!isBasic && isZero(solver->getReducedCost()[col])) {
        cleanedDualDegen++;
      }
    }
    for (int row = 0; row < solver->getNumRows(); row++) {
      const bool isBasic = isBasicSlack(solver, row);
      const double val = std::abs(solver->getRowActivity()[row] - solver->getRightHandSide()[row]);
      if (isBasic && isZero(val)) {
        cleanedPrimalDegen++;
      }
      if (!isBasic && isZero(solver->getRowPrice()[row])) {
        cleanedDualDegen++;
      }
    }

    // Cbc
    if (use_bb_option(strategy, BB_Strategy_Options::cbc)) {
      if (DO_DEFAULT) {
        // Default
  #ifdef TRACE
        printf("Performing branch-and-bound with Cbc (default) on cleaned instance.\n");
  #endif
        params.set(BB_STRATEGY, default_yesusercuts_strategy);
        doBranchAndBoundNoCuts(params, solver, cleanedBBInfoCbcDefault);
      }

      if (DO_CUTSON) {
        // Cuts on
#ifdef TRACE
        printf("Performing branch-and-bound with Cbc (cuts on) on cleaned instance.\n");
#endif
        params.set(BB_STRATEGY, cutson_yesusercuts_strategy);
        doBranchAndBoundNoCuts(params, solver, cleanedBBInfoCbcCutsOn);
      }

      if (DO_CUTSOFF) {
        // Cuts off
#ifdef TRACE
        printf("Performing branch-and-bound with Cbc (cuts off) on cleaned instance.\n");
#endif
        params.set(BB_STRATEGY, cutsoff_yesusercuts_strategy);
        doBranchAndBoundNoCuts(params, solver, cleanedBBInfoCbcCutsOff);
      }
    } /* cbc */

#ifdef VPC_USE_CPLEX
    // Cplex
    if (strategy & BB_Strategy_Options::cplex) {
      if (DO_DEFAULT) {
        // Default
  #ifdef TRACE
        printf("Performing branch-and-bound with Cplex (default) on cleaned instance.\n");
  #endif
        params.set(BB_STRATEGY, default_yesusercuts_strategy);
        doBranchAndBoundWithCplexCallable(cleaned_name.c_str(), cleanedBBInfoCplexDefault);
      }

      if (DO_CUTSON) {
        // Cuts on
#ifdef TRACE
        printf("Performing branch-and-bound with Cplex (cuts on) on cleaned instance.\n");
#endif
        params.set(BB_STRATEGY, cutson_yesusercuts_strategy);
        doBranchAndBoundWithCplexCallable(cleaned_name.c_str(), cleanedBBInfoCplexCutsOn);
      }

      if (DO_CUTSOFF) {
        // Cuts off
#ifdef TRACE
        printf("Performing branch-and-bound with Cplex (cuts off) on cleaned instance.\n");
#endif
        params.set(BB_STRATEGY, cutsoff_yesusercuts_strategy);
        doBranchAndBoundWithCplexCallable(cleaned_name.c_str(), cleanedBBInfoCplexCutsOff);
      }
    } /* cplex */
#endif

#ifdef USE_GUROBI
    // Gurobi
    if (use_bb_option(strategy, BB_Strategy_Options::gurobi)) {
      if (DO_DEFAULT) {
        // Default
  #ifdef TRACE
        printf("Performing branch-and-bound with Gurobi (default) on cleaned instance.\n");
  #endif
        doBranchAndBoundWithGurobi(params, default_yesusercuts_strategy, cleaned_name.c_str(), cleanedBBInfoGurobiDefault);
      }

      if (DO_CUTSON) {
        // Cuts on
#ifdef TRACE
        printf("Performing branch-and-bound with Gurobi (cuts on) on cleaned instance.\n");
#endif
        doBranchAndBoundWithGurobi(params, cutson_yesusercuts_strategy, cleaned_name.c_str(), cleanedBBInfoGurobiCutsOn);
      }

      if (DO_CUTSOFF) {
        // Cuts off
#ifdef TRACE
        printf("Performing branch-and-bound with Gurobi (cuts off) on cleaned instance.\n");
#endif
        doBranchAndBoundWithGurobi(params, cutsoff_yesusercuts_strategy, cleaned_name.c_str(), cleanedBBInfoGurobiCutsOff);
      }
    } /* gurobi */
#endif // use_gurobi
  } /* check if any cleaning was performed */
  else if (cleanedNumNonZero == 0) {
    presolvedLPOptCplex = ip_obj;
    presolvedLPOptGurobi = ip_obj;
    cleanedSICOpt = ip_obj;
    cleanedSICOptRd2 = ip_obj;
    cleanedSICStrOpt = ip_obj;
    cleanedSICStrOptRd2 = ip_obj;

    // Cbc
    initializeBBInfo(cleanedBBInfoCbcDefault, ip_obj);
    initializeBBInfo(cleanedBBInfoCbcCutsOn, ip_obj);
    initializeBBInfo(cleanedBBInfoCbcCutsOff, ip_obj);
    initializeBBInfo(cleanedSICStrBBInfoCbcDefault, ip_obj);
    initializeBBInfo(cleanedSICStrBBInfoCbcCutsOn, ip_obj);
    initializeBBInfo(cleanedSICStrBBInfoCbcCutsOff, ip_obj);

    // Cplex
    initializeBBInfo(cleanedBBInfoCplexDefault, ip_obj);
    initializeBBInfo(cleanedBBInfoCplexCutsOn, ip_obj);
    initializeBBInfo(cleanedBBInfoCplexCutsOff, ip_obj);
    initializeBBInfo(cleanedSICStrBBInfoCplexDefault, ip_obj);
    initializeBBInfo(cleanedSICStrBBInfoCplexCutsOn, ip_obj);
    initializeBBInfo(cleanedSICStrBBInfoCplexCutsOff, ip_obj);

    // Gurobi
    initializeBBInfo(cleanedBBInfoGurobiDefault, ip_obj);
    initializeBBInfo(cleanedBBInfoGurobiCutsOn, ip_obj);
    initializeBBInfo(cleanedBBInfoGurobiCutsOff, ip_obj);
  } /* check if cleaning yielded an empty problem */
  else {
    cleanedNumInteger = numInteger;
    cleanedNumBinary = numBinary;
    cleanedPrimalDegen = origPrimalDegen;
    cleanedDualDegen = origDualDegen;

    // Cbc
    cleanedBBInfoCbcDefault = origBBInfoCbcDefault;
    cleanedBBInfoCbcCutsOn = origBBInfoCbcCutsOn;
    cleanedBBInfoCbcCutsOff= origBBInfoCbcCutsOff;

//    cleanedSICStrBBInfoCbcDefault = origSICStrBBInfoCbcDefault;
//    cleanedSICStrBBInfoCbcCutsOn = origSICStrBBInfoCbcCutsOn;
//    cleanedSICStrBBInfoCbcCutsOff = origSICStrBBInfoCbcCutsOff;

    // Cplex
    cleanedBBInfoCplexDefault = origBBInfoCplexDefault;
    cleanedBBInfoCplexCutsOn = origBBInfoCplexCutsOn;
    cleanedBBInfoCplexCutsOff = origBBInfoCplexCutsOff;
//
//    cleanedSICStrBBInfoCplexDefault = origSICStrBBInfoCplexDefault;
//    cleanedSICStrBBInfoCplexCutsOn = origSICStrBBInfoCplexCutsOn;
//    cleanedSICStrBBInfoCplexCutsOff = origSICStrBBInfoCplexCutsOff;

    // Gurobi
    cleanedBBInfoGurobiDefault = origBBInfoGurobiDefault;
    cleanedBBInfoGurobiCutsOn = origBBInfoGurobiCutsOn;
    cleanedBBInfoGurobiCutsOff = origBBInfoGurobiCutsOff;

//    cleanedSICStrBBInfoGurobiDefault = origSICStrBBInfoGurobiDefault;
//    cleanedSICStrBBInfoGurobiCutsOn = origSICStrBBInfoGurobiCutsOn;
//    cleanedSICStrBBInfoGurobiCutsOff = origSICStrBBInfoGurobiCutsOff;

//    cleanedNumSICs = structSICsOriginal.sizeCuts();
//    cleanedNumSICsRd2 = structSICsOriginalRd2.sizeCuts();
//    cleanedNumSICsStr = structSICsOriginalStr.sizeCuts();
//    cleanedNumSICsStrRd2 = structSICsOriginalRd2Str.sizeCuts();
//    cleanedSICOpt = origSICOpt;
//    cleanedSICOptRd2 = origSICOptRd2;
//    cleanedSICStrOpt = origSICStrOpt;
//    cleanedSICStrOptRd2 = origSICStrOptRd2;
  } /* strong branching did nothing so do not repeat the experiments */

  // Reset B&B strategy
  params.set(BB_STRATEGY, strategy);

//  //// Original problem
//  writeEntryToLog(numOrigRows, params.logfile);
//  writeEntryToLog(numOrigCols, params.logfile);
//  writeEntryToLog(numNonZero, params.logfile);
//  writeEntryToLog(numInteger, params.logfile);
//  writeEntryToLog(numBinary, params.logfile);
//  writeEntryToLog(origLPOpt, params.logfile);
//  writeEntryToLog(origSBLB, params.logfile);
//  writeEntryToLog(origPrimalDegen, params.logfile);
//  writeEntryToLog(origDualDegen, params.logfile);
//  writeEntryToLog(structSICsOriginal.sizeCuts(), params.logfile);
//  writeEntryToLog(origSICOpt, params.logfile);
//  writeEntryToLog(structSICsOriginalRd2.sizeCuts(), params.logfile);
//  writeEntryToLog(origSICOptRd2, params.logfile);
//  writeEntryToLog(structSICsOriginalStr.sizeCuts(), params.logfile);
//  writeEntryToLog(origSICStrOpt, params.logfile);
//  writeEntryToLog(structSICsOriginalRd2Str.sizeCuts(), params.logfile);
//  writeEntryToLog(origSICStrOptRd2, params.logfile);
//  // Cbc
//  if (param.getParamVal(ParamIndices::BB_STRATEGY_PARAM_IND) & BB_Strategy_Options::cbc) {
//    printBBInfo(origBBInfoCbcDefault, params.logfile, !DO_DEFAULT);
//    printBBInfo(origBBInfoCbcCutsOn, params.logfile, !DO_CUTSON);
//    printBBInfo(origBBInfoCbcCutsOff, params.logfile, !DO_CUTSOFF);
//    printBBInfo(origSICStrBBInfoCbcDefault, params.logfile,
//        !(DO_BRANCHING_WITH_SICS && DO_DEFAULT));
//    printBBInfo(origSICStrBBInfoCbcCutsOn, params.logfile,
//        !(DO_BRANCHING_WITH_SICS && DO_CUTSON));
//    printBBInfo(origSICStrBBInfoCbcCutsOff, params.logfile,
//        !(DO_BRANCHING_WITH_SICS && DO_CUTSOFF));
//  } /* cbc */
//  else {
//    for (int i = 0; i < (int) BB_INFO_CONTENTS.size() * 6; i++) {
//      writeEntryToLog("", params.logfile);
//    }
//  }
//  // Cplex
//  if (param.getParamVal(ParamIndices::BB_STRATEGY_PARAM_IND) & BB_Strategy_Options::cplex) {
//    printBBInfo(origBBInfoCplexDefault, params.logfile, !DO_DEFAULT);
//    printBBInfo(origBBInfoCplexCutsOn, params.logfile, !DO_CUTSON);
//    printBBInfo(origBBInfoCplexCutsOff, params.logfile, !DO_CUTSOFF);
//    printBBInfo(origSICStrBBInfoCplexDefault, params.logfile,
//        !(DO_BRANCHING_WITH_SICS && DO_DEFAULT));
//    printBBInfo(origSICStrBBInfoCplexCutsOn, params.logfile,
//        !(DO_BRANCHING_WITH_SICS && DO_CUTSON));
//    printBBInfo(origSICStrBBInfoCplexCutsOff, params.logfile,
//        !(DO_BRANCHING_WITH_SICS && DO_CUTSOFF));
//  } /* cplex */
//  else {
//    for (int i = 0; i < (int) BB_INFO_CONTENTS.size() * 6; i++) {
//      writeEntryToLog("", params.logfile);
//    }
//  }
//  // Gurobi
//  if (param.getParamVal(ParamIndices::BB_STRATEGY_PARAM_IND) & BB_Strategy_Options::gurobi) {
//    printBBInfo(origBBInfoGurobiDefault, params.logfile, !DO_DEFAULT);
//    printBBInfo(origBBInfoGurobiCutsOn, params.logfile, !DO_CUTSON);
//    printBBInfo(origBBInfoGurobiCutsOff, params.logfile,
//        !DO_CUTSOFF);
//    printBBInfo(origSICStrBBInfoGurobiDefault, params.logfile,
//        !(DO_BRANCHING_WITH_SICS && DO_DEFAULT));
//    printBBInfo(origSICStrBBInfoGurobiCutsOn, params.logfile,
//        !(DO_BRANCHING_WITH_SICS && DO_CUTSON));
//    printBBInfo(origSICStrBBInfoGurobiCutsOff, params.logfile,
//        !(DO_BRANCHING_WITH_SICS && DO_CUTSOFF));
//  } /* gurobi */
//  else {
//    for (int i = 0; i < (int) BB_INFO_CONTENTS.size() * 6; i++) {
//      writeEntryToLog("", params.logfile);
//    }
//  }
//  //// Cleaned problem
//  writeEntryToLog(solver->getNumRows(), params.logfile);
//  writeEntryToLog(solver->getNumCols(), params.logfile);
//  writeEntryToLog(cleanedNumNonZero, params.logfile);
//  writeEntryToLog(cleanedNumInteger, params.logfile);
//  writeEntryToLog(cleanedNumBinary, params.logfile);
//  writeEntryToLog(numBoundsChanged, params.logfile);
//  writeEntryToLog(numSBFixed, params.logfile);
//  if (param.getParamVal(ParamIndices::BB_STRATEGY_PARAM_IND) & BB_Strategy_Options::cplex) {
//    writeEntryToLog(presolvedLPOptCplex, params.logfile);
//  } else {
//    writeEntryToLog("", params.logfile);
//  }
//  if (param.getParamVal(ParamIndices::BB_STRATEGY_PARAM_IND) & BB_Strategy_Options::gurobi) {
//    writeEntryToLog(presolvedLPOptGurobi, params.logfile);
//  } else {
//    writeEntryToLog("", params.logfile);
//  }
//  writeEntryToLog(cleanedLPOpt, params.logfile);
//  writeEntryToLog(cleanedSBLB, params.logfile);
//  writeEntryToLog(cleanedPrimalDegen, params.logfile);
//  writeEntryToLog(cleanedDualDegen, params.logfile);
//  writeEntryToLog(cleanedNumSICs, params.logfile);
//  writeEntryToLog(cleanedSICOpt, params.logfile);
//  writeEntryToLog(cleanedNumSICsRd2, params.logfile);
//  writeEntryToLog(cleanedSICOptRd2, params.logfile);
//  writeEntryToLog(cleanedNumSICsStr, params.logfile);
//  writeEntryToLog(cleanedSICStrOpt, params.logfile);
//  writeEntryToLog(cleanedNumSICsStrRd2, params.logfile);
//  writeEntryToLog(cleanedSICStrOptRd2, params.logfile);
//  // Cbc
//  if (param.getParamVal(ParamIndices::BB_STRATEGY_PARAM_IND) & BB_Strategy_Options::cbc) {
//    printBBInfo(cleanedBBInfoCbcDefault, params.logfile, !DO_DEFAULT);
//    printBBInfo(cleanedBBInfoCbcCutsOn, params.logfile, !DO_CUTSON);
//    printBBInfo(cleanedBBInfoCbcCutsOff, params.logfile,
//        !DO_CUTSOFF);
//    printBBInfo(cleanedSICStrBBInfoCbcDefault, params.logfile,
//        !(DO_BRANCHING_WITH_SICS && DO_DEFAULT));
//    printBBInfo(cleanedSICStrBBInfoCbcCutsOn, params.logfile,
//        !(DO_BRANCHING_WITH_SICS && DO_CUTSON));
//    printBBInfo(cleanedSICStrBBInfoCbcCutsOff, params.logfile,
//        !(DO_BRANCHING_WITH_SICS && DO_CUTSOFF));
//  } /* cbc */
//  else {
//    for (int i = 0; i < (int) BB_INFO_CONTENTS.size() * 6; i++) {
//      writeEntryToLog("", params.logfile);
//    }
//  }
//  // Cplex
//  if (param.getParamVal(ParamIndices::BB_STRATEGY_PARAM_IND) & BB_Strategy_Options::cplex) {
//    printBBInfo(cleanedBBInfoCplexDefault, params.logfile,
//        !DO_DEFAULT);
//    printBBInfo(cleanedBBInfoCplexCutsOn, params.logfile,
//        !DO_CUTSON);
//    printBBInfo(cleanedBBInfoCplexCutsOff, params.logfile,
//        !DO_CUTSOFF);
//    printBBInfo(cleanedSICStrBBInfoCplexDefault, params.logfile,
//        !(DO_BRANCHING_WITH_SICS && DO_DEFAULT));
//    printBBInfo(cleanedSICStrBBInfoCplexCutsOn, params.logfile,
//        !(DO_BRANCHING_WITH_SICS && DO_CUTSON));
//    printBBInfo(cleanedSICStrBBInfoCplexCutsOff, params.logfile,
//        !(DO_BRANCHING_WITH_SICS && DO_CUTSOFF));
//  } /* cplex */
//  else {
//    for (int i = 0; i < (int) BB_INFO_CONTENTS.size() * 6; i++) {
//      writeEntryToLog("", params.logfile);
//    }
//  }
//  // Gurobi
//  if (param.getParamVal(ParamIndices::BB_STRATEGY_PARAM_IND) & BB_Strategy_Options::gurobi) {
//    printBBInfo(cleanedBBInfoGurobiDefault, params.logfile,
//        !DO_DEFAULT);
//    printBBInfo(cleanedBBInfoGurobiCutsOn, params.logfile,
//        !DO_CUTSON);
//    printBBInfo(cleanedBBInfoGurobiCutsOff, params.logfile,
//        !DO_CUTSOFF);
//    printBBInfo(cleanedSICStrBBInfoGurobiDefault, params.logfile,
//        !(DO_BRANCHING_WITH_SICS && DO_DEFAULT));
//    printBBInfo(cleanedSICStrBBInfoGurobiCutsOn, params.logfile,
//        !(DO_BRANCHING_WITH_SICS && DO_CUTSON));
//    printBBInfo(cleanedSICStrBBInfoGurobiCutsOff, params.logfile,
//        !(DO_BRANCHING_WITH_SICS && DO_CUTSOFF));
//  } /* gurobi */
//  else {
//    for (int i = 0; i < (int) BB_INFO_CONTENTS.size() * 6; i++) {
//      writeEntryToLog("", params.logfile);
//    }
//  }
//  writeEntryToLog(GlobalVariables::timeStats.get_total_time(INIT_SOLVE_TIME),
//      params.logfile);
} /* performCleaning */

/**
 * Makes sure no variable bounds can be tightened,
 * including via strong branching
 */
bool cleanProblem(const VPCParameters& params, OsiClpSolverInterface* solver,
    int& numBoundsChanged, int& numSBFixed) {
  bool is_clean = true;

  const int numCols = solver->getNumCols();
//  std::vector<double> origObjCoeff(solver->getObjCoefficients(),
//      solver->getObjCoefficients() + numCols);
  std::vector<double> origColSolution(solver->getColSolution(),
        solver->getColSolution() + numCols);
  std::vector<double> obj(numCols, 0.);

  // Set up solver for checking bounds
  OsiClpSolverInterface* boundSolver =
      dynamic_cast<OsiClpSolverInterface*>(solver->clone());
  boundSolver->setObjective(obj.data());
  boundSolver->setObjSense(solver->getObjSense());
  double objOffset = 0.;
  solver->getDblParam(OsiDblParam::OsiObjOffset, objOffset);
  boundSolver->setDblParam(OsiDblParam::OsiObjOffset, objOffset);

  // For strong branching
  setupClpForStrongBranching(solver);
  solver->enableFactorization();
  solver->markHotStart();
  for (int col = 0; col < numCols; col++) {
    // If col is integral, ensure bounds are integral
    // If the integer variable is fractional, we also strong branch
    if (solver->isInteger(col)) {
      // Is the variable fractional in the solution?
      if (!isVal(origColSolution[col], std::floor(origColSolution[col]),
          params.get(doubleConst::AWAY))
          && !isVal(origColSolution[col], std::ceil(origColSolution[col]),
              params.get(doubleConst::AWAY))) {
        const double origLB = solver->getColLower()[col];
        const double origUB = solver->getColUpper()[col];
        bool downBranchFeasible = true, upBranchFeasible = true;

        // Check down branch
        solver->setColUpper(col, std::floor(origColSolution[col]));
        //solver->solveFromHotStart();
        solveFromHotStart(solver, col, true, origUB, std::floor(origColSolution[col]));
        if (solver->isProvenPrimalInfeasible()) {
          downBranchFeasible = false;
        }
        solver->setColUpper(col, origUB);

        // Return to original state
        solver->solveFromHotStart();

        // Check up branch
        solver->setColLower(col, std::ceil(origColSolution[col]));
        //solver->solveFromHotStart();
        solveFromHotStart(solver, col, false, origLB, std::ceil(origColSolution[col]));
        if (solver->isProvenPrimalInfeasible()) {
          upBranchFeasible = false;
        }
        solver->setColLower(col, origLB);

        // Return to original state
        solver->solveFromHotStart();

        // Check if some side of the split is infeasible
        if (!downBranchFeasible || !upBranchFeasible) {
          if (!downBranchFeasible && !upBranchFeasible) {
            // Infeasible problem
            error_msg(errorstring,
                "Infeasible problem due to integer variable %d (value %e).\n",
                col, origColSolution[col]);
            writeErrorToLog(errorstring, params.logfile);
            exit(1);
          }
          numSBFixed++;
          is_clean = false;
          if (!downBranchFeasible) {
            boundSolver->setColLower(col, std::ceil(origColSolution[col]));
          } else {
            boundSolver->setColUpper(col, std::floor(origColSolution[col]));
          }
        }
      } /* end strong branching */

      // Check integrality of LB
      if (!isVal(solver->getColLower()[col],
          std::floor(solver->getColLower()[col]))) {
        numBoundsChanged++;
        is_clean = false;
        boundSolver->setColLower(col, std::ceil(solver->getColLower()[col]));
      }

      // Check integrality of UB
      if (!isVal(solver->getColUpper()[col],
          std::ceil(solver->getColUpper()[col]))) {
        numBoundsChanged++;
        is_clean = false;
        boundSolver->setColUpper(col, std::floor(solver->getColUpper()[col]));
      }
    } /* check if integer */

    // Check if LB can be tightened
    boundSolver->setObjCoeff(col, 1.);
    boundSolver->resolve();
    if (checkSolverOptimality(boundSolver, false)) {
      if (greaterThanVal(boundSolver->getColSolution()[col],
          solver->getColLower()[col])) {
        numBoundsChanged++;
        is_clean = false;
        boundSolver->setColLower(col, boundSolver->getColSolution()[col]);
      }
    }

    // Check if UB can be tightened
    boundSolver->setObjCoeff(col, -1.);
    boundSolver->resolve();
    if (checkSolverOptimality(boundSolver, false)) {
      if (lessThanVal(boundSolver->getColSolution()[col],
          solver->getColUpper()[col])) {
        numBoundsChanged++;
        is_clean = false;
        boundSolver->setColUpper(col, boundSolver->getColSolution()[col]);
      }
    }

    // Reset
    boundSolver->setObjCoeff(col, 0.);
  } /* end iterating over columns */
  solver->unmarkHotStart();
  solver->disableFactorization();

  // Change upper and lower bounds wherever needed
  for (int col = 0; col < solver->getNumCols(); col++) {
    solver->setColLower(col, boundSolver->getColLower()[col]);
    solver->setColUpper(col, boundSolver->getColUpper()[col]);
  } /* end iterating over columns to update solver */

  if (boundSolver) {
    delete boundSolver;
  }

  return is_clean;
} /* cleanProblem */
