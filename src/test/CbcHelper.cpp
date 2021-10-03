/**
 * @file CbcHelper.cpp
 * @author A. M. Kazachkov
 * @date 2021-Oct-3
 */
#include "CbcHelper.hpp"

// Project files
#include "BBHelper.hpp"
#include "CutHelper.hpp" // applyCuts
#include "SolverHelper.hpp"
#include "SolverInterface.hpp"
#include "VPCParameters.hpp"
using namespace VPCParametersNamespace;

// COIN-OR
#include <CoinTime.hpp>
#include <OsiCuts.hpp>

#ifdef USE_CBC
#include <CbcModel.hpp>

// Variable selection
#include <CbcBranchDefaultDecision.hpp>
#include <OsiChooseVariable.hpp>

// General strategy
#include <CbcStrategy.hpp>
#endif // USE_CBC

#ifdef USE_CBC
void setStrategyForBBTestCbc(
    const VPCParametersNamespace::VPCParameters& params,
    CbcModel* const cbc_model,
    int seed = -1) {
  if (seed < 0) seed = params.get(intParam::RANDOM_SEED);
  // Parameters that should always be set
  cbc_model->setMaximumSeconds(params.get(doubleParam::BB_TIMELIMIT)); // time limit
  cbc_model->setRandomSeed(seed); // random seed

  int strategy = params.get(intParam::BB_STRATEGY);
  if (strategy <= 0) {
    // Default strategy
    CbcStrategyDefault strategy;
    cbc_model->setStrategy(strategy);
    /*
    CbcStrategyDefault strategy(-1);
    strategy.setupPreProcessing(-1,0);
    cbc_model->setStrategy(strategy);
    cbc_model->setMaximumCutPassesAtRoot(0);
    cbc_model->setMaximumCutPasses(0);
    cbc_model->setWhenCuts(0);
    */
  } else {
    if (use_bb_option(strategy, BB_Strategy_Options::user_cuts)) {
      // Make sure dual reductions are off
    }

    // Turn off all cuts
    if (use_bb_option(strategy, BB_Strategy_Options::all_cuts_off)) {
      cbc_model->setMaximumCutPassesAtRoot(0);
      cbc_model->setMaximumCutPasses(0);
      cbc_model->setWhenCuts(0);
    }

    // Presolve
    if (use_bb_option(strategy, BB_Strategy_Options::presolve_off)) {
      cbc_model->setTypePresolve(0);
    }
  }

  // Check if we should use strong branching
  // Not sure this works when using StrategyDefault as well...
  if (use_bb_option(std::abs(strategy), BB_Strategy_Options::strong_branching_on)) {
    CbcBranchDefaultDecision branch;
    OsiChooseStrong choose;
    choose.setNumberStrong(cbc_model->solver()->getNumCols());
    choose.setNumberBeforeTrusted(std::numeric_limits<int>::max());
    branch.setChooseMethod(choose);
    cbc_model->setBranchingMethod(branch);
  }
} /* setStrategyForBBTestCbc */

/**
 * Perform branch-and-bound without cuts
 */
void doBranchAndBoundNoCuts(const VPCParametersNamespace::VPCParameters& params,
    const OsiSolverInterface* const solver, BBInfo& info) {
#ifdef TRACE
  printf("\nBB with no cuts.\n");
#endif
  const bool test_using_main = false;

  // Set up solver
  SolverInterface* BBSolver;
  BBSolver = dynamic_cast<SolverInterface*>(solver->clone());
  setLPSolverParameters(BBSolver, params.get(intParam::VERBOSITY), params.get(doubleParam::TIMELIMIT));

  // Set up model
  CbcModel* cbc_model = new CbcModel;
  cbc_model->swapSolver(BBSolver);
  cbc_model->setModelOwnsSolver(true);
  setIPSolverParameters(cbc_model, params.get(VERBOSITY));

  // Set up options and run B&B
  if (!test_using_main) {
    setStrategyForBBTestCbc(params, cbc_model);
#ifdef TRACE
    cbc_model->branchAndBound(3);
#else
    cbc_model->branchAndBound(0);
#endif
  } else {
#ifndef CBC_VERSION_210PLUS
    CbcMain0(*cbc_model);
    std::string name, logLevel, presolveOnOff, preprocessOnOff, cutsOnOff, heurOnOff, solveOption;
    name = "BBHelper_doBranchAndBoundNoCuts";
    presolveOnOff = "-presolve=off";
    preprocessOnOff = "-preprocess=off";
    cutsOnOff = "-cuts=off";
    heurOnOff = "-heur=off";
    solveOption = "-solve";
#ifdef TRACE
    logLevel = "-loglevel=3";
#else
    logLevel = "-loglevel=0";
#endif

    int argc = 0;
    const char** cbc_options = new const char*[20];
    cbc_options[argc++] = name.c_str();
    cbc_options[argc++] = logLevel.c_str();
    cbc_options[argc++] = presolveOnOff.c_str();
    cbc_options[argc++] = preprocessOnOff.c_str();
    cbc_options[argc++] = cutsOnOff.c_str();
    cbc_options[argc++] = heurOnOff.c_str();
    cbc_options[argc++] = solveOption.c_str();

    CbcMain1(argc, cbc_options, *cbc_model);
    delete[] cbc_options;
#endif
  }

  // Collect statistics
  info.time = CoinCpuTime()
      - cbc_model->getDblParam(CbcModel::CbcStartSeconds);
  info.obj = cbc_model->getObjValue();
  info.iters = cbc_model->getIterationCount();
  info.nodes = cbc_model->getNodeCount();
  //info.bound = cbc_model->getCutoff();
  //info.bound = cbc_model->getDblParam(CbcModel::CbcDblParam::CbcCurrentCutoff);

  // Free
  if (cbc_model) {
    delete cbc_model;
    cbc_model = NULL;
  }
} /* doBranchAndBoundNoCuts */

/************************************************************/
/**
 * Perform branch-and-bound using the given cuts, perhaps doing cut selection
 */
void doBranchAndBoundYesCuts(const VPCParametersNamespace::VPCParameters& params,
    const OsiSolverInterface* const solver, BBInfo& info, const OsiCuts& structCuts,
    const bool doCutSelection, const int numCutsToAddPerRound,
    const int maxRounds, const std::string logstring) {
//#ifdef TRACE
  printf("%s", logstring.c_str());
//#endif
  const bool test_using_main = false;

  // Check that there are cuts
  const int numCuts = structCuts.sizeCuts();
  if (numCuts == 0) {
    info.nodes = 0;
    info.time = 0.;
    return;
  }

  // Set up solver
  SolverInterface* BBSolver;
  BBSolver = dynamic_cast<SolverInterface*>(solver->clone());
  setLPSolverParameters(BBSolver, params.get(intParam::VERBOSITY), params.get(doubleParam::TIMELIMIT));

  // Apply cuts
//  if (doCutSelection) {
//    applyCutsInRounds(BBSolver, structCuts, numCutsToAddPerRound, maxRounds);
//  } else {
    applyCutsCustom(BBSolver, structCuts, params.logfile);
//  }

  // Set up model
  CbcModel* cbc_model = new CbcModel;
  cbc_model->swapSolver(BBSolver);
  cbc_model->setModelOwnsSolver(true);
  setIPSolverParameters(cbc_model, params.get(VERBOSITY));

  // Set up options and run B&B
  if (!test_using_main) {
    setStrategyForBBTestCbc(params, cbc_model);
#ifdef TRACE
    cbc_model->branchAndBound(3);
#else
    cbc_model->branchAndBound(0);
#endif
  } else {
#ifndef CBC_VERSION_210PLUS
    CbcMain0(*cbc_model);
    std::string name, logLevel, presolveOnOff, preprocessOnOff, cutsOnOff, heurOnOff, solveOption;
    name = "BBHelper_doBranchAndBoundYesCuts";
    presolveOnOff = "-presolve=off";
    preprocessOnOff = "-preprocess=off";
    cutsOnOff = "-cuts=off";
    heurOnOff = "-heur=off";
    solveOption = "-solve";
#ifdef TRACE
    logLevel = "-loglevel=3";
#else
    logLevel = "-loglevel=0";
#endif

    int argc = 0;
    const char** cbc_options = new const char*[20];
    cbc_options[argc++] = name.c_str();
    cbc_options[argc++] = logLevel.c_str();
    cbc_options[argc++] = presolveOnOff.c_str();
    cbc_options[argc++] = preprocessOnOff.c_str();
    cbc_options[argc++] = cutsOnOff.c_str();
    cbc_options[argc++] = heurOnOff.c_str();
    cbc_options[argc++] = solveOption.c_str();

    CbcMain1(argc, cbc_options, *cbc_model);
    delete[] cbc_options;
#endif
  }

  // Collect statistics
  info.time = CoinCpuTime()
      - cbc_model->getDblParam(CbcModel::CbcStartSeconds);
  info.obj = cbc_model->getObjValue();
  info.iters = cbc_model->getIterationCount();
  info.nodes = cbc_model->getNodeCount();
  //info.bound = cbc_model->getCutoff();

  // Free
  if (cbc_model) {
    delete cbc_model;
    cbc_model = NULL;
  }
} /* doBranchAndBoundYesCuts */
#endif // USE_CBC
