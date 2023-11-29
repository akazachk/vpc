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
#include "VPCParameters.hpp"
using namespace VPCParametersNamespace;

// COIN-OR
#include <CoinTime.hpp>
#include <OsiCuts.hpp>

#ifdef USE_CBC
#include <CbcModel.hpp>
#include <CbcSolver.hpp>  // CbcParam, CbcSolverUsefulData
#include <CbcSolverHeuristics.hpp> // doHeuristics

// Variable selection
#include <CbcBranchDefaultDecision.hpp>
#include <OsiChooseVariable.hpp>

// General strategy
#include <CbcStrategy.hpp>

// Cut Store
#include "CglStoredVpc.hpp" // CglStoredVpc
#endif // USE_CBC

#ifdef USE_CBC
void setStrategyForBBTestCbc(
    const VPCParametersNamespace::VPCParameters& params,
    const int given_strategy,
    CbcModel* const cbc_model,
    const double best_bound,
    int seed = -1) {
  if (seed < 0) seed = params.get(intParam::RANDOM_SEED);
  // Parameters that should always be set
  cbc_model->setMaximumSeconds(params.get(doubleParam::BB_TIMELIMIT)); // time limit
  cbc_model->setRandomSeed(seed); // random seed

  int strategy = given_strategy;

  CbcStrategyDefault cbc_strategy;
  cbc_model->setStrategy(cbc_strategy);

  if (use_bb_option(strategy, BB_Strategy_Options::heuristics_on)){
    // set parameters repsonsible for turning on default heuristics (according to CbcSolver)
    CbcSolverUsefulData cbcData;
    cbcData[CbcParam::ROUNDING]->setVal("on");
    cbcData[CbcParam::FPUMP]->setVal("on");
    cbcData[CbcParam::GREEDY]->setVal("on");
    cbcData[CbcParam::FPUMPITS]->setVal(30);
    cbcData[CbcParam::FPUMPTUNE]->setVal(1005043);
    // turn on heuristics - doesn't actually do them until solve
    doHeuristics(cbc_model, 1, cbcData, cbcData.noPrinting(), 1005043);
  }

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

  // set bound
  if (use_bb_option(strategy, BB_Strategy_Options::use_best_bound)) {
    verify(!isInfinity(std::abs(best_bound)), "Best bound must be finite if applied.\n");
    // sets primal bound - add a small tolerance so CBC finds the eventual solution
    cbc_model->setCutoff(best_bound * (1 + 1e-5));
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

/// @brief This is so user can trap events and capture useful statistics
/// CbcModel model_ is available as well as anything else you care to pass in.
class CbcUserCutEventHandler : public CbcEventHandler {
  public:
    VPCParameters params;
    int num_vars;
    double obj_offset;
    const OsiCuts* cuts;
    bool addAsLazy;
    double first_lp_opt;
    BBInfo info;

    /**@name Overrides */
    //@{
    /// @brief Custom event handler that keeps history of the tree
    /// @details Returns CbcAction based on one of the following events
    ///   enum CbcEvent {
    ///     *! Processing of the current node is complete. amk: Called within CbcModel::doOneNode.
    ///     node = 200,
    ///     *! A tree status interval has arrived.
    ///     *! amk: Called within CbcModel::branchAndBound.
    ///     *! amk: At this point, we know the next step will be picking the next leaf node to branch on using tree_->bestNode(cutoff);
    ///     *! amk: We may not reach a treeStatus interval if stoppingCriterionReached() is true after doOneNode is called (below the treeStatus event).
    ///     *! amk: If stoppingCriterionReached(), the tree is destroyed in tree_->cleanTree().
    ///     treeStatus,
    ///     *! A solution has been found.
    ///     solution,
    ///     *! A heuristic solution has been found.
    ///     heuristicSolution,
    ///     *! A solution will be found unless user takes action (first check).
    ///     beforeSolution1,
    ///     *! A solution will be found unless user takes action (thorough check).
    ///     beforeSolution2,
    ///     *! After failed heuristic.
    ///     afterHeuristic,
    ///     *! On entry to small branch and bound.
    ///     smallBranchAndBound,
    ///     *! After a pass of heuristic.
    ///     heuristicPass,
    ///     *! When converting constraints to cuts.
    ///     convertToCuts,
    ///     *! Having generated cuts, allows user to think. // amk: added in Cbc 2.10
    ///     generatedCuts,
    ///     *! End of search.
    ///     endSearch
    ///   };
    virtual CbcAction event(CbcEvent whichEvent) {
      // don't track anything if we're not in the main branch and bound loop
      if ((model_->specialOptions() & 2048) != 0){
        return noAction;
      }
      const int num_nodes = model_->getNodeCount();
      // model_->getBestPossibleObjValue() is wrong when getting objective after cuts
      const double objValue = model_->solver()->getObjValue();
      this->info.bound = objValue;
      if (num_nodes == 0) {
        if (this->info.root_passes == 1) {
          this->info.first_cut_pass = objValue;
        }
      }

      // unsure what this is supposed to do, but when you add a cut store, it's tripped
//      { /// DEBUG DEBUG
//        if (model_->currentNumberCuts() > 0 || model_->globalCuts()->sizeRowCuts() > 0) {
//          exit(1);
//        }
//      } /// DEBUG DEBUG

      if (whichEvent == CbcEventHandler::treeStatus) {
      }
      else if (whichEvent == CbcEventHandler::node) {
        if (num_nodes <= 1) {
          this->info.root_iters = model_->getIterationCount();
          if (this->info.root_passes == 0) {
            // In case no cuts are ever generated
            this->info.first_cut_pass = objValue;
            this->info.last_cut_pass = objValue;
          }
        }
      }
      else if (whichEvent == CbcEventHandler::beforeSolution2) {
        this->info.last_sol_time = model_->getCurrentSeconds();
      }
#ifdef CBC_VERSION_210PLUS
      else if (whichEvent == CbcEventHandler::generatedCuts) {
        if (num_nodes == 0) {
          // When root_passes = 0, we have not applied any cuts yet
          this->info.root_passes++;
          this->info.last_cut_pass = objValue;
          this->info.root_time = model_->getCurrentSeconds();
        }
      }
#endif // CBC_VERSION_210PLUS
      return CbcEventHandler::CbcAction::noAction;
    } /* event */
    //@}

    /**@name Constructors, destructor, etc. */
    //@{
    /// @brief Default constructor
    CbcUserCutEventHandler() {
      initialize(NULL);
    } /* constructor (default) */

    /// @brief Special constructors
    CbcUserCutEventHandler(
        const VPCParameters& params,
        int num_vars,
        double obj_offset,
        const OsiCuts* cutsToAdd,
        const bool lazyFlag = false
        ) : params(params), num_vars(num_vars), obj_offset(obj_offset), cuts(cutsToAdd), addAsLazy(lazyFlag) {
    } /* constructor (custom) */

    /// @brief Constructor with pointer to model (redundant as setEventHandler does)
    CbcUserCutEventHandler(CbcModel* model) : CbcEventHandler(model) {
      try {
        initialize(dynamic_cast<CbcUserCutEventHandler*>(model->getEventHandler()));
      } catch (std::exception& e) {
        initialize(NULL);
      }
    } /* constructor (pointer to CbcModel) */

    /// @brief Destructor
    virtual ~CbcUserCutEventHandler() {
    } /* destructor */

    /// @brief The copy constructor
    CbcUserCutEventHandler(const CbcUserCutEventHandler& rhs) : CbcEventHandler(rhs) {
      initialize(&rhs);
    } /* copy constructor */

    /// @brief Assignment
    CbcUserCutEventHandler& operator=(const CbcUserCutEventHandler& rhs) {
      if (this != &rhs) {
        CbcEventHandler::operator=(rhs);
      }
      initialize(&rhs);
      return *this;
    } /* assignment operator */

    /// @brief Clone
    virtual CbcEventHandler* clone() const {
      return new CbcUserCutEventHandler(*this);
    } /* clone */
    //@}

  protected:
    ///@{
    /// @name Helper methods

    /// @brief Copy our stuff
    void initialize(const CbcUserCutEventHandler* const rhs) {
      if (rhs) {
        this->params = rhs->params;
        this->num_vars = rhs->num_vars;
        this->obj_offset = rhs->obj_offset;
        this->cuts = rhs->cuts;
        this->addAsLazy = rhs->addAsLazy;
        this->first_lp_opt = rhs->first_lp_opt;
        this->info = rhs->info;
      } else {
        this->num_vars = -1;
        this->obj_offset = 0.;
        this->addAsLazy = false;
        this->info.obj = params.get(doubleParam::INF);
        this->first_lp_opt = -1 * params.get(doubleParam::INF);
        this->info.root_passes = 0;
        this->info.first_cut_pass = -1 * params.get(doubleParam::INF);
        this->info.last_cut_pass = -1 * params.get(doubleParam::INF);
        this->info.root_iters = 0;
        this->info.root_time = 0.;
        this->info.last_sol_time = 0.;
      }
    } /* initialize */
    //@}
}; /* CbcUserCutEventHandler */

void doBranchAndBoundWithCbc(
    const VPCParameters& params,
    int strategy,
    CbcModel& model,
    BBInfo& info,
    const double best_bound,
    std::vector<double>* const solution = NULL) {
  printf("\tRunning B&B with Cbc. Strategy: %d. Random seed: %d.\n",
      strategy, params.get(intParam::RANDOM_SEED));
  setStrategyForBBTestCbc(params, strategy, &model, best_bound);

  model.branchAndBound(params.get(intParam::VERBOSITY));

  // Status of problem - 0 finished, 1 stopped, 2 difficulties
  const int optimstatus = model.status();
  if ((optimstatus == 0 || optimstatus == 1 || optimstatus == 2) &&
      (model.isProvenInfeasible() || model.isContinuousUnbounded() || model.isProvenDualInfeasible())) {
    error_msg(errorstring, "CBC: Failed to optimize MIP.\n");
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }

  // Collect statistics
  info.time = CoinCpuTime()
    - model.getDblParam(CbcModel::CbcStartSeconds);
  info.obj = model.getObjValue();
  info.bound = model.getBestPossibleObjValue();
  info.iters = model.getIterationCount();
  info.nodes = model.getNodeCount();
  //info.last_cut_pass = model.rootObjectiveAfterCuts();

#ifdef TRACE
    printf("Cbc: Solution value: %s.\n", stringValue(info.obj, "%1.6f", params.get(doubleParam::INF)).c_str());
    printf("Cbc: Best bound: %s.\n", stringValue(info.bound, "%1.6f", params.get(doubleParam::INF)).c_str());
    printf("Cbc: Number iterations: %s.\n", stringValue(info.iters, "%ld").c_str());
    printf("Cbc: Number nodes: %s.\n", stringValue(info.nodes, "%ld").c_str());
    printf("Cbc: Time: %f.\n", info.time);
#endif
} /* doBranchAndBoundWithCbc (CbcModel) */

void doBranchAndBoundWithUserCutsCbc(
    const VPCParameters& params,
    int strategy,
    CbcModel& model,
    const OsiCuts* cuts,
    BBInfo& info,
    const double best_bound,
    const bool addAsLazy,
    std::vector<double>* const solution = NULL) {
  // Ensure that user cuts setting is enabled
  if (!use_bb_option(strategy, BB_Strategy_Options::user_cuts)) {
    warning_msg(warnstring,
        "Need to use user_cuts option; strategy currently: %d.\n", strategy);
    strategy = enable_bb_option(strategy, BB_Strategy_Options::user_cuts);
  }

  // Add custom event handler
  double objOffset = 0.;
  model.solver()->getDblParam(OsiDblParam::OsiObjOffset, objOffset);
  CbcUserCutEventHandler* cb = new CbcUserCutEventHandler(params,
      model.getNumCols(), objOffset, cuts, addAsLazy);

  // Apply cuts
  CglStoredVpc store;
  if (cuts) {
    // if cuts are on, make sure vpcs are used as a cut generator instead of just being appended to the model
    // sometimes vpcs don't play nicely with default cuts and we want to give CBC a means of handling that
    if (!use_bb_option(strategy, BB_Strategy_Options::all_cuts_off)){
      for (int i = 0; i < cuts->sizeRowCuts(); i++) {
        store.addCut(cuts->rowCut(i));
      }
      model.addCutGenerator(&store, 1, "StoredVPCs");
    } else {
    model.solver()->applyCuts(*cuts);
    }
    //applyCutsCustom(model.solver(), *cuts);
  }

  if (!cb) {
    error_msg(errstr, "Could not get event handler.\n");
    writeErrorToLog(errstr, params.logfile);
    exit(1);
  }

  // CbcModel makes clone of event handler
  model.passInEventHandler(cb);
  delete cb;
  cb = NULL;

  // Finally, run B&B
  doBranchAndBoundWithCbc(params, strategy, model, info, best_bound, solution);

  // Get back eventHandler
  try {
    cb = dynamic_cast<CbcUserCutEventHandler*>(model.getEventHandler());
  } catch (std::exception& e) {
  }
  if (!cb) {
    error_msg(errstr, "Could not get event handler.\n");
    writeErrorToLog(errstr, params.logfile);
    exit(1);
  }

  // Save information
  info.root_passes = cb->info.root_passes;
  if (info.root_passes > 0) {
    info.first_cut_pass = cb->info.first_cut_pass; // second because first is lp opt val
    info.last_cut_pass = cb->info.last_cut_pass;
    //info.last_cut_pass = model.rootObjectiveAfterCuts();
    info.root_iters = cb->info.root_iters;
    info.root_time = cb->info.root_time;
    info.last_sol_time = cb->info.last_sol_time;
  } else {
    // I guess this can happen if we solve during presolve,
    // or if we do no root passes of cuts
    if (info.nodes == 0) {
      info.first_cut_pass = info.obj;
      info.last_cut_pass = info.obj;
      info.root_iters = info.iters; // all iters were at the root
      info.root_time = info.time; // all time was spent at the root
      info.last_sol_time = cb->info.last_sol_time; // roughly the same as total time in this case
    } else {
      info.first_cut_pass = cb->first_lp_opt;
      info.last_cut_pass = cb->first_lp_opt;
      info.root_iters = cb->info.root_iters;
      info.root_time = cb->info.root_time; // all time was spent at the root
      info.last_sol_time = cb->info.last_sol_time; // roughly the same as total time in this case
    }
  }
} /* doBranchAndBoundWithUserCutsCbc (CbcModel) */

void doBranchAndBoundWithCbc(const VPCParameters& params, int strategy,
    const char* f_name, BBInfo& info, const double best_bound,
    std::vector<double>* const solution) {
#ifdef TRACE
  printf("\n## Reading from file into Cbc. ##\n");
#endif
  // Read in file
  OsiSolverInterface* BBSolver;
  initializeSolver(BBSolver, params.get(FILENAME), params.get(VERBOSITY), params.get(TIMELIMIT), params.logfile);

  // Set up model
  CbcModel* cbc_model = new CbcModel;
  cbc_model->swapSolver(BBSolver);
  cbc_model->setModelOwnsSolver(true);
  setIPSolverParameters(cbc_model, params.get(VERBOSITY));

  if (use_bb_option(strategy, BB_Strategy_Options::user_cuts)) {
    doBranchAndBoundWithUserCutsCbc(params, strategy, *cbc_model, NULL, info, best_bound, false, solution);
  } else {
    doBranchAndBoundWithCbc(params, strategy, *cbc_model, info, best_bound, solution);
  }

  if (BBSolver && !cbc_model->modelOwnsSolver()) { delete BBSolver; }
  if (cbc_model) { delete cbc_model; }
} /* doBranchAndBoundWithCbc (filename) */

void doBranchAndBoundWithCbc(const VPCParameters& params, int strategy,
    const OsiSolverInterface* const solver, BBInfo& info,
    const double best_bound, std::vector<double>* const solution) {
  //std::string f_name;
  //createTmpFileCopy(params, solver, f_name);
  //doBranchAndBoundWithCbc(params, strategy, f_name.c_str(), info, best_bound, solution);
  //remove(f_name.c_str()); // remove temporary file

  OsiSolverInterface* BBSolver = solver->clone();
  setLPSolverParameters(BBSolver, params.get(VERBOSITY), params.get(TIMELIMIT));

  // Set up model
  CbcModel* cbc_model = new CbcModel;
  cbc_model->swapSolver(BBSolver);
  cbc_model->setModelOwnsSolver(true);
  setIPSolverParameters(cbc_model, params.get(VERBOSITY));

  if (use_bb_option(strategy, BB_Strategy_Options::user_cuts)) {
    doBranchAndBoundWithUserCutsCbc(params, strategy, *cbc_model, NULL, info, best_bound, false, solution);
  } else {
    doBranchAndBoundWithCbc(params, strategy, *cbc_model, info, best_bound, solution);
  }

  if (BBSolver && !cbc_model->modelOwnsSolver()) { delete BBSolver; }
  if (cbc_model) { delete cbc_model; }
} /* doBranchAndBoundWithCbc (Osi) */

void doBranchAndBoundWithUserCutsCbc(const VPCParameters& params,
    int strategy, const char* f_name, const OsiCuts* cuts, BBInfo& info,
    const double best_bound, const bool addAsLazy) {
#ifdef TRACE
  printf("\n## Reading from file into Cbc and adding user cuts. ##\n");
#endif
  // Read in file
  OsiSolverInterface* BBSolver;
  initializeSolver(BBSolver, params.get(FILENAME), params.get(VERBOSITY), params.get(TIMELIMIT), params.logfile);

  // Set up model
  CbcModel* cbc_model = new CbcModel;
  cbc_model->swapSolver(BBSolver);
  cbc_model->setModelOwnsSolver(true);
  setIPSolverParameters(cbc_model, params.get(VERBOSITY));

  doBranchAndBoundWithUserCutsCbc(params, strategy, *cbc_model, cuts, info, best_bound, addAsLazy);

  if (BBSolver && !cbc_model->modelOwnsSolver()) { delete BBSolver; }
  if (cbc_model) { delete cbc_model; }
} /* doBranchAndBoundWithUserCutsCbc (filename) */

void doBranchAndBoundWithUserCutsCbc(const VPCParameters& params, int strategy,
    const OsiSolverInterface* const solver, const OsiCuts* cuts, BBInfo& info,
    const double best_bound, const bool addAsLazy) {
  //std::string f_name;
  //createTmpFileCopy(params, solver, f_name);
  //doBranchAndBoundWithUserCutsCbc(params, strategy, f_name.c_str(), cuts, info, best_bound, addAsLazy);
  //remove(f_name.c_str()); // remove temporary file

  OsiSolverInterface* BBSolver = solver->clone();
  setLPSolverParameters(BBSolver, params.get(VERBOSITY), params.get(TIMELIMIT));

  // Set up model
  CbcModel* cbc_model = new CbcModel;
  cbc_model->swapSolver(BBSolver);
  cbc_model->setModelOwnsSolver(true);
  setIPSolverParameters(cbc_model, params.get(VERBOSITY));

  doBranchAndBoundWithUserCutsCbc(params, strategy, *cbc_model, cuts, info, best_bound, addAsLazy);

  if (BBSolver && !cbc_model->modelOwnsSolver()) { delete BBSolver; }
  if (cbc_model) { delete cbc_model; }
} /* doBranchAndBoundWithUserCutsCbc (Osi) */

///**
// * Perform branch-and-bound without cuts
// */
//void doBranchAndBoundNoCuts(const VPCParametersNamespace::VPCParameters& params,
//    const OsiSolverInterface* const solver, BBInfo& info) {
//#ifdef TRACE
//  printf("\nBB with no cuts.\n");
//#endif
//  const bool test_using_main = false;
//
//  // Set up solver
//  SolverInterface* BBSolver;
//  BBSolver = dynamic_cast<SolverInterface*>(solver->clone());
//  setLPSolverParameters(BBSolver, params.get(intParam::VERBOSITY), params.get(doubleParam::TIMELIMIT));
//
//  // Set up model
//  CbcModel* cbc_model = new CbcModel;
//  cbc_model->swapSolver(BBSolver);
//  cbc_model->setModelOwnsSolver(true);
//  setIPSolverParameters(cbc_model, params.get(VERBOSITY));
//
//  // Set up options and run B&B
//  if (!test_using_main) {
//    setStrategyForBBTestCbc(params, cbc_model);
//#ifdef TRACE
//    cbc_model->branchAndBound(3);
//#else
//    cbc_model->branchAndBound(0);
//#endif
//  } else {
//#ifndef CBC_VERSION_210PLUS
//    CbcMain0(*cbc_model);
//    std::string name, logLevel, presolveOnOff, preprocessOnOff, cutsOnOff, heurOnOff, solveOption;
//    name = "BBHelper_doBranchAndBoundNoCuts";
//    presolveOnOff = "-presolve=off";
//    preprocessOnOff = "-preprocess=off";
//    cutsOnOff = "-cuts=off";
//    heurOnOff = "-heur=off";
//    solveOption = "-solve";
//#ifdef TRACE
//    logLevel = "-loglevel=3";
//#else
//    logLevel = "-loglevel=0";
//#endif
//
//    int argc = 0;
//    const char** cbc_options = new const char*[20];
//    cbc_options[argc++] = name.c_str();
//    cbc_options[argc++] = logLevel.c_str();
//    cbc_options[argc++] = presolveOnOff.c_str();
//    cbc_options[argc++] = preprocessOnOff.c_str();
//    cbc_options[argc++] = cutsOnOff.c_str();
//    cbc_options[argc++] = heurOnOff.c_str();
//    cbc_options[argc++] = solveOption.c_str();
//
//    CbcMain1(argc, cbc_options, *cbc_model);
//    delete[] cbc_options;
//#endif
//  }
//
//  // Collect statistics
//  info.time = CoinCpuTime()
//    - cbc_model->getDblParam(CbcModel::CbcStartSeconds);
//  info.obj = cbc_model->getObjValue();
//  info.iters = cbc_model->getIterationCount();
//  info.nodes = cbc_model->getNodeCount();
//  //info.bound = cbc_model->getCutoff();
//  //info.bound = cbc_model->getDblParam(CbcModel::CbcDblParam::CbcCurrentCutoff);
//
//  // Free
//  if (cbc_model) {
//    delete cbc_model;
//    cbc_model = NULL;
//  }
//} /* doBranchAndBoundNoCuts */
//
///**
// * Perform branch-and-bound using the given cuts, perhaps doing cut selection
// */
//void doBranchAndBoundYesCuts(const VPCParametersNamespace::VPCParameters& params,
//    const OsiSolverInterface* const solver, BBInfo& info, const OsiCuts& structCuts,
//    const bool doCutSelection, const int numCutsToAddPerRound,
//    const int maxRounds, const std::string logstring) {
//  //#ifdef TRACE
//  printf("%s", logstring.c_str());
//  //#endif
//  const bool test_using_main = false;
//
//  // Check that there are cuts
//  const int numCuts = structCuts.sizeCuts();
//  if (numCuts == 0) {
//    info.nodes = 0;
//    info.time = 0.;
//    return;
//  }
//
//  // Set up solver
//  SolverInterface* BBSolver;
//  BBSolver = dynamic_cast<SolverInterface*>(solver->clone());
//  setLPSolverParameters(BBSolver, params.get(intParam::VERBOSITY), params.get(doubleParam::TIMELIMIT));
//
//  // Apply cuts
//  //  if (doCutSelection) {
//  //    applyCutsInRounds(BBSolver, structCuts, numCutsToAddPerRound, maxRounds);
//  //  } else {
//  applyCutsCustom(BBSolver, structCuts, params.logfile);
//  //  }
//
//  // Set up model
//  CbcModel* cbc_model = new CbcModel;
//  cbc_model->swapSolver(BBSolver);
//  cbc_model->setModelOwnsSolver(true);
//  setIPSolverParameters(cbc_model, params.get(VERBOSITY));
//
//  // Set up options and run B&B
//  if (!test_using_main) {
//    setStrategyForBBTestCbc(params, cbc_model);
//#ifdef TRACE
//    cbc_model->branchAndBound(3);
//#else
//    cbc_model->branchAndBound(0);
//#endif
//  } else {
//#ifndef CBC_VERSION_210PLUS
//    CbcMain0(*cbc_model);
//    std::string name, logLevel, presolveOnOff, preprocessOnOff, cutsOnOff, heurOnOff, solveOption;
//    name = "BBHelper_doBranchAndBoundYesCuts";
//    presolveOnOff = "-presolve=off";
//    preprocessOnOff = "-preprocess=off";
//    cutsOnOff = "-cuts=off";
//    heurOnOff = "-heur=off";
//    solveOption = "-solve";
//#ifdef TRACE
//    logLevel = "-loglevel=3";
//#else
//    logLevel = "-loglevel=0";
//#endif
//
//    int argc = 0;
//    const char** cbc_options = new const char*[20];
//    cbc_options[argc++] = name.c_str();
//    cbc_options[argc++] = logLevel.c_str();
//    cbc_options[argc++] = presolveOnOff.c_str();
//    cbc_options[argc++] = preprocessOnOff.c_str();
//    cbc_options[argc++] = cutsOnOff.c_str();
//    cbc_options[argc++] = heurOnOff.c_str();
//    cbc_options[argc++] = solveOption.c_str();
//
//    CbcMain1(argc, cbc_options, *cbc_model);
//    delete[] cbc_options;
//#endif
//  }
//
//  // Collect statistics
//  info.time = CoinCpuTime()
//    - cbc_model->getDblParam(CbcModel::CbcStartSeconds);
//  info.obj = cbc_model->getObjValue();
//  info.iters = cbc_model->getIterationCount();
//  info.nodes = cbc_model->getNodeCount();
//  //info.bound = cbc_model->getCutoff();
//
//  // Free
//  if (cbc_model) {
//    delete cbc_model;
//    cbc_model = NULL;
//  }
//} /* doBranchAndBoundYesCuts */
#endif // USE_CBC
