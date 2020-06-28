/**
 * @file GurobiHelper.cpp
 * @author A. M. Kazachkov
 * @date 2019-Feb-28
 */
#include "GurobiHelper.hpp"

//#include <cstdio> // for tmpnam

// Project files
#include "BBHelper.hpp"
#include "CutHelper.hpp" // applyCuts
#include "SolverHelper.hpp"
#include "VPCParameters.hpp"
using namespace VPCParametersNamespace;

// COIN-OR
#include <CoinTime.hpp>
#include <OsiCuts.hpp>

// Gurobi
#ifdef USE_GUROBI
#include "gurobi_c++.h"

/**
 * Creates temporary file (in /tmp) so that it can be read by a different solver
 * It does not delete the file
 */
void createTmpFileCopy(const VPCParameters& params, GRBModel& model, std::string& f_name, const std::string add_ext = ".mps.gz") {
  if (f_name.empty()) {
    try {
      createTmpFilename(f_name, add_ext);
    } catch (const std::exception &e) {
      error_msg(errorstring, "Could not generate temp file: %s.\n", e.what());
      writeErrorToLog(errorstring, params.logfile);
      exit(1);
    }
  } else {
    // Ensure f_name has proper extension
    if (!add_ext.empty() && f_name.compare(f_name.size()-add_ext.size(), add_ext.size(), add_ext) != 0) {
      f_name += add_ext;
    }
  }
  model.write(f_name.c_str());
} /* createTmpFileCopy (Gurobi) */

void setStrategyForBBTestGurobi(const VPCParameters& params, const int strategy,
    GRBModel& model, const double best_bound, int seed = -1) {
  if (seed < 0) seed = params.get(intParam::RANDOM_SEED);
  // Parameters that should always be set
  model.set(GRB_DoubleParam_TimeLimit, params.get(doubleConst::BB_TIMELIMIT)); // time limit
  model.set(GRB_IntParam_Threads, 1); // single-threaded
  model.set(GRB_IntParam_Seed, seed); // random seed
//  model.set(GRB_DoubleParam_MIPGap, param.getEPS()); // I guess the default 1e-4 is okay, though it messes up for large objective values

  if (params.get(VERBOSITY) == 0)
    model.set(GRB_IntParam_OutputFlag, 0); // turn off output

  if (strategy <= 0) {
    // Default strategy
  } else {
    if (use_bb_option(strategy, BB_Strategy_Options::user_cuts)) {
      model.set(GRB_IntParam_DualReductions, 0); // disable dual reductions
      model.set(GRB_IntParam_PreCrush, 1); // must be enabled when using user cuts
    }

    // Turn off all cuts
    if (use_bb_option(strategy, BB_Strategy_Options::all_cuts_off)) {
      model.set(GRB_IntParam_Cuts, 0); // turn off all cuts
    }

    // Presolve
    if (use_bb_option(strategy, BB_Strategy_Options::presolve_off)) {
      model.set(GRB_IntParam_Presolve, 0); // disable presolve
    }

    // Heuristics
    if (use_bb_option(strategy, BB_Strategy_Options::heuristics_off)) {
      model.set(GRB_DoubleParam_Heuristics, 0.); // disable heuristics
      //model.set(GRB_IntParam_PumpPasses, 0); // disable feasibility pump; is this turned off by setting heuristic frequency to 0.?
      //model.set(GRB_IntParam_RINS, 0); // disable RINS; is this turned off by setting heuristic frequency to 0.?
    }

    // Feed the solver the best bound provided
    if (use_bb_option(strategy, BB_Strategy_Options::use_best_bound)) {
      if (!isInfinity(std::abs(best_bound))) {
        // BestObjStop: stop when primal bound <= z
        // BestBdStop: stop when dual bound >= z
        // Cutoff: prune subtrees with objective value > z
        //model.set(GRB_DoubleParam_BestObjStop, best_bound + 1e-3); // give the solver the best IP objective value (it is a minimization problem) with a tolerance
        model.set(GRB_DoubleParam_BestBdStop, best_bound - 1e-7); // give the solver the best IP objective value (it is a minimization problem) with a tolerance
        //model.set(GRB_DoubleParam_Cutoff, best_bound + 1e-3); // give the solver the best IP objective value (it is a minimization problem) with a tolerance
      }
    }
  } /* strategy > 0 */

  // Check if we should use strong branching
  if (use_bb_option(std::abs(strategy), BB_Strategy_Options::strong_branching_on)) {
    model.set(GRB_IntParam_VarBranch, 3); // turn on strong branching
  }
} /* setStrategyForBBTestGurobi */

/**
 * User cut callback
 */
class GurobiUserCutCallback : public GRBCallback {
  public:
    VPCParameters params;
    int num_vars;
    double obj_offset;
    GRBVar* grb_vars;
    const OsiCuts* cuts;
    bool addAsLazy;
    double first_lp_opt;
    BBInfo info;

    GurobiUserCutCallback(const VPCParameters& params, int num_vars, double obj_offset,
        GRBVar* vars, const OsiCuts* cutsToAdd, const bool lazyFlag = false)
        : params(params), num_vars(num_vars), obj_offset(obj_offset), grb_vars(vars),
        cuts(cutsToAdd), addAsLazy(lazyFlag) {
      this->info.obj = params.get(doubleConst::INF);
      this->first_lp_opt = -1 * params.get(doubleConst::INF);
      this->info.root_passes = 0;
      this->info.first_cut_pass = -1 * params.get(doubleConst::INF);
      this->info.last_cut_pass = -1 * params.get(doubleConst::INF);
      this->info.root_time = 0.;
      this->info.last_sol_time = 0.;
//      this->numTimesApplied = 0;
    } /* constructor */

    /**
     * Destructor
     */
    ~GurobiUserCutCallback() {
      if (this->grb_vars) {
        delete[] grb_vars;
      }
    } /* destructor */

  protected:
    double getObjValue() {
      double* vals = GRBCallback::getNodeRel(grb_vars, num_vars);
      double obj_val = obj_offset;
      for (int i = 0; i < num_vars; i++) {
        obj_val += vals[i] * grb_vars[i].get(GRB_DoubleAttr_Obj);
      }
      if (vals) {
        delete[] vals;
      }
      return obj_val;
    } /* getObjValue */

    void callback() {
      try {
        if (where == GRB_CB_MIPSOL) {
          const double obj = GRBCallback::getDoubleInfo(GRB_CB_MIPSOL_OBJ);
          if (obj < this->info.obj) { // integer-feasible solution with better solution has been found
            this->info.last_sol_time = GRBCallback::getDoubleInfo(GRB_CB_RUNTIME);
            this->info.obj = obj;
          }
        } else if (where == GRB_CB_MIPNODE) {
          const int num_nodes = GRBCallback::getDoubleInfo(GRB_CB_MIPNODE_NODCNT);
          if (num_nodes > 0) {
            return;
          }
          this->info.root_passes++;
          this->info.root_time = GRBCallback::getDoubleInfo(GRB_CB_RUNTIME);

          if (GRBCallback::getIntInfo(GRB_CB_MIPNODE_STATUS) != GRB_OPTIMAL) {
            return;
          }

          const double objValue = getObjValue();
          this->info.last_cut_pass = objValue;

          // Make sure this is our first entry into the root
          if (this->info.root_passes > 1) {
            if (this->info.root_passes == 2) {
              this->info.first_cut_pass = objValue;
            }
            return;
          }

          this->first_lp_opt = objValue;
          this->info.first_cut_pass  = objValue;

          // Add cuts to the model, one at a time
          if (cuts) {
            for (int cut_ind = 0; cut_ind < cuts->sizeCuts(); cut_ind++) {
              const OsiRowCut* curr_cut = cuts->rowCutPtr(cut_ind);
              const int num_el = curr_cut->row().getNumElements();
              const int* ind = curr_cut->row().getIndices();
              const double* vals = curr_cut->row().getElements();

              // Cannot add it all at once due to the C++ interface being worse than the C one
              //GRBaddconstr(model, num_el, ind, vals, GRB_GREATER_EQUAL, curr_cut->rhs(), NULL);
              GRBLinExpr lhs = 0;
              for (int i = 0; i < num_el; i++) {
                const int curr_var = ind[i];
                const double curr_value = vals[i];
                lhs += curr_value * grb_vars[curr_var];
              }
              // model.addConstr(lhs, GRB_GREATER_EQUAL, curr_cut->rhs());
              addCut(lhs, GRB_GREATER_EQUAL, curr_cut->rhs());
              if (addAsLazy) {
                addCut(lhs, GRB_GREATER_EQUAL, curr_cut->rhs());
              }
            } /* add cuts to the model */
          } /* check that cuts is not null */
        } /* where == MIPNODE (make sure we are in the right "where") */
//        else if (where == GRB_CB_MIP) {
//          // Get the number of times cuts were applied
//          num_total_cuts_applied = GRBCallback::getIntInfo(GRB_CB_MIP_CUTCNT);
//        } /* where == MIP (make sure we are in the right "where") */
//        else if (where == GRB_CB_MESSAGE) {
//          // Number user cuts
//          std::string msg = GRBCallback::getStringInfo(GRB_CB_MSG_STRING);
//        } /* where == MESSAGE */
      } catch (GRBException& e) {
        error_msg(errorstring, "Gurobi: Error during callback: %s\n", e.getMessage().c_str());
        writeErrorToLog(errorstring, params.logfile);
        exit(1);
      } catch (...) {
        error_msg(errorstring, "Gurobi: Error during callback.\n");
        writeErrorToLog(errorstring, params.logfile);
        exit(1);
      }
    } /* callback */
}; /* class GurobiUserCutCallback */

void presolveModelWithGurobi(const VPCParameters& params, int strategy, 
    GRBModel& model, double& presolved_lp_opt, std::string& presolved_name,
    const double best_bound) {
//#ifdef TRACE
  printf("\n## Gurobi: Presolving model ##\n");
//#endif
  try {
    GRBModel presolved_model = model.presolve();
    GRBModel presolved_model_mip = presolved_model;
    GRBVar* vars = presolved_model.getVars();
    const int num_vars = presolved_model.get(GRB_IntAttr_NumVars);
    for (int i = 0; i < num_vars; i++) {
      vars[i].set(GRB_CharAttr_VType, GRB_CONTINUOUS);
    }

    strategy = static_cast<int>(BB_Strategy_Options::presolve_on);
    setStrategyForBBTestGurobi(params, strategy, presolved_model, best_bound);

    presolved_model.optimize();

    // Save optimal value
    int optimstatus = presolved_model.get(GRB_IntAttr_Status);
    if (optimstatus == GRB_OPTIMAL) {
        presolved_lp_opt = presolved_model.get(GRB_DoubleAttr_ObjVal);
        size_t slashindex = presolved_name.find_last_of("/\\");
        presolved_model_mip.set(GRB_StringAttr_ModelName, presolved_name.substr(slashindex+1));
        printf("Saving Gurobi-presolved model to %s.mps.gz.\n", presolved_name.c_str());
        createTmpFileCopy(params, presolved_model_mip, presolved_name, ".mps.gz"); // adds .mps.gz ext
        if (vars) {
          delete[] vars;
        }
    } else {
      error_msg(errorstring, "Gurobi: Was not able to solve presolved model to optimality. Status: %d.\n", optimstatus);
      writeErrorToLog(errorstring, params.logfile);
      exit(1);
    }
  } catch (GRBException& e) {
    error_msg(errorstring, "Gurobi: Exception caught: %s\n", e.getMessage().c_str());
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  } catch (...) {
    error_msg(errorstring, "Gurobi: Unknown exception caught.\n");
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }
} /* presolveModelWithGurobi (GRBModel) */

void presolveModelWithGurobi(const VPCParameters& params, int strategy, const char* f_name,
    double& presolved_lp_opt, std::string& presolved_name, const double best_bound) {
  try {
    GRBEnv env = GRBEnv();
    GRBModel model = GRBModel(env, f_name);
    presolveModelWithGurobi(params, strategy, model, presolved_lp_opt, presolved_name, best_bound);
  } catch (GRBException& e) {
    error_msg(errorstring, "Gurobi: Exception caught: %s\n", e.getMessage().c_str());
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  } catch (...) {
    error_msg(errorstring, "Gurobi: Unknown exception caught.\n");
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }
} /* presolveModelWithGurobi (filename) */

void presolveModelWithGurobi(const VPCParameters& params, int strategy,
    const OsiSolverInterface* const solver, double& presolved_lp_opt,
    std::string& presolved_name, const double best_bound) {
  std::string f_name = "";
  createTmpFileCopy(params, solver, f_name);
  presolveModelWithGurobi(params, strategy, f_name.c_str(), presolved_lp_opt, presolved_name, best_bound);
  remove(f_name.c_str()); // remove temporary file
} /* presolveModelWithGurobi (Osi) */

void doBranchAndBoundWithGurobi(const VPCParameters& params, int strategy,
    GRBModel& model, BBInfo& info, const double best_bound, 
    std::vector<double>* const solution = NULL) {
//#ifdef TRACE
  printf("\n## Running B&B with Gurobi. Strategy: %d. Random seed: %d. ##\n",
      strategy, params.get(intParam::RANDOM_SEED));
//#endif
  try {
    setStrategyForBBTestGurobi(params, strategy, model, best_bound);

    model.optimize();

    int optimstatus = model.get(GRB_IntAttr_Status);

    if (optimstatus == GRB_INF_OR_UNBD) {
      const int presolve_flag = model.get(GRB_IntParam_Presolve);
      if (presolve_flag) {
        model.set(GRB_IntParam_Presolve, 0);
        model.optimize();
        optimstatus = model.get(GRB_IntAttr_Status);
      }
    }

    if (optimstatus == GRB_INFEASIBLE || optimstatus == GRB_UNBOUNDED || optimstatus == GRB_INF_OR_UNBD) {
      error_msg(errorstring, "Gurobi: Failed to optimize MIP.\n");
      writeErrorToLog(errorstring, params.logfile);
      exit(1);
    }

    switch (optimstatus) {
      case GRB_CUTOFF:
      case GRB_ITERATION_LIMIT:
      case GRB_NODE_LIMIT:
      case GRB_TIME_LIMIT: {
        info.obj = model.get(GRB_DoubleAttr_ObjVal);
        info.bound = model.get(GRB_DoubleAttr_ObjBound);
        break;
      }
      case GRB_USER_OBJ_LIMIT: {
        info.obj = model.get(GRB_DoubleAttr_ObjVal);
        info.bound = model.get(GRB_DoubleAttr_ObjBound);
        break;
      }
      case GRB_OPTIMAL: {
        info.obj = model.get(GRB_DoubleAttr_ObjVal);
        info.bound = model.get(GRB_DoubleAttr_ObjBound);
        break;
      }
      default: {
        error_msg(errorstring, "Gurobi: Other status after solve: %d.\n", optimstatus);
        writeErrorToLog(errorstring, params.logfile);
        exit(1);
      }
    } // switch optimistatus
    info.iters = (long) model.get(GRB_DoubleAttr_IterCount);
    info.nodes = (long) model.get(GRB_DoubleAttr_NodeCount);
    info.time = model.get(GRB_DoubleAttr_Runtime);

#ifdef TRACE
    printf("Gurobi: Solution value: %1.6f.\n", info.obj);
    printf("Gurobi: Best bound: %1.6f.\n", info.bound);
    printf("Gurobi: Number iterations: %ld.\n", info.iters);
    printf("Gurobi: Number nodes: %ld.\n", info.nodes);
    printf("Gurobi: Time: %f.\n", info.time);
#endif

   // Save the solution if needed
   if (solution) {
     // Get variables
     GRBVar* vars = model.getVars();
     const int num_vars = model.get(GRB_IntAttr_NumVars);
     (*solution).resize(num_vars);
     for (int i = 0; i < num_vars; i++) {
       (*solution)[i] = vars[i].get(GRB_DoubleAttr_X);
     }
     if (vars) {
       delete[] vars;
     }
   }
  } catch(GRBException& e) {
    error_msg(errorstring, "Gurobi: Exception caught: %s\n", e.getMessage().c_str());
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  } catch (...) {
    error_msg(errorstring, "Gurobi: Unknown exception caught.\n");
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }
} /* doBranchAndBoundWithGurobi (GRBModel) */

void doBranchAndBoundWithUserCutsGurobi(const VPCParameters& params,
    int strategy, GRBModel& model, const OsiCuts* cuts, BBInfo& info,
    const double best_bound, const bool addAsLazy,
    std::vector<double>* const solution = NULL) {
  // Ensure that user cuts setting is enabled
  if (!use_bb_option(strategy, BB_Strategy_Options::user_cuts)) {
    warning_msg(warnstring,
        "Need to use user_cuts option; strategy currently: %d.\n", strategy);
    strategy = enable_bb_option(strategy, BB_Strategy_Options::user_cuts);
  }

  try {
    // Due to using the C++ interface, we need to access the variable list first
    //GRBVar* grb_vars = model.getVars();
    GurobiUserCutCallback cb = GurobiUserCutCallback(params,
        model.get(GRB_IntAttr_NumVars), model.get(GRB_DoubleAttr_ObjCon),
        model.getVars(), cuts, addAsLazy);
    model.setCallback(&cb);

    // Update the model
    model.update();
    /*
    const int retcode = GRBupdatemodel(&model);
    if (retcode) {
      error_msg(errorstring, "Gurobi: Error updating the model; error code: %d.\n", retcode);
      writeErrorToLog(errorstring, params.logfile);
      exit(1);
    }
    */

    // Finally, run B&B
    doBranchAndBoundWithGurobi(params, strategy, model, info, best_bound, solution);

    // Save information
    info.root_passes = cb.info.root_passes;
    if (info.root_passes > 0) {
      info.first_cut_pass = cb.info.first_cut_pass; // second because first is lp opt val
      info.last_cut_pass = cb.info.last_cut_pass;
      info.root_time = cb.info.root_time;
      info.last_sol_time = cb.info.last_sol_time;
    } else {
      info.first_cut_pass = info.obj;
      info.last_cut_pass = info.obj;
      info.root_time = info.time; // all time was spent at the root
      info.last_sol_time = cb.info.last_sol_time; // roughly the same as total time in this case
    }
//#ifdef TRACE
//    printf("Gurobi: Times cuts applied: %d.\n", cb.numTimesApplied);
//#endif
  } catch (GRBException& e) {
    error_msg(errorstring, "Gurobi: Exception caught: %s\n", e.getMessage().c_str());
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  } catch (...) {
    error_msg(errorstring, "Gurobi: Unknown exception caught.\n");
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }
} /* doBranchAndBoundWithUserCutsGurobi (GRBModel) */

void doBranchAndBoundWithGurobi(const VPCParameters& params, int strategy,
    const char* f_name, BBInfo& info, const double best_bound,
    std::vector<double>* const solution) {
#ifdef TRACE
  printf("\n## Reading from file into Gurobi. ##\n");
#endif
  try {
    GRBEnv env = GRBEnv();
    GRBModel model = GRBModel(env, f_name);
    if (use_bb_option(strategy, BB_Strategy_Options::user_cuts)) {
      doBranchAndBoundWithUserCutsGurobi(params, strategy, model, NULL, info, best_bound, false, solution);
    } else {
      doBranchAndBoundWithGurobi(params, strategy, model, info, best_bound, solution);
    }
  } catch (GRBException& e) {
    error_msg(errorstring, "Gurobi: Exception caught: %s\n", e.getMessage().c_str());
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  } catch (...) {
    error_msg(errorstring, "Gurobi: Unknown exception caught.\n");
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }
} /* doBranchAndBoundWithGurobi (filename) */

void doBranchAndBoundWithGurobi(const VPCParameters& params, int strategy,
    const OsiSolverInterface* const solver, BBInfo& info,
    const double best_bound, std::vector<double>* const solution) {
  std::string f_name;
  createTmpFileCopy(params, solver, f_name);
  doBranchAndBoundWithGurobi(params, strategy, f_name.c_str(), info, best_bound,
      solution);
  remove(f_name.c_str()); // remove temporary file
} /* doBranchAndBoundWithGurobi (Osi) */

void doBranchAndBoundWithUserCutsGurobi(const VPCParameters& params,
    int strategy, const char* f_name, const OsiCuts* cuts, BBInfo& info,
    const double best_bound, const bool addAsLazy) {
#ifdef TRACE
  printf("\n## Reading from file into Gurobi and adding user cuts. ##\n");
#endif
  try {
    GRBEnv env = GRBEnv();
    GRBModel model = GRBModel(env, f_name);
    doBranchAndBoundWithUserCutsGurobi(params, strategy, model, cuts, info, best_bound, addAsLazy);
  } catch (GRBException& e) {
    error_msg(errorstring, "Gurobi: Exception caught: %s\n", e.getMessage().c_str());
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  } catch (...) {
    error_msg(errorstring, "Gurobi: Unknown exception caught.\n");
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }
} /* doBranchAndBoundWithUserCutsGurobi (filename) */

void doBranchAndBoundWithUserCutsGurobi(const VPCParameters& params, int strategy,
    const OsiSolverInterface* const solver, const OsiCuts* cuts, BBInfo& info,
    const double best_bound, const bool addAsLazy) {
  std::string f_name;
  createTmpFileCopy(params, solver, f_name);
  doBranchAndBoundWithUserCutsGurobi(params, strategy, f_name.c_str(), cuts, info, best_bound, addAsLazy);
  remove(f_name.c_str()); // remove temporary file
} /* doBranchAndBoundWithUserCutsGurobi (Osi) */
#endif /* USE_GUROBI */
