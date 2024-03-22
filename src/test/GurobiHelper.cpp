/**
 * @file GurobiHelper.cpp
 * @author A. M. Kazachkov
 * @date 2019-Feb-28
 */
#include "GurobiHelper.hpp"

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
#include <gurobi_c++.h>

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
  model.set(GRB_DoubleParam_TimeLimit, params.get(doubleParam::BB_TIMELIMIT)); // time limit
  model.set(GRB_IntParam_Threads, 1); // single-threaded
  if (seed >= 0) {
    model.set(GRB_IntParam_Seed, seed); // random seed
  }

  if (params.get(VERBOSITY) == 0) {
    model.set(GRB_IntParam_OutputFlag, 0); // turn off output
  } else {
    model.set(GRB_IntParam_OutputFlag, 1); // turn on output
    model.set(GRB_IntParam_DisplayInterval, 1); // set display frequency
  }

  if (strategy <= 0) {
    // Default strategy
  } else {
    if (use_bb_option(strategy, BB_Strategy_Options::user_cuts)) {
      //model.set(GRB_IntParam_DualReductions, 0); // disable dual reductions
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
        model.set(GRB_DoubleParam_BestBdStop, best_bound * (1 - 1e-5)); // give the solver the best dual bound (it is a minimization problem) with a tolerance
        //model.set(GRB_DoubleParam_Cutoff, best_bound + 1e-3); // give the solver the best IP objective value (it is a minimization problem) with a tolerance
      }
      // Check if user provides mip start or solution file
      std::string solfile = params.get(stringParam::SOLFILE);
      std::string ext1 = "_gurobi.sol.gz";
      std::string ext2 = "_gurobi.sol";
      std::string ext3 = "_gurobi.mst.gz";
      std::string ext4 = "_gurobi.mst";
      bool user_provides_start = false;
      user_provides_start |= (solfile.size() > ext1.size()) && (solfile.compare(solfile.size() - ext1.size(), ext1.size(), ext1) == 0);
      user_provides_start |= (solfile.size() > ext2.size()) && (solfile.compare(solfile.size() - ext2.size(), ext2.size(), ext2) == 0);
      user_provides_start |= (solfile.size() > ext3.size()) && (solfile.compare(solfile.size() - ext3.size(), ext3.size(), ext3) == 0);
      user_provides_start |= (solfile.size() > ext4.size()) && (solfile.compare(solfile.size() - ext4.size(), ext4.size(), ext4) == 0);
      if (user_provides_start) {
        model.read(solfile);
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
    bool setFirstCutPass;
    double first_lp_opt;
    BBInfo info;

    GurobiUserCutCallback(const VPCParameters& params, int num_vars, double obj_offset,
        GRBVar* vars, const OsiCuts* cutsToAdd, const bool lazyFlag = false)
        : params(params), num_vars(num_vars), obj_offset(obj_offset), grb_vars(vars),
        cuts(cutsToAdd), addAsLazy(lazyFlag) {
      this->info.obj = params.get(doubleParam::INF);
      this->first_lp_opt = -1 * params.get(doubleParam::INF);
      this->setFirstCutPass = false;
      this->info.root_passes = 0;
      this->info.first_cut_pass = -1 * params.get(doubleParam::INF);
      this->info.last_cut_pass = -1 * params.get(doubleParam::INF);
      this->info.root_iters = 0;
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
        }
        else if (where == GRB_CB_MIP) {
          const int num_nodes = GRBCallback::getDoubleInfo(GRB_CB_MIP_NODCNT);
          if (num_nodes > 0) {
            return;
          }
          this->info.root_iters = GRBCallback::getDoubleInfo(GRB_CB_MIP_ITRCNT);
          if (this->info.root_passes == 0) {
            // In case we never enter the root node, keep the first lp opt here
            // This should be updated at the end of presolve?
            // Example in which this is needed: neos-796608_presolved
            const double objValue = GRBCallback::getDoubleInfo(GRB_CB_MIP_OBJBND);
            this->first_lp_opt = objValue;
            this->info.first_cut_pass = objValue;
            this->info.last_cut_pass = objValue;
          }
          else if (this->info.root_passes == 1 && !this->setFirstCutPass) {
            // If we have exactly one round of cuts, the very next MIP cb after the MIPNODE cb
            // should contain the value after that first round of cuts was added
            // However, the MIP cb may be reached several times after the MIPNODE cb
            // (for example, if heuristics are called)
            // For simplicity, we will store the value of the last time the MIP cb is entered with root_passes == 1
            // NB: the obj bound from this call is *rounded* to account for integrality
            const double objValue = GRBCallback::getDoubleInfo(GRB_CB_MIP_OBJBND);
            this->info.first_cut_pass = objValue;
            this->info.last_cut_pass = objValue;
            this->setFirstCutPass = true;

            //const int cutCount = GRBCallback::getIntInfo(GRB_CB_MIP_CUTCNT);
            //printf("### DEBUG: MIP: (pass %ld) Obj value = %f.\n", this->info.root_passes, objValue);
            //printf("### DEBUG: MIP: (pass %ld) Cut count = %d.\n", this->info.root_passes, cutCount);
          }
        }
        else if (where == GRB_CB_MIPNODE) {
          // Note that the MIPNODE callback will be called once for each cut pass during the root node solve.
          // The MIPNODE_NODCNT value will remain at 0 until the root node is complete.
          // If you query relaxation values from during the root node,
          // the first MIPNODE callback will give the relaxation with no cutting planes,
          // and the last will give the relaxation after all root cuts have been applied.
          // ISSUE: if there is exactly one round of cuts, how do we get the relaxation value after that round?
          const int num_nodes = GRBCallback::getDoubleInfo(GRB_CB_MIPNODE_NODCNT);
          if (num_nodes > 0) {
            return;
          }
          this->info.root_passes++;
          //printf("### DEBUG: MIPNODE: Cut pass %ld.\n", this->info.root_passes);
          //printf("### DEBUG: MIPNODE: Obj value = %f.\n", getObjValue());
          this->info.root_time = GRBCallback::getDoubleInfo(GRB_CB_RUNTIME);

          if (GRBCallback::getIntInfo(GRB_CB_MIPNODE_STATUS) != GRB_OPTIMAL) {
            return;
          }

          const double objValue = getObjValue();
          this->info.last_cut_pass = objValue;

#ifdef SAVE_CUT_ROUNDS
          {
            printf("*** INFO: Root pass %03ld: Obj value = %.6f\n", this->info.root_passes, objValue);
            const std::string logname = params.get(stringParam::LOGFILE);
            std::string fileWithCuts = "cutrounds.csv";
            std::string dir, instname, in_file_ext;
            if (!logname.empty()) {
              parseFilename(dir, instname, in_file_ext, logname, params.logfile);
              fileWithCuts = dir + "/" + instname + "_" + fileWithCuts;
            }
            FILE* cutrounds = fopen(fileWithCuts.c_str(), "a");
            if (this->info.root_passes == 1) {
              parseFilename(dir, instname, in_file_ext, params.get(stringParam::FILENAME), params.logfile);
              fprintf(cutrounds, "\n%s,", instname.c_str());
            }
            fprintf(cutrounds, "%.20f,", objValue);
            fclose(cutrounds);
          }
#endif

          // Make sure this is our first entry into the root
          if (this->info.root_passes > 1) {
            if (this->info.root_passes == 2) {
              this->info.first_cut_pass = objValue;
              this->setFirstCutPass = true;
            }
            return;
          }

          this->first_lp_opt = objValue; // root LP relaxation
          this->info.first_cut_pass = objValue; // stores LP relaxation value because we only get here when root_passes == 1

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
        // ignore warnings that a cut was skipped - that's not critical
        if (e.getMessage() != "addCut"){
          error_msg(errorstring, "Gurobi: Error during callback: %s\n", e.getMessage().c_str());
          writeErrorToLog(errorstring, params.logfile);
          exit(1);
        }
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

    printf("Solving Gurobi-presolved model to get new LP optimum value.\n");
    presolved_model.optimize();

    // Save optimal value
    int optimstatus = presolved_model.get(GRB_IntAttr_Status);
    if (optimstatus == GRB_OPTIMAL) {
        presolved_lp_opt = presolved_model.get(GRB_DoubleAttr_ObjVal);
        size_t slashindex = presolved_name.find_last_of("/\\");
        presolved_model_mip.set(GRB_StringAttr_ModelName, presolved_name.substr(slashindex+1));
        printf("Saving Gurobi-presolved model to \"%s.mps.gz\".\n", presolved_name.c_str());
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
  printf("\tRunning B&B with Gurobi. Strategy: %d. Random seed: %d.\n",
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
    printf("Gurobi: Solution value: %s.\n", stringValue(info.obj, "%1.6f", GRB_INFINITY).c_str());
    printf("Gurobi: Best bound: %s.\n", stringValue(info.bound, "%1.6f", GRB_INFINITY).c_str());
    printf("Gurobi: Number iterations: %s.\n", stringValue(info.iters, "%ld").c_str());
    printf("Gurobi: Number nodes: %s.\n", stringValue(info.nodes, "%ld").c_str());
    printf("Gurobi: Time: %f.\n", info.time);
#endif

   // Save the solution if needed
   if (solution && model.get(GRB_IntAttr_SolCount) > 0) {
     // Get variables and double-check objective matches
     double obj = model.get(GRB_DoubleAttr_ObjCon);
     GRBVar* vars = model.getVars();
     const int num_vars = model.get(GRB_IntAttr_NumVars);
     (*solution).resize(num_vars);
     for (int i = 0; i < num_vars; i++) {
       (*solution)[i] = vars[i].get(GRB_DoubleAttr_X); // will throw an error if no feasible solution is available
       obj += (*solution)[i] * vars[i].get(GRB_DoubleAttr_Obj);
     }
     if (vars) {
       delete[] vars;
     }
     if (!isVal(obj, info.obj)) {
       error_msg(errorstring, "Gurobi: %.6e, computed objective value from solution, does not match solver's obj value %.6e.\n", obj, info.obj);
       writeErrorToLog(errorstring, params.logfile);
       exit(1);
     }
   } // save ip solution

   if (use_temp_option(params.get(intParam::TEMP), TempOptions::SAVE_IP_OPT) && model.get(GRB_IntAttr_SolCount) > 0) {
     std::string dir, instname, in_file_ext;
     std::string filename = params.get(stringParam::FILENAME);
     parseFilename(dir, instname, in_file_ext, filename, params.logfile);
     std::string f_name = dir + "/" + instname + "_gurobi.mst.gz";
     if (!fexists(f_name.c_str())) {
       printf("Saving Gurobi MIP start file to \"%s\".\n", f_name.c_str());
       model.write(f_name);
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
      info.root_iters = cb.info.root_iters;
      info.root_time = cb.info.root_time;
      info.last_sol_time = cb.info.last_sol_time;
    } else {
      // I guess this can happen if we solve during presolve,
      // or if we do no root passes of cuts
      if (info.nodes == 0) {
        info.first_cut_pass = info.obj;
        info.last_cut_pass = info.obj;
        info.root_iters = info.iters; // all iters were at the root
        info.root_time = info.time; // all time was spent at the root
        info.last_sol_time = cb.info.last_sol_time; // roughly the same as total time in this case
      } else {
        info.first_cut_pass = cb.first_lp_opt;
        info.last_cut_pass = cb.first_lp_opt;
        info.root_iters = cb.info.root_iters;
        info.root_time = cb.info.root_time; // all time was spent at the root
        info.last_sol_time = cb.info.last_sol_time; // roughly the same as total time in this case
      }
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
  // create a new strings to match what was actually created
  std::string f_name_no_ext = f_name.substr(0, f_name.size() - 4);
  std::string f_name_gz = f_name + ".gz";
  doBranchAndBoundWithUserCutsGurobi(params, strategy, f_name.c_str(), cuts, info, best_bound, addAsLazy);
  // remove temporary files
  remove(f_name.c_str());
  remove(f_name_gz.c_str());
  remove(f_name_no_ext.c_str());
} /* doBranchAndBoundWithUserCutsGurobi (Osi) */
#endif /* USE_GUROBI */
