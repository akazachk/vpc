#include "main_round_orchestrator.hpp"

#include <cmath>

#include <CglGMI.hpp>

#include "CutHelper.hpp"
#include "main_round_actions.hpp"
#include "utility.hpp"

#ifdef CALC_COND_NUM
#include "condition_number.hpp"
#endif

using namespace VPCParametersNamespace;

namespace main_round_orchestrator {

enum OverallTimeStatsLocal {
  INIT_SOLVE_TIME,
  GOMORY_GEN_TIME,
  GOMORY_APPLY_TIME,
  CUT_TOTAL_TIME,
  VPC_GEN_TIME,
  VPC_APPLY_TIME,
  VPC_INIT_SOLVE_TIME,
  VPC_DISJ_SETUP_TIME,
  VPC_DISJ_GEN_TIME,
  VPC_PRLP_SETUP_TIME,
  VPC_PRLP_SOLVE_TIME,
  VPC_GEN_CUTS_TIME,
  BB_TIME,
  TOTAL_APPLY_TIME,
  TOTAL_TIME,
  NUM_TIME_STATS
};

int runCutRounds(
    int num_rounds,
    bool compute_cond_num,
    FILE* cutrounds_logfile,
    int num_cutround_info,
    std::vector<int>& disjOptions,
    std::vector<OsiCuts>& vpcs_by_round,
    std::vector<SummaryCutInfo>& cutInfoVec,
    std::vector<SummaryCutInfo>& cutInfoGMICVec,
    std::vector<SummaryBoundInfo>& boundInfoVec,
    std::vector<double>& gmic_gen_time_by_round,
    std::vector<double>& gmic_apply_time_by_round,
    std::vector<double>& gmic_cond_num2_by_round,
    std::vector<double>& vpc_gen_time_by_round,
    std::vector<double>& vpc_apply_time_by_round,
    std::vector<double>& vpc_cond_num2_by_round,
    const std::vector<double>& ip_solution,
    OsiCuts& gmics,
    OsiCuts& vpcs,
    MainAppState& state) {
  state.boundInfo.num_vpc = 0;
  state.boundInfo.num_gmic = 0;

  int round_ind = 0;
  for (round_ind = 0; round_ind < num_rounds; ++round_ind) {
    if (num_rounds > 1) {
      printf("\n## Starting round %d/%d. ##\n", round_ind + 1, num_rounds);
    }
    int num_disj = 0;

    boundInfoVec[round_ind].lp_obj = state.solver->getObjValue();
    boundInfoVec[round_ind].num_vpc = 0;
    boundInfoVec[round_ind].num_gmic = 0;

    state.timer.start_timer(OverallTimeStatsLocal::CUT_TOTAL_TIME);

    OsiCuts currGMICs;
    if (state.params.get(GOMORY) == -1 || state.params.get(GOMORY) == 1) {
      const double gmic_gen_start_time = state.timer.get_total_time(OverallTimeStatsLocal::GOMORY_GEN_TIME);
      state.timer.start_timer(OverallTimeStatsLocal::GOMORY_GEN_TIME);

      CglGMI GMIGen;
      GMIGen.getParam().setMAX_SUPPORT(state.params.get(intConst::MIN_SUPPORT_THRESHOLD));
      GMIGen.getParam().setMAX_SUPPORT_REL(state.params.get(doubleParam::MAX_SUPPORT_REL));
      GMIGen.getParam().setMAXDYN(state.solver->getInfinity());
      GMIGen.getParam().setMINVIOL(0.0);
      GMIGen.generateCuts(*state.solver, currGMICs);

      state.timer.end_timer(OverallTimeStatsLocal::GOMORY_GEN_TIME);
      const double gmic_gen_end_time = state.timer.get_total_time(OverallTimeStatsLocal::GOMORY_GEN_TIME);
      gmic_gen_time_by_round[round_ind] = gmic_gen_end_time - gmic_gen_start_time;

      state.boundInfo.num_gmic += currGMICs.sizeCuts();
      boundInfoVec[round_ind].num_gmic += currGMICs.sizeCuts();
      gmics.insert(currGMICs);

      state.timer.start_timer(OverallTimeStatsLocal::GOMORY_APPLY_TIME);
      const double gmic_apply_start_time = state.timer.get_total_time(OverallTimeStatsLocal::GOMORY_APPLY_TIME);
      applyCutsCustom(state.GMICSolver, currGMICs, state.params.logfile);
      const double gmic_apply_end_time = state.timer.get_total_time(OverallTimeStatsLocal::GOMORY_APPLY_TIME);
      state.timer.end_timer(OverallTimeStatsLocal::GOMORY_APPLY_TIME);
      gmic_apply_time_by_round[round_ind] = gmic_apply_end_time - gmic_apply_start_time;

      state.boundInfo.gmic_obj = state.GMICSolver->getObjValue();
      boundInfoVec[round_ind].gmic_obj = state.boundInfo.gmic_obj;

#ifdef CALC_COND_NUM
      if (compute_cond_num) {
        gmic_cond_num2_by_round[round_ind] = compute_condition_number(state.GMICSolver, 2);
      }
#endif

      if (state.params.get(GOMORY) > 0) {
        applyCutsCustom(state.solver, currGMICs, state.params.logfile);
      }
      applyCutsCustom(state.allCutsSolver, currGMICs, state.params.logfile);

      updateGMICInfo(cutInfoGMICVec[round_ind], &currGMICs, state.params.get(EPS) / 2.);
    }

    if ((state.params.get(VPCParametersNamespace::intParam::MODE) != static_cast<int>(VPCParametersNamespace::VPCMode::DISJ_SET_PBB))
        && ((int) disjOptions.size() > round_ind)) {
      state.params.set(intParam::DISJ_TERMS, disjOptions[round_ind]);
    }
    const bool SHOULD_GENERATE_VPCS = (state.params.get(intParam::DISJ_TERMS) != 0 || disjOptions.size() > 0);
    const bool USE_CUSTOM = (state.params.get(intParam::MODE) == static_cast<int>(VPCParametersNamespace::VPCMode::CUSTOM));

    if (SHOULD_GENERATE_VPCS) {
      const double vpc_start_time = state.timer.get_total_time(OverallTimeStatsLocal::VPC_GEN_TIME);
      state.timer.start_timer(OverallTimeStatsLocal::VPC_GEN_TIME);

      CglVPC gen(state.params, round_ind);
      gen.timer.add_value(
          CglVPC::VPCTimeStatsName[static_cast<int>(CglVPC::VPCTimeStats::INIT_SOLVE_TIME)],
          state.timer.get_value(OverallTimeStatsLocal::INIT_SOLVE_TIME));

      if (USE_CUSTOM) {
        main_round_actions::doCustomRoundOfCuts(round_ind, vpcs_by_round[round_ind], gen, num_disj, state, cutInfoVec);
      } else {
        gen.generateCuts(*state.solver, vpcs_by_round[round_ind]);

        state.exitReason = gen.exitReason;
        if (gen.disj() || gen.disjSet()) {
          DisjunctionSet* disjSet = gen.disjSet();
          if (gen.disjSet() == NULL) {
            disjSet = new DisjunctionSet;
            disjSet->addDisjunction(gen.disj());
          }

          state.boundInfo.num_vpc += gen.num_cuts;
          boundInfoVec[round_ind].num_vpc += gen.num_cuts;

          for (const Disjunction* const disj : disjSet->disjunctions) {
            num_disj++;

            if (state.boundInfo.best_disj_obj < disj->best_obj) {
              state.boundInfo.best_disj_obj = disj->best_obj;
            }
            if (boundInfoVec[round_ind].best_disj_obj < disj->best_obj) {
              boundInfoVec[round_ind].best_disj_obj = disj->best_obj;
            }
            if (state.boundInfo.worst_disj_obj < disj->worst_obj) {
              state.boundInfo.worst_disj_obj = disj->worst_obj;
            }
            if (boundInfoVec[round_ind].worst_disj_obj < disj->worst_obj) {
              boundInfoVec[round_ind].worst_disj_obj = disj->worst_obj;
            }

            if (round_ind == 0) {
              state.boundInfo.num_root_bounds_changed = disj->common_changed_var.size();
              state.boundInfo.root_obj = disj->root_obj;
            }
            boundInfoVec[round_ind].num_root_bounds_changed = disj->common_changed_var.size();
            boundInfoVec[round_ind].root_obj = disj->root_obj;
          }

          if (gen.disjSet() == NULL) {
            delete disjSet;
            disjSet = NULL;
          }
        }
        updateDisjInfo(state.disjInfo, num_disj, gen);
        updateCutInfo(cutInfoVec[round_ind], gen, &vpcs_by_round[round_ind], state.params.get(EPS) / 2.);
      }

      const TimeStats& vpctimer = gen.timer;
      std::vector<CglVPC::VPCTimeStats> vpc_stats = {
        CglVPC::VPCTimeStats::INIT_SOLVE_TIME,
        CglVPC::VPCTimeStats::DISJ_SETUP_TIME,
        CglVPC::VPCTimeStats::DISJ_GEN_TIME,
        CglVPC::VPCTimeStats::PRLP_SETUP_TIME,
        CglVPC::VPCTimeStats::PRLP_SOLVE_TIME,
        CglVPC::VPCTimeStats::GEN_CUTS_TIME,
      };
      std::vector<OverallTimeStatsLocal> overall_stats = {
        OverallTimeStatsLocal::VPC_INIT_SOLVE_TIME,
        OverallTimeStatsLocal::VPC_DISJ_SETUP_TIME,
        OverallTimeStatsLocal::VPC_DISJ_GEN_TIME,
        OverallTimeStatsLocal::VPC_PRLP_SETUP_TIME,
        OverallTimeStatsLocal::VPC_PRLP_SOLVE_TIME,
        OverallTimeStatsLocal::VPC_GEN_CUTS_TIME
      };
      for (int i = 0; i < (int) vpc_stats.size(); i++) {
        const CglVPC::VPCTimeStats stat = vpc_stats[i];
        const OverallTimeStatsLocal overall_stat = overall_stats[i];
        const clock_t currtimevalue = vpctimer.get_value(CglVPC::VPCTimeStatsName[static_cast<int>(stat)]);
        state.timer.add_value(overall_stat, currtimevalue);
      }

      state.timer.end_timer(OverallTimeStatsLocal::VPC_GEN_TIME);
      const double vpc_end_time = state.timer.get_total_time(OverallTimeStatsLocal::VPC_GEN_TIME);
      vpc_gen_time_by_round[round_ind] = vpc_end_time - vpc_start_time;
    } else {
      setCutInfo(cutInfoVec[round_ind], 0, NULL);
    }

    if (state.params.get(TEMP) == static_cast<int>(TempOptions::CHECK_CUTS_AGAINST_BB_OPT) && !ip_solution.empty()) {
      const int num_violated = checkCutsAgainstFeasibleSolution(vpcs_by_round[round_ind], ip_solution);
      printf("\n## Number of cuts violating the IP solution: %d ##\n", num_violated);
    }

    vpcs.insert(vpcs_by_round[round_ind]);

    state.timer.end_timer(OverallTimeStatsLocal::CUT_TOTAL_TIME);

    printf("\n## Applying disjunctive cuts (# cuts = %d). ##\n", (int) vpcs_by_round[round_ind].sizeCuts());
    state.timer.start_timer(OverallTimeStatsLocal::TOTAL_APPLY_TIME);

    if (vpcs_by_round[round_ind].sizeCuts() > 0) {
      if (state.params.get(GOMORY) != 0) {
        const double vpc_apply_start_time = state.timer.get_total_time(OverallTimeStatsLocal::VPC_APPLY_TIME);
        state.timer.start_timer(OverallTimeStatsLocal::VPC_APPLY_TIME);
        applyCutsCustom(state.VPCSolver, vpcs_by_round[round_ind], state.params.logfile);
        state.timer.end_timer(OverallTimeStatsLocal::VPC_APPLY_TIME);
        const double vpc_apply_end_time = state.timer.get_total_time(OverallTimeStatsLocal::VPC_APPLY_TIME);
        vpc_apply_time_by_round[round_ind] = vpc_apply_end_time - vpc_apply_start_time;

        applyCutsCustom(state.solver, vpcs_by_round[round_ind], state.params.logfile);
        applyCutsCustom(state.allCutsSolver, vpcs_by_round[round_ind], state.params.logfile);
      } else {
        const double vpc_apply_start_time = state.timer.get_total_time(OverallTimeStatsLocal::VPC_APPLY_TIME);
        state.timer.start_timer(OverallTimeStatsLocal::VPC_APPLY_TIME);
        applyCutsCustom(state.solver, vpcs_by_round[round_ind], state.params.logfile);
        state.timer.end_timer(OverallTimeStatsLocal::VPC_APPLY_TIME);
        const double vpc_apply_end_time = state.timer.get_total_time(OverallTimeStatsLocal::VPC_APPLY_TIME);
        vpc_apply_time_by_round[round_ind] = vpc_apply_end_time - vpc_apply_start_time;
      }
    }
    state.timer.end_timer(OverallTimeStatsLocal::TOTAL_APPLY_TIME);

    state.boundInfo.vpc_obj      = (state.params.get(GOMORY) != 0) ? state.VPCSolver->getObjValue() : state.solver->getObjValue();
    state.boundInfo.gmic_vpc_obj = (state.params.get(GOMORY) != 0) ? state.allCutsSolver->getObjValue() : state.solver->getObjValue();
    state.boundInfo.all_cuts_obj = state.boundInfo.gmic_vpc_obj;

    boundInfoVec[round_ind].gmic_vpc_obj = state.boundInfo.gmic_vpc_obj;
    boundInfoVec[round_ind].vpc_obj      = state.boundInfo.vpc_obj;
    boundInfoVec[round_ind].all_cuts_obj = state.boundInfo.all_cuts_obj;

#ifdef CALC_COND_NUM
    if (compute_cond_num) {
      vpc_cond_num2_by_round[round_ind] = compute_condition_number((state.params.get(GOMORY) != 0) ? state.VPCSolver : state.solver, 2);
    }
#endif

    printf("\n*** INFO:\n");
    printf(
        "Round %d/%d: Completed round of VPC generation (exit reason: %s).\n",
        round_ind + 1, state.params.get(ROUNDS),
        CglVPC::ExitReasonName[static_cast<int>(state.exitReason)].c_str());
    printf(
        "VPCs generated = %d. GMICs generated = %d.\n",
        vpcs_by_round[round_ind].sizeCuts(),
        boundInfoVec[round_ind].num_gmic);
    printf("Obj value: %s.\n",
        stringValue(state.solver->getObjValue(), "%1.6f").c_str());
    printf("Initial obj value: %s.\n",
        stringValue(state.boundInfo.lp_obj, "%1.6f").c_str());
    if (round_ind > 0) {
      printf("Round start obj value: %s.\n",
          stringValue(boundInfoVec[round_ind - 1].all_cuts_obj, "%1.6f").c_str());
    }
    printf("Disj lb: %s.\n",
        stringValue(state.boundInfo.best_disj_obj, "%1.6f").c_str());
    printf("Elapsed time: %1.2f seconds. Time limit = %1.2f seconds.\n",
        state.timer.get_total_time(state.overallTimeStatsName[OverallTimeStatsLocal::TOTAL_TIME]),
        state.params.get(TIMELIMIT));
    printf("***\n");
    fflush(stdout);

    if (use_temp_option(state.params.get(intParam::TEMP), TempOptions::PRINT_BOUND_BY_ROUND)) {
      int count = 0;
      fprintf(cutrounds_logfile, "%d,", round_ind + 1); count++;
      fprintf(cutrounds_logfile, "%s,", stringValue(boundInfoVec[round_ind].gmic_obj, "%.20f").c_str()); count++;
      fprintf(cutrounds_logfile, "%s,", stringValue(boundInfoVec[round_ind].vpc_obj, "%.20f").c_str()); count++;
      fprintf(cutrounds_logfile, "%s,", stringValue(boundInfoVec[round_ind].all_cuts_obj, "%.20f").c_str()); count++;
      fprintf(cutrounds_logfile, "%d,", boundInfoVec[round_ind].num_gmic); count++;
      fprintf(cutrounds_logfile, "%s,",
        stringValue(
          (boundInfoVec[round_ind].num_gmic > 0)
            ? 100. * cutInfoGMICVec[round_ind].min_support / state.solver->getNumCols()
            : 0.,
          "%.6f",
          1e100
        ).c_str()
      ); count++;
      fprintf(cutrounds_logfile, "%s,",
        stringValue(
          (boundInfoVec[round_ind].num_gmic > 0)
            ? 100. * cutInfoGMICVec[round_ind].max_support / state.solver->getNumCols()
            : 0.,
          "%.6f",
          1e100
        ).c_str()
      ); count++;
      fprintf(cutrounds_logfile, "%s,",
        stringValue(
          (boundInfoVec[round_ind].num_gmic > 0)
            ? 100. * cutInfoGMICVec[round_ind].avg_support / state.solver->getNumCols()
            : 0.,
          "%.6f",
          1e100
        ).c_str()
      ); count++;
      fprintf(cutrounds_logfile, "%.20f,", gmic_gen_time_by_round[round_ind]); count++;
      fprintf(cutrounds_logfile, "%.20f,", gmic_apply_time_by_round[round_ind]); count++;
#ifdef CALC_COND_NUM
      if (compute_cond_num) { fprintf(cutrounds_logfile, "%s,", stringValue(gmic_cond_num2_by_round[round_ind], "%.20f").c_str()); count++; }
#endif
      fprintf(cutrounds_logfile, "%d,", boundInfoVec[round_ind].num_vpc); count++;
      fprintf(cutrounds_logfile, "%s,",
        stringValue(
          (boundInfoVec[round_ind].num_vpc > 0)
            ? 100. * cutInfoVec[round_ind].min_support / state.solver->getNumCols()
            : 0.,
          "%.6f",
          1e100
        ).c_str()
      ); count++;
      fprintf(cutrounds_logfile, "%s,",
        stringValue(
          (boundInfoVec[round_ind].num_vpc > 0)
            ? 100. * cutInfoVec[round_ind].max_support / state.solver->getNumCols()
            : 0.,
          "%.6f",
          1e100
        ).c_str()
      ); count++;
      fprintf(cutrounds_logfile, "%s,",
        stringValue(
          (boundInfoVec[round_ind].num_vpc > 0)
            ? 100. * cutInfoVec[round_ind].avg_support / state.solver->getNumCols()
            : 0.,
          "%.6f",
          1e100
        ).c_str()
      ); count++;
      fprintf(cutrounds_logfile, "%.20f,", vpc_gen_time_by_round[round_ind]); count++;
      fprintf(cutrounds_logfile, "%.20f,", vpc_apply_time_by_round[round_ind]); count++;
#ifdef CALC_COND_NUM
      if (compute_cond_num) { fprintf(cutrounds_logfile, "%s,", stringValue(vpc_cond_num2_by_round[round_ind], "%.20f").c_str()); count++; }
#endif
      fprintf(cutrounds_logfile, "\n");
      assert(count == num_cutround_info);
    }

    if (state.timer.reachedTimeLimit(state.overallTimeStatsName[OverallTimeStatsLocal::TOTAL_TIME], state.params.get(TIMELIMIT))) {
      printf("\n*** INFO: Reached time limit (current time = %f, time limit = %f).\n",
          state.timer.get_total_time(state.overallTimeStatsName[OverallTimeStatsLocal::TOTAL_TIME]), state.params.get(TIMELIMIT));
      break;
    }

    if ((SHOULD_GENERATE_VPCS && vpcs_by_round[round_ind].sizeCuts() == 0)
        || (!SHOULD_GENERATE_VPCS && boundInfoVec[round_ind].num_gmic == 0)
        || !state.solver->isProvenOptimal()
        || isInfinity(std::abs(state.solver->getObjValue()))
        || state.exitReason == CglVPC::ExitReason::OPTIMAL_SOLUTION_FOUND_EXIT) {
      break;
    }
  }

  if (round_ind < num_rounds) {
    num_rounds = round_ind + 1;
  }

  return num_rounds;
}

} // namespace main_round_orchestrator
