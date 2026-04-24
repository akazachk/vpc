#include "main_round_actions.hpp"

#include <cmath>
#include <cstdlib>

#include "CutHelper.hpp"
#include "DisjunctionHelper.hpp"
#include "utility.hpp"

namespace main_round_actions {

void doCustomRoundOfCuts(
    int round_ind,
    OsiCuts& vpcs,
    CglVPC& gen,
    int& num_disj,
    MainAppState& state,
    std::vector<SummaryCutInfo>& cutInfoVec) {
  printf("\n## Setting up disjunction(s) ##\n");

  gen.timer.start_timer(CglVPC::VPCTimeStatsName[static_cast<int>(CglVPC::VPCTimeStats::DISJ_GEN_TIME)]);
  DisjunctionSet* disjSet = new DisjunctionSet();
  CglVPC::ExitReason setDisjExitReason = setDisjunctions(disjSet, state.solver, state.params, VPCParametersNamespace::VPCMode::SPLITS);
  gen.timer.end_timer(CglVPC::VPCTimeStatsName[static_cast<int>(CglVPC::VPCTimeStats::DISJ_GEN_TIME)]);

  const int numDisj = disjSet->size();

  if (setDisjExitReason == CglVPC::ExitReason::OPTIMAL_SOLUTION_FOUND_EXIT) {
    warning_msg(warnstr,
        "An integer (optimal) solution was found prior while getting disjunction. "
        "We will generate between n and 2n cuts, restricting the value of each variable.\n");
    const double* solution = disjSet->integer_sol.data();
    if (!solution) {
      error_msg(errorstring,
          "Though status is that optimal integer solution found, unable to find this solution.\n");
      writeErrorToLog(errorstring, state.params.logfile);
      exit(1);
    }
    for (int col = 0; col < state.solver->getNumCols(); col++) {
      const double val = solution[col];

      for (int b = 0; b < 2; b++) {
        if ((b == 0 && greaterThanVal(val, state.solver->getColLower()[col]))
            || (b == 1 && lessThanVal(val, state.solver->getColUpper()[col]))) {
          const double mult = (b == 0) ? 1. : -1.;
          const double el = mult * 1.;

          OsiRowCut currCut;
          currCut.setLb(mult * val);
          currCut.setRow(1, &col, &el, false);
          gen.addCut(currCut, vpcs,
              CglVPC::CutType::OPTIMALITY_CUT,
              CglVPC::ObjectiveType::ONE_SIDED,
              0);
        }
      }
    }
    state.exitReason = CglVPC::ExitReason::OPTIMAL_SOLUTION_FOUND_EXIT;
    state.boundInfo.num_vpc += gen.num_cuts;
    updateCutInfo(cutInfoVec[round_ind], gen, &vpcs, state.params.get(VPCParametersNamespace::EPS) / 2.);
  } else if (setDisjExitReason == CglVPC::ExitReason::SUCCESS_EXIT && numDisj > 0) {
    const int cutLimit = std::ceil(
        gen.getCutLimit(state.params.get(VPCParametersNamespace::CUTLIMIT),
            state.solver->getFractionalIndices().size()) / numDisj);
    gen.setupRepeatedUse(true);
    for (Disjunction* disj : disjSet->disjunctions) {
      if (!disj)
        continue;
#ifdef TRACE
      printf("\n## Generating cuts from disj %s ##\n", disj->name.c_str());
#endif
      num_disj++;
      gen.params.set(VPCParametersNamespace::CUTLIMIT, cutLimit);
      gen.setDisjunction(disj, false);
      gen.generateCuts(*state.solver, vpcs);
      state.exitReason = gen.exitReason;
      state.boundInfo.num_vpc += gen.num_cuts;
      if (state.boundInfo.best_disj_obj < gen.disj()->best_obj)
        state.boundInfo.best_disj_obj = gen.disj()->best_obj;
      if (state.boundInfo.worst_disj_obj < gen.disj()->worst_obj)
        state.boundInfo.worst_disj_obj = gen.disj()->worst_obj;
      updateDisjInfo(state.disjInfo, num_disj, gen);
      updateCutInfo(cutInfoVec[round_ind], gen, &vpcs, state.params.get(VPCParametersNamespace::EPS) / 2.);
    }
    gen.setupRepeatedUse(false);
  } else if (setDisjExitReason == CglVPC::ExitReason::NO_DISJUNCTION_EXIT) {
    state.exitReason = CglVPC::ExitReason::NO_DISJUNCTION_EXIT;
  } else {
    error_msg(errorstring, "Unknown exit reason (%s) from setDisjunctions.\n", CglVPC::ExitReasonName[static_cast<int>(setDisjExitReason)].c_str());
    writeErrorToLog(errorstring, state.params.logfile);
    exit(1);
  }

  if (disjSet) {
    delete disjSet;
  }
}

} // namespace main_round_actions
