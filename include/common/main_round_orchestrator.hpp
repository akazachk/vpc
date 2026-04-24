#pragma once

#include <cstdio>
#include <vector>

#include <OsiCuts.hpp>

#include "analysis.hpp"
#include "main_state.hpp"

namespace main_round_orchestrator {

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
    MainAppState& state);

} // namespace main_round_orchestrator
