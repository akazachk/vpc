#pragma once

#include <vector>

#include <OsiCuts.hpp>

#include "CglVPC.hpp"
#include "main_state.hpp"

namespace main_round_actions {

void doCustomRoundOfCuts(
    int round_ind,
    OsiCuts& vpcs,
    CglVPC& gen,
    int& num_disj,
    MainAppState& state,
    std::vector<SummaryCutInfo>& cutInfoVec);

} // namespace main_round_actions
