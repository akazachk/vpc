// Name:     analysis.hpp
// Author:   A. M. Kazachkov
// Date:     2019-Mar-04
//-----------------------------------------------------------------------------
#pragma once

#include <string>
#include <vector>

class OsiCuts;
struct VPCParameters;
struct SummaryBBInfo;

struct RoundInfo {
  double old_obj, new_obj;
  int num_frac_vars;
  double min_frac, max_frac;
  int num_disj_terms;
  int num_obj_tried, num_cuts, num_failures;
  std::vector<int> numFails;
};

void analyzeStrength(const VPCParameters& params, double init_obj,
    double gmic_obj, double vpc_obj, double best_disj_obj, double ip_obj,
    const int num_vpc_total, const int num_gmic_total,
//    OsiCuts* vpcs, OsiCuts* gmics,
    std::string& output);
void analyzeBB(const VPCParameters& params, SummaryBBInfo& info_nocuts,
    SummaryBBInfo& info_mycuts, SummaryBBInfo& info_allcuts, std::string& output);
