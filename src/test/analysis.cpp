// Name:     analysis.cpp
// Author:   A. M. Kazachkov
// Date:     2019-Mar-04
//-----------------------------------------------------------------------------
#include "analysis.hpp"

// COIN-OR
#include <OsiCuts.hpp>

// Project files
#include "BBHelper.hpp"
#include "VPCParameters.hpp"
#include "utility.hpp" // isInfinity, stringValue

void printHeader(FILE* logfile) {

} /* printHeader */

/**
 * The cut properties we want to look at are:
 * 1. Gap closed
 * 2. Activity (after adding cuts)
 * 3. Density
 */
void analyzeStrength(const VPCParameters& params, double init_obj,
    double gmic_obj, double vpc_obj, double best_disj_obj, double ip_obj,
    const int num_vpc_total, const int num_gmic_total,
//    OsiCuts* vpcs, OsiCuts* gmics,
    std::string& output) {
//  const int num_vpc_total = (vpcs != NULL) ? vpcs->sizeCuts() : 0;
//  const int num_gmic_total = (gmics != NULL) ? gmics->sizeCuts() : 0;

  // Print results from adding cuts
  int NAME_WIDTH = 25;
  int NUM_DIGITS_BEFORE_DEC = 7;
  int NUM_DIGITS_AFTER_DEC = 7;
  char tmpstring[300];

  snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "\n## Results from adding cuts ##\n");
  output += tmpstring;
  snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%-*.*s%s\n", NAME_WIDTH, NAME_WIDTH, "LP: ",
      stringValue(init_obj, "% -*.*f", NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC).c_str());
  output += tmpstring;
  if (!isInfinity(std::abs(gmic_obj))) {
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%-*.*s%s (%d cuts)\n", NAME_WIDTH, NAME_WIDTH, "GMICs: ",
        stringValue(gmic_obj, "% -*.*f", NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC).c_str(),
        num_gmic_total);
    output += tmpstring;
  }
  if (!isInfinity(std::abs(vpc_obj))) {
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%-*.*s%s (%d cuts)\n", NAME_WIDTH, NAME_WIDTH, "VPCs: ",
        stringValue(vpc_obj, "% -*.*f", NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC).c_str(),
        num_vpc_total);
    output += tmpstring;
  }
  if (!isInfinity(std::abs(best_disj_obj))) {
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%-*.*s%s\n", NAME_WIDTH, NAME_WIDTH, "Disjunctive lb: ",
        stringValue(best_disj_obj, "% -*.*f", NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC).c_str());
    output += tmpstring;
  }
  if (!isInfinity(std::abs(ip_obj))) {
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%-*.*s%s\n", NAME_WIDTH, NAME_WIDTH, "IP: ",
        stringValue(ip_obj, "% -*.*f", NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC).c_str());
    output += tmpstring;
  }
} /* analyzeStrength */

void analyzeBB(const VPCParameters& params, SummaryBBInfo& info_nocuts,
    SummaryBBInfo& info_mycuts, SummaryBBInfo& info_allcuts, std::string& output) {
  if (params.get(BB_RUNS) == 0) {
    return;
  }
  // B&B mode: ones bit = no_cuts, tens bit = w/vpcs, hundreds bit = w/gmics
  const int mode_param = params.get(intParam::BB_MODE);
  const int mode_ones = mode_param % 10;
  const int mode_tens = (mode_param % 100 - (mode_param % 10)) / 10;
  const int mode_hundreds = (mode_param % 1000 - (mode_param % 100)) / 100;
  const bool branch_with_no_cuts = (mode_ones > 0);
  const bool branch_with_vpcs = (mode_tens > 0) && (info_mycuts.num_cuts > 0);
  const bool branch_with_gmics = (mode_hundreds > 0) && (info_allcuts.num_cuts > 0);
  if (branch_with_no_cuts + branch_with_vpcs + branch_with_gmics == 0) {
    return;
  }

  // Save results to string and also print to the logfile
  printBBInfo(info_mycuts.avg_bb_info, info_allcuts.avg_bb_info, params.logfile);
  int NAME_WIDTH = 10; //25
  int NUM_DIGITS_BEFORE_DEC = 15; //10
  int NUM_DIGITS_AFTER_DEC = 2; //2
  char tmpstring[300];

  snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "\n## Branch-and-bound results ##\n"); output += tmpstring;
  snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%-*s", NAME_WIDTH, ""); output += tmpstring;
  snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%-*s", NUM_DIGITS_BEFORE_DEC, "Obj"); output += tmpstring;
  snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%-*s", NUM_DIGITS_BEFORE_DEC, "Bound"); output += tmpstring;
  snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%-*s", NUM_DIGITS_BEFORE_DEC, "Iters"); output += tmpstring;
  snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%-*s", NUM_DIGITS_BEFORE_DEC, "Nodes"); output += tmpstring;
  snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%-*s", NUM_DIGITS_BEFORE_DEC, "Root passes"); output += tmpstring;
  snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%-*s", NUM_DIGITS_BEFORE_DEC, "First cut pass"); output += tmpstring;
  snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%-*s", NUM_DIGITS_BEFORE_DEC, "Last cut pass"); output += tmpstring;
  snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%-*s", NUM_DIGITS_BEFORE_DEC, "Root time"); output += tmpstring;
  snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%-*s", NUM_DIGITS_BEFORE_DEC, "Last sol time"); output += tmpstring;
  snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%-*s", NUM_DIGITS_BEFORE_DEC, "Time"); output += tmpstring;
  snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "\n"); output += tmpstring;
  if (branch_with_no_cuts) {
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%-*s", NAME_WIDTH, "Gur"); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%s", stringValue(info_nocuts.avg_bb_info.obj, "%-*.*f", NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC).c_str()); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%s", stringValue(info_nocuts.avg_bb_info.bound, "%-*.*f", NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC).c_str()); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%-*ld", NUM_DIGITS_BEFORE_DEC, info_nocuts.avg_bb_info.iters); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%-*ld", NUM_DIGITS_BEFORE_DEC, info_nocuts.avg_bb_info.nodes); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%-*ld", NUM_DIGITS_BEFORE_DEC, info_nocuts.avg_bb_info.root_passes); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%s", stringValue(info_nocuts.avg_bb_info.first_cut_pass, "%-*.*f", NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC).c_str()); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%s", stringValue(info_nocuts.avg_bb_info.last_cut_pass, "%-*.*f", NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC).c_str()); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%s", stringValue(info_nocuts.avg_bb_info.root_time, "%-*.*f", NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC).c_str()); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%s", stringValue(info_nocuts.avg_bb_info.last_sol_time, "%-*.*f", NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC).c_str()); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%s", stringValue(info_nocuts.avg_bb_info.time, "%-*.*f", NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC).c_str()); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "\n"); output += tmpstring;
  } // no cuts
  if (branch_with_vpcs) {
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%-*s", NAME_WIDTH, "Gur+V"); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%s", stringValue(info_mycuts.avg_bb_info.obj, "%-*.*f", NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC).c_str()); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%s", stringValue(info_mycuts.avg_bb_info.bound, "%-*.*f", NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC).c_str()); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%-*ld", NUM_DIGITS_BEFORE_DEC, info_mycuts.avg_bb_info.iters); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%-*ld", NUM_DIGITS_BEFORE_DEC, info_mycuts.avg_bb_info.nodes); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%-*ld", NUM_DIGITS_BEFORE_DEC, info_mycuts.avg_bb_info.root_passes); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%s", stringValue(info_mycuts.avg_bb_info.first_cut_pass, "%-*.*f", NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC).c_str()); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%s", stringValue(info_mycuts.avg_bb_info.last_cut_pass, "%-*.*f", NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC).c_str()); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%s", stringValue(info_mycuts.avg_bb_info.root_time, "%-*.*f", NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC).c_str()); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%s", stringValue(info_mycuts.avg_bb_info.last_sol_time, "%-*.*f", NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC).c_str()); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%s", stringValue(info_mycuts.avg_bb_info.time, "%-*.*f", NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC).c_str()); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "\n"); output += tmpstring;
  } // vpcs
  if (branch_with_gmics) {
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%-*s", NAME_WIDTH, "Gur+V+G"); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%s", stringValue(info_allcuts.avg_bb_info.obj, "%-*.*f", NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC).c_str()); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%s", stringValue(info_allcuts.avg_bb_info.bound, "%-*.*f", NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC).c_str()); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%-*ld", NUM_DIGITS_BEFORE_DEC, info_allcuts.avg_bb_info.iters); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%-*ld", NUM_DIGITS_BEFORE_DEC, info_allcuts.avg_bb_info.nodes); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%-*ld", NUM_DIGITS_BEFORE_DEC, info_allcuts.avg_bb_info.root_passes); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%s", stringValue(info_allcuts.avg_bb_info.first_cut_pass, "%-*.*f", NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC).c_str()); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%s", stringValue(info_allcuts.avg_bb_info.last_cut_pass, "%-*.*f", NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC).c_str()); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%s", stringValue(info_allcuts.avg_bb_info.root_time, "%-*.*f", NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC).c_str()); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%s", stringValue(info_allcuts.avg_bb_info.last_sol_time, "%-*.*f", NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC).c_str()); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%s", stringValue(info_allcuts.avg_bb_info.time, "%-*.*f", NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC).c_str()); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "\n"); output += tmpstring;
  } // gmics
} /* analyzeBB */


//  char tmpstring[300];
//  std::string output = "";
//  sprintf(tmpstring, "\n## Results from adding cuts ##\n");
//  output += tmpstring;
//  sprintf(tmpstring, "%-*.*s%s\n", NAME_WIDTH, NAME_WIDTH, "LP: ",
//      stringValue(init_obj, "% -*.*f", NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC).c_str());
//  output += tmpstring;
//  if (!isInfinity(std::abs(gmic_obj))) {
//    sprintf(tmpstring, "%-*.*s%s (%d cuts)\n", NAME_WIDTH, NAME_WIDTH, "GMICs: ",
//        stringValue(gmic_obj, "% -*.*f", NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC).c_str(),
//        num_gmic_total);
//    output += tmpstring;
//  }
//  if (!isInfinity(std::abs(vpc_obj))) {
//    sprintf(tmpstring, "%-*.*s%s (%d cuts)\n", NAME_WIDTH, NAME_WIDTH, "VPCs: ",
//        stringValue(vpc_obj, "% -*.*f", NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC).c_str(),
//        num_vpc_total);
//    output += tmpstring;
//  }
//  if (!isInfinity(std::abs(best_disj_obj))) {
//    sprintf(tmpstring, "%-*.*s%s\n", NAME_WIDTH, NAME_WIDTH, "Disjunctive lb: ",
//        stringValue(best_disj_obj, "% -*.*f", NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC).c_str());
//    output += tmpstring;
//  }
//  if (!isInfinity(std::abs(ip_obj))) {
//    sprintf(tmpstring, "%-*.*s%s\n", NAME_WIDTH, NAME_WIDTH, "IP: ",
//        stringValue(ip_obj, "% -*.*f", NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC).c_str());
//    output += tmpstring;
//  }
//  printf("%s", output.c_str());
