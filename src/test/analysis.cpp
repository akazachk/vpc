// Name:     analysis.cpp
// Author:   A. M. Kazachkov
// Date:     2019-Mar-04
//-----------------------------------------------------------------------------
#include "analysis.hpp"

// COIN-OR
#include <OsiCuts.hpp>

// Project files
#include "BBHelper.hpp"
#include "CglVPC.hpp"
#include "VPCParameters.hpp"
#include "utility.hpp" // isInfinity, stringValue

const int countBoundInfoEntries = 11;
const int countGapInfoEntries = 4;
const int countBBInfoEntries = static_cast<int>(BB_INFO_CONTENTS.size()) * 4 * 2;
const int countOrigProbEntries = 14;
const int countPostCutProbEntries = 6;
const int countDisjInfoEntries = 17;
const int countCutInfoEntries = 10;
const int countObjInfoEntries = 1 + 4 * static_cast<int>(CglVPC::CutHeuristic::NUM_CUT_HEUR);
const int countFailInfoEntries = 1 + static_cast<int>(CglVPC::FailureType::NUM_FAILURES);
const int countParamInfoEntries = intParam::NUM_INT_PARAMS + doubleParam::NUM_DOUBLE_PARAMS;
int countTimeInfoEntries = 0; // set in printHeader

void printHeader(const VPCParameters& params,
    const std::vector<std::string>& time_name,
    const char SEP) {
  FILE* logfile = params.logfile;
  if (logfile == NULL)
    return;

  countTimeInfoEntries = time_name.size();

  // First line of the header details the categories of information displayed
  std::string tmpstring = "";
  fprintf(logfile, "%c", SEP);
  fprintf(logfile, "%s", "BOUND INFO");
  tmpstring.assign(countBoundInfoEntries, SEP);
  fprintf(logfile, "%s", tmpstring.c_str());
  fprintf(logfile, "%s", "GAP INFO");
  tmpstring.assign(countGapInfoEntries, SEP);
  fprintf(logfile, "%s", tmpstring.c_str());
  fprintf(logfile, "%s", "BB INFO");
  tmpstring.assign(countBBInfoEntries, SEP);
  fprintf(logfile, "%s", tmpstring.c_str());
  fprintf(logfile, "%s", "ORIG PROB");
  tmpstring.assign(countOrigProbEntries, SEP);
  fprintf(logfile, "%s", tmpstring.c_str());
  fprintf(logfile, "%s", "POST-CUT PROB");
  tmpstring.assign(countPostCutProbEntries, SEP);
  fprintf(logfile, "%s", tmpstring.c_str());
  fprintf(logfile, "%s", "DISJ INFO");
  tmpstring.assign(countDisjInfoEntries, SEP);
  fprintf(logfile, "%s", tmpstring.c_str());
  fprintf(logfile, "%s", "CUT INFO");
  tmpstring.assign(countCutInfoEntries, SEP);
  fprintf(logfile, "%s", tmpstring.c_str());
  fprintf(logfile, "%s", "OBJ INFO");
  tmpstring.assign(countObjInfoEntries, SEP);
  fprintf(logfile, "%s", tmpstring.c_str());
  fprintf(logfile, "%s", "FAIL INFO");
  tmpstring.assign(countFailInfoEntries, SEP);
  fprintf(logfile, "%s", tmpstring.c_str());
  fprintf(logfile, "%s", "PARAM INFO");
  tmpstring.assign(countParamInfoEntries, SEP);
  fprintf(logfile, "%s", tmpstring.c_str());
  fprintf(logfile, "%s", "TIME INFO");
  tmpstring.assign(countTimeInfoEntries, SEP);
  fprintf(logfile, "%s", tmpstring.c_str());
  fprintf(logfile, "%s", "END");
  fprintf(logfile, "\n");

  fprintf(logfile, "%s%c", "INSTANCE", SEP);
  { // BOUND INFO
    int count = 0;
    fprintf(logfile, "%s%c", "LP OBJ", SEP); count++;  // 1
    fprintf(logfile, "%s%c", "BEST DISJ OBJ", SEP); count++; // 2
    fprintf(logfile, "%s%c", "WORST DISJ OBJ", SEP); count++; // 3
    fprintf(logfile, "%s%c", "IP OBJ", SEP); count++; // 4
    fprintf(logfile, "%s%c", "NUM GMIC", SEP); count++; // 5
    fprintf(logfile, "%s%c", "GMIC OBJ", SEP); count++; // 6
    fprintf(logfile, "%s%c", "NUM L&PC", SEP); count++; // 7
    fprintf(logfile, "%s%c", "L&PC OBJ", SEP); count++; // 8
    fprintf(logfile, "%s%c", "NUM VPC", SEP); count++; // 9
    fprintf(logfile, "%s%c", "VPC OBJ", SEP); count++; // 10
    fprintf(logfile, "%s%c", "GMIC+VPC BOUND", SEP); count++; // 11
    assert(count == countBoundInfoEntries);
  } // BOUND INFO
  { // GAP INFO
    int count = 0;
    fprintf(logfile, "%s%c", "GMIC % GAP CLOSED", SEP); count++; // 1
    fprintf(logfile, "%s%c", "L&PC % GAP CLOSED", SEP); count++; // 2
    fprintf(logfile, "%s%c", "VPC % GAP CLOSED", SEP); count++; // 3
    fprintf(logfile, "%s%c", "GMIC+VPC % GAP CLOSED", SEP); count++; // 4
    assert(count == countGapInfoEntries);
  } // GAP INFO
  { // BB INFO
    int count = 0;
    for (std::string name : BB_INFO_CONTENTS) {
      fprintf(logfile, "%s%c", ("GUR1 " + name).c_str(), SEP); count++;
      fprintf(logfile, "%s%c", ("GUR1+V " + name).c_str(), SEP); count++;
      fprintf(logfile, "%s%c", ("BEST GUR " + name).c_str(), SEP); count++;
      fprintf(logfile, "%s%c", ("BEST GUR+V " + name).c_str(), SEP); count++;
      fprintf(logfile, "%s%c", ("AVG GUR " + name).c_str(), SEP); count++;
      fprintf(logfile, "%s%c", ("AVG GUR+V " + name).c_str(), SEP); count++;
    }
    for (std::string name : BB_INFO_CONTENTS) {
      fprintf(logfile, "%s%c", ("ALL GUR " + name).c_str(), SEP); count++;
    }
    for (std::string name : BB_INFO_CONTENTS) {
      fprintf(logfile, "%s%c", ("ALL GUR+V " + name).c_str(), SEP); count++;
    }
    assert(count == countBBInfoEntries);
  } // BB INFO
  { // ORIG PROB
    int count = 0;
    fprintf(logfile, "%s%c", "NUM ROWS", SEP); count++;
    fprintf(logfile, "%s%c", "NUM COLS", SEP); count++;
    fprintf(logfile, "%s%c", "NUM FRAC", SEP); count++;
    fprintf(logfile, "%s%c", "MIN FRACTIONALITY", SEP); count++;
    fprintf(logfile, "%s%c", "MAX FRACTIONALITY", SEP); count++;
    fprintf(logfile, "%s%c", "NUM FIXED ROWS", SEP); count++;
    fprintf(logfile, "%s%c", "NUM EQ ROWS", SEP); count++;
    fprintf(logfile, "%s%c", "NUM INEQ ROWS", SEP); count++;
    fprintf(logfile, "%s%c", "NUM BOUND ROWS", SEP); count++;
    fprintf(logfile, "%s%c", "NUM ASSIGN ROWS", SEP); count++;
    fprintf(logfile, "%s%c", "GEN INT", SEP); count++;
    fprintf(logfile, "%s%c", "BINARY", SEP); count++;
    fprintf(logfile, "%s%c", "CONTINUOUS", SEP); count++;
    fprintf(logfile, "%s%c", "A-DENSITY", SEP); count++;
    assert(count == countOrigProbEntries);
  } // ORIG PROB
  { // POST-CUT PROB
    int count = 0;
    fprintf(logfile, "%s%c", "NEW NUM FRAC", SEP); count++;
    fprintf(logfile, "%s%c", "NEW MIN FRACTIONALITY", SEP); count++;
    fprintf(logfile, "%s%c", "NEW MAX FRACTIONALITY", SEP); count++;
    fprintf(logfile, "%s%c", "NEW A-DENSITY", SEP); count++;
    fprintf(logfile, "%s%c", "ACTIVE GMIC", SEP); count++;
    fprintf(logfile, "%s%c", "ACTIVE VPC", SEP); count++;
    assert(count == countPostCutProbEntries);
  } // POST-CUT PROB
  { // DISJ INFO
    int count = 0;
    fprintf(logfile, "%s%c", "NUM DISJ TERMS", SEP); count++;
    fprintf(logfile, "%s%c", "MIN DENSITY PRLP", SEP); count++;
    fprintf(logfile, "%s%c", "MAX DENSITY PRLP", SEP); count++;
    fprintf(logfile, "%s%c", "MIN NUM ROWS PRLP", SEP); count++;
    fprintf(logfile, "%s%c", "MAX NUM ROWS PRLP", SEP); count++;
    fprintf(logfile, "%s%c", "MIN NUM COLS PRLP", SEP); count++;
    fprintf(logfile, "%s%c", "MAX NUM COLS PRLP", SEP); count++;
    fprintf(logfile, "%s%c", "MIN NUM POINTS PRLP", SEP); count++;
    fprintf(logfile, "%s%c", "MAX NUM POINTS PRLP", SEP); count++;
    fprintf(logfile, "%s%c", "TOTAL NUM POINTS PRLP", SEP); count++;
    fprintf(logfile, "%s%c", "MIN NUM RAYS PRLP", SEP); count++;
    fprintf(logfile, "%s%c", "MAX NUM RAYS PRLP", SEP); count++;
    fprintf(logfile, "%s%c", "TOTAL NUM RAYS PRLP", SEP); count++;
    fprintf(logfile, "%s%c", "NUM PARTIAL BB NODES", SEP); count++;
    fprintf(logfile, "%s%c", "NUM PRUNED NODES", SEP); count++;
    fprintf(logfile, "%s%c", "MIN DEPTH", SEP); count++;
    fprintf(logfile, "%s%c", "MAX DEPTH", SEP); count++;
    assert(count == countDisjInfoEntries);
  } // DISJ INFO
  { // CUT INFO
    int count = 0;
    fprintf(logfile, "%s%c", "NUM CUTS", SEP); count++;
    fprintf(logfile, "%s%c", "NUM ONE SIDED", SEP); count++;
    fprintf(logfile, "%s%c", "MIN SUPPORT GOMORY", SEP); count++;
    fprintf(logfile, "%s%c", "MAX SUPPORT GOMORY", SEP); count++;
    fprintf(logfile, "%s%c", "AVG SUPPORT GOMORY", SEP); count++;
    fprintf(logfile, "%s%c", "STD DEV SUPPORT GOMORY", SEP); count++;
    fprintf(logfile, "%s%c", "MIN SUPPORT VPC", SEP); count++;
    fprintf(logfile, "%s%c", "MAX SUPPORT VPC", SEP); count++;
    fprintf(logfile, "%s%c", "AVG SUPPORT VPC", SEP); count++;
    fprintf(logfile, "%s%c", "STD DEV SUPPORT VPC", SEP); count++;
    assert(count == countCutInfoEntries);
  } // CUT INFO
  { // OBJ INFO
    // For each objective: num obj, num fails, num active
    int count = 0;
    fprintf(logfile, "%s%c", "NUM OBJ", SEP); count++;
    for (int fail_ind = 0; fail_ind < static_cast<int>(CglVPC::CutHeuristic::NUM_CUT_HEUR); fail_ind++) {
      fprintf(logfile, "NUM OBJ %s%c", CglVPC::CutHeuristicName[fail_ind].c_str(), SEP); count++;
      fprintf(logfile, "NUM FAILS %s%c", CglVPC::CutHeuristicName[fail_ind].c_str(), SEP); count++;
      fprintf(logfile, "NUM CUTS %s%c", CglVPC::CutHeuristicName[fail_ind].c_str(), SEP); count++;
      fprintf(logfile, "NUM ACTIVE %s%c", CglVPC::CutHeuristicName[fail_ind].c_str(), SEP); count++;
    }
    assert(count == countObjInfoEntries);
  } // OBJ INFO
  { // FAIL INFO
    int count = 0;
    fprintf(logfile, "%s%c", "NUM FAILS", SEP); count++;
    for (int fail_ind = 0; fail_ind < static_cast<int>(CglVPC::FailureType::NUM_FAILURES); fail_ind++) {
      fprintf(logfile, "%s%c", CglVPC::FailureTypeName[fail_ind].c_str(), SEP); count++;
    }
    assert(count == countFailInfoEntries);
  } // FAIL INFO
  { // PARAM INFO
    printParams(params, logfile, 1); // only int/double param names
  } // PARAM INFO
  { // TIME INFO
    int count = 0;
    for (int t = 0; t < (int) time_name.size(); t++) {
      fprintf(logfile, "%s%c", time_name[t].c_str(), SEP); count++;
    }
    assert(count == countTimeInfoEntries);
  } // TIME INFO

  fprintf(logfile, "\n");
  fflush(logfile);
} /* printHeader */

void printBoundAndGapInfo(const SummaryBoundInfo& boundInfo, FILE* logfile, const char SEP) {
  if (!logfile)
    return;

  { // BOUND INFO
    int count = 0;
    fprintf(logfile, "%s%c", stringValue(boundInfo.lp_obj, "%2.20f").c_str(), SEP); count++;
    fprintf(logfile, "%s%c", stringValue(boundInfo.best_disj_obj, "%2.20f").c_str(), SEP); count++;
    fprintf(logfile, "%s%c", stringValue(boundInfo.worst_disj_obj, "%2.20f").c_str(), SEP); count++;
    if (!isInfinity(std::abs(boundInfo.ip_obj))) {
      fprintf(logfile, "%s%c", stringValue(boundInfo.ip_obj, "%2.20f").c_str(), SEP); count++;
    } else {
      fprintf(logfile, "%c", SEP); count++;
    }
    fprintf(logfile, "%s%c", stringValue(boundInfo.num_gmic).c_str(), SEP); count++;
    if (boundInfo.num_gmic > 0) {
      fprintf(logfile, "%s%c", stringValue(boundInfo.gmic_obj, "%2.20f").c_str(), SEP); count++;
    } else {
      fprintf(logfile, "%c", SEP); count++;
    }
    fprintf(logfile, "%s%c", stringValue(boundInfo.num_lpc).c_str(), SEP); count++;
    if (boundInfo.num_lpc > 0) {
      fprintf(logfile, "%s%c", stringValue(boundInfo.lpc_obj, "%2.20f").c_str(), SEP); count++;
    } else {
      fprintf(logfile, "%c", SEP); count++;
    }
    fprintf(logfile, "%s%c", stringValue(boundInfo.num_vpc).c_str(), SEP); count++;
    if (boundInfo.num_vpc > 0) {
      fprintf(logfile, "%s%c", stringValue(boundInfo.vpc_obj, "%2.20f").c_str(), SEP); count++;
    } else {
      fprintf(logfile, "%c", SEP); count++;
    }
    if (!isInfinity(std::abs(boundInfo.gmic_vpc_obj))) {
      fprintf(logfile, "%s%c", stringValue(boundInfo.gmic_vpc_obj, "%2.20f").c_str(), SEP); count++;
    } else {
      fprintf(logfile, "%c", SEP); count++;
    }
    assert(count == countBoundInfoEntries);
  } // BOUND INFO

  { // GAP INFO
    int count = 0;
    if (!isInfinity(std::abs(boundInfo.ip_obj))) {
      if (!isInfinity(std::abs(boundInfo.gmic_obj))) {
        double val = 100. * (boundInfo.gmic_obj - boundInfo.lp_obj)
            / (boundInfo.ip_obj - boundInfo.lp_obj);
        fprintf(logfile, "%s%c", stringValue(val, "%2.6f").c_str(), SEP); count++;
      } else {
        fprintf(logfile, "%c", SEP); count++; // gmic
      }
      if (!isInfinity(std::abs(boundInfo.lpc_obj))) {
        double val = 100. * (boundInfo.lpc_obj - boundInfo.lp_obj)
            / (boundInfo.ip_obj - boundInfo.lp_obj);
        fprintf(logfile, "%s%c", stringValue(val, "%2.6f").c_str(), SEP); count++;
      } else {
        fprintf(logfile, "%c", SEP); count++; // lpc
      }
      if (!isInfinity(std::abs(boundInfo.vpc_obj))) {
        double val = 100. * (boundInfo.vpc_obj - boundInfo.lp_obj)
            / (boundInfo.ip_obj - boundInfo.lp_obj);
        fprintf(logfile, "%s%c", stringValue(val, "%2.6f").c_str(), SEP); count++;
      } else {
        fprintf(logfile, "%c", SEP); count++; // vpc
      }
      if (!isInfinity(std::abs(boundInfo.gmic_vpc_obj))) {
        double val = 100. * (boundInfo.gmic_vpc_obj - boundInfo.lp_obj)
            / (boundInfo.ip_obj - boundInfo.lp_obj);
        fprintf(logfile, "%s%c", stringValue(val, "%2.6f").c_str(), SEP); count++;
      } else {
        fprintf(logfile, "%c", SEP); count++; // gmic_vpc
      }
    } else {
      fprintf(logfile, "%c", SEP); count++; // gmic
      fprintf(logfile, "%c", SEP); count++; // lpc
      fprintf(logfile, "%c", SEP); count++; // vpc
      fprintf(logfile, "%c", SEP); count++; // gmic+vpc
    }
    assert(count == countGapInfoEntries);
  }
} /* printBoundAndGapInfo */

void printBBInfo(const std::vector<SummaryBBInfo>& info_vec, FILE* logfile,
    const bool print_blanks, const char SEP) {
  if (!logfile)
    return;

//  const std::vector<bool> did_branch(info_vec.size());
//  for (unsigned i = 0; i < info_vec.size(); i++) {
//    did_branch[i] = info_vec[i].vec_bb_info.size() > 0;
//  }

  int count = 0;
  if (!print_blanks) {
    for (auto& info : info_vec) {
      fprintf(logfile, "%s%c", stringValue(info.first_bb_info.obj, "%2.20f").c_str(), SEP);
      count++;
    }
    for (auto& info : info_vec) {
      fprintf(logfile, "%s%c", stringValue(info.best_bb_info.obj, "%2.20f").c_str(), SEP);
      count++;
    }
    for (auto& info : info_vec) {
      fprintf(logfile, "%s%c", stringValue(info.avg_bb_info.obj, "%2.20f").c_str(), SEP);
      count++;
    }
    for (auto info : info_vec) {
      fprintf(logfile, "%s%c", stringValue(info.first_bb_info.bound, "%2.20f").c_str(), SEP);
      count++;
    }
    for (auto info : info_vec) {
      fprintf(logfile, "%s%c", stringValue(info.best_bb_info.bound, "%2.20f").c_str(), SEP);
      count++;
    }
    for (auto info : info_vec) {
      fprintf(logfile, "%s%c", stringValue(info.avg_bb_info.bound, "%2.20f").c_str(), SEP);
      count++;
    }
    for (auto info : info_vec) {
      fprintf(logfile, "%ld%c", info.first_bb_info.iters, SEP);
      count++;
    }
    for (auto info : info_vec) {
      fprintf(logfile, "%ld%c", info.best_bb_info.iters, SEP);
      count++;
    }
    for (auto info : info_vec) {
      fprintf(logfile, "%ld%c", info.avg_bb_info.iters, SEP);
      count++;
    }
    for (auto info : info_vec) {
      fprintf(logfile, "%ld%c", info.first_bb_info.nodes, SEP);
      count++;
    }
    for (auto info : info_vec) {
      fprintf(logfile, "%ld%c", info.best_bb_info.nodes, SEP);
      count++;
    }
    for (auto info : info_vec) {
      fprintf(logfile, "%ld%c", info.avg_bb_info.nodes, SEP);
      count++;
    }
    for (auto info : info_vec) {
      fprintf(logfile, "%ld%c", info.first_bb_info.root_passes, SEP);
      count++;
    }
    for (auto info : info_vec) {
      fprintf(logfile, "%ld%c", info.best_bb_info.root_passes, SEP);
      count++;
    }
    for (auto info : info_vec) {
      fprintf(logfile, "%ld%c", info.avg_bb_info.root_passes, SEP);
      count++;
    }
    for (auto info : info_vec) {
      fprintf(logfile, "%2.20f%c", info.first_bb_info.first_cut_pass, SEP);
      count++;
    }
    for (auto info : info_vec) {
      fprintf(logfile, "%2.20f%c", info.best_bb_info.first_cut_pass, SEP);
      count++;
    }
    for (auto info : info_vec) {
      fprintf(logfile, "%2.20f%c", info.avg_bb_info.first_cut_pass, SEP);
      count++;
    }
    for (auto info : info_vec) {
      fprintf(logfile, "%2.20f%c", info.first_bb_info.last_cut_pass, SEP);
      count++;
    }
    for (auto info : info_vec) {
      fprintf(logfile, "%2.20f%c", info.best_bb_info.last_cut_pass, SEP);
      count++;
    }
    for (auto info : info_vec) {
      fprintf(logfile, "%2.20f%c", info.avg_bb_info.last_cut_pass, SEP);
      count++;
    }
    for (auto info : info_vec) {
      fprintf(logfile, "%2.3f%c", info.first_bb_info.root_time, SEP);
      count++;
    }
    for (auto info : info_vec) {
      fprintf(logfile, "%2.3f%c", info.best_bb_info.root_time, SEP);
      count++;
    }
    for (auto info : info_vec) {
      fprintf(logfile, "%2.3f%c", info.avg_bb_info.root_time, SEP);
      count++;
    }
    for (auto info : info_vec) {
      fprintf(logfile, "%2.3f%c", info.first_bb_info.last_sol_time, SEP);
      count++;
    }
    for (auto info : info_vec) {
      fprintf(logfile, "%2.3f%c", info.best_bb_info.last_sol_time, SEP);
      count++;
    }
    for (auto info : info_vec) {
      fprintf(logfile, "%2.3f%c", info.avg_bb_info.last_sol_time, SEP);
      count++;
    }
    for (auto info : info_vec) {
      fprintf(logfile, "%2.3f%c", info.first_bb_info.time, SEP);
      count++;
    }
    for (auto info : info_vec) {
      fprintf(logfile, "%2.3f%c", info.best_bb_info.time, SEP);
      count++;
    }
    for (auto info : info_vec) {
      fprintf(logfile, "%2.3f%c", info.avg_bb_info.time, SEP);
      count++;
    }

    // Finally, all
//    for (unsigned i = 0; i < info_vec.size(); i++) {
    for (auto info : info_vec) {
      std::vector<std::string> vec_str;
      createStringFromBBInfoVec(info.vec_bb_info, vec_str);
      for (unsigned i = 0; i < vec_str.size(); i++) {
        fprintf(logfile, "%s%c", vec_str[i].c_str(), SEP);
        count++;
      }
    }
  } else {
    for (int i = 0; i < (int) BB_INFO_CONTENTS.size() * info_vec.size(); i++) {
      fprintf(logfile, "%c", SEP); count++;
    }
  }
  assert(count == countBBInfoEntries);
} /* printBBInfo */
//
///**
// * amountToPrint:
// *  0 = all (newline-separated),
// *  1 = only header (comma-separated),
// *  2 = only values (comma-separated)
// */
//void printBBInfo(const SummaryBBInfo& info_nocuts,
//    const SummaryBBInfo& info_mycuts,
////    const SummaryBBInfo& info_allcuts,
//    FILE *logfile, const int amountToPrint, const char SEP) {
//  if (logfile == NULL) {
//    return;
//  }
//
//  //////////////////// BB INFO
//  switch (amountToPrint) {
//    case 0: {
//      std::cerr << "printBBInfo: Case 0 of amountToPrint not implmented\n";
//      exit(1);
////      for (std::string name : BB_INFO_CONTENTS) {
////        fprintf(myfile, "%s%c", ("FIRST BB " + name + " VPC").c_str(), SEP);
////        fprintf(myfile, "%s%c", ("FIRST BB " + name + " VPC+GOMORY").c_str(), SEP);
////      }
////      for (std::string name : BB_INFO_CONTENTS) {
////        fprintf(myfile, "%s%c", ("BEST BB " + name + " VPC").c_str(), SEP);
////        fprintf(myfile, "%s%c", ("BEST BB " + name + " VPC+GOMORY").c_str(), SEP);
////      }
////      for (std::string name : BB_INFO_CONTENTS) {
////        fprintf(myfile, "%s%c", ("AVG BB " + name + " VPC").c_str(), SEP);
////        fprintf(myfile, "%s%c", ("AVG BB " + name + " VPC+GOMORY").c_str(), SEP);
////      }
////      for (std::string name : BB_INFO_CONTENTS) {
////        fprintf(myfile, "%s%c", ("ALL BB " + name + " VPC").c_str(), SEP);
////        fprintf(myfile, "%s%c", ("ALL BB " + name + " VPC+GOMORY").c_str(), SEP);
////      }
//      break;
//    }
//    case 1: {
//      int count = 0;
//      for (std::string name : BB_INFO_CONTENTS) {
//        fprintf(logfile, "%s%c", ("GUR1 " + name).c_str(), SEP); count++;
//        fprintf(logfile, "%s%c", ("GUR1+V " + name).c_str(), SEP); count++;
////        fprintf(logfile, "%s%c", ("GUR1+V+G " + name).c_str(), SEP); count++;
//      }
//      for (std::string name : BB_INFO_CONTENTS) {
//        fprintf(logfile, "%s%c", ("BEST GUR " + name).c_str(), SEP); count++;
//        fprintf(logfile, "%s%c", ("BEST GUR+V " + name).c_str(), SEP); count++;
////        fprintf(logfile, "%s%c", ("BEST GUR+V+G " + name).c_str(), SEP); count++;
//      }
//      for (std::string name : BB_INFO_CONTENTS) {
//        fprintf(logfile, "%s%c", ("AVG GUR " + name).c_str(), SEP); count++;
//        fprintf(logfile, "%s%c", ("AVG GUR+V " + name).c_str(), SEP); count++;
////        fprintf(logfile, "%s%c", ("AVG GUR+V+G " + name).c_str(), SEP); count++;
//      }
//      for (std::string name : BB_INFO_CONTENTS) {
//        fprintf(logfile, "%s%c", ("ALL GUR " + name).c_str(), SEP); count++;
//        fprintf(logfile, "%s%c", ("ALL GUR+V " + name).c_str(), SEP); count++;
////        fprintf(logfile, "%s%c", ("ALL GUR+V+G " + name).c_str(), SEP); count++;
//      }
//      assert(count == countBBInfoEntries);
//      break;
//    }
//    case 2: {
//      if (info_mycuts.vec_bb_info.size() == 0) {
//        for (int i = 0; i < countBBInfoEntries; i++) {
//          fprintf(logfile, "%c", SEP);
//        }
//        fflush(logfile);
//        break;
//      }
//      // First
//      printBBInfo(info_nocuts.vec_bb_info[0], info_mycuts.vec_bb_info[0],
////          info_allcuts.vec_bb_info[0],
//          logfile, false, SEP);
//
//      // Min
//      printBBInfo(info_nocuts.best_bb_info, info_mycuts.best_bb_info,
////          info_allcuts.best_bb_info,
//          logfile, false, SEP);
//
//      // Average
//      printBBInfo(info_nocuts.avg_bb_info, info_mycuts.avg_bb_info,
////          info_allcuts.avg_bb_info,
//          logfile, false, SEP);
//
//      // All
//      std::vector<std::string> vec_str_nocuts, vec_str_mycuts;
////      std::vector<std::string> vec_str_allcuts;
//      createStringFromBBInfoVec(info_nocuts.vec_bb_info, vec_str_nocuts);
//      createStringFromBBInfoVec(info_mycuts.vec_bb_info, vec_str_mycuts);
////      createStringFromBBInfoVec(info_allcuts.vec_bb_info, vec_str_allcuts);
//      for (int i = 0; i < (int) vec_str_mycuts.size(); i++) {
//        fprintf(logfile, "%s%c", vec_str_nocuts[i].c_str(), SEP);
//        fprintf(logfile, "%s%c", vec_str_mycuts[i].c_str(), SEP);
////        fprintf(logfile, "%s%c", vec_str_allcuts[i].c_str(), SEP);
//      }
//      break;
//    }
//    default: {
//      // nothing
//    }
//  } // switch amountToPrint
//  fflush(logfile);
//} /* printBBInfo */

/**
 * The cut properties we want to look at are:
 * 1. Gap closed
 * 2. Activity (after adding cuts)
 * 3. Density
 */
void analyzeStrength(const VPCParameters& params, const SummaryBoundInfo& boundInfo,
    const OsiCuts* const vpc, std::string& output) {
  // Print results from adding cuts
  int NAME_WIDTH = 25;
  int NUM_DIGITS_BEFORE_DEC = 7;
  int NUM_DIGITS_AFTER_DEC = 7;
  char tmpstring[300];

  snprintf(tmpstring, sizeof(tmpstring) / sizeof(char),
      "\n## Results from adding cuts ##\n");
  output += tmpstring;
  snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%-*.*s%s\n",
      NAME_WIDTH, NAME_WIDTH, "LP: ",
      stringValue(boundInfo.lp_obj, "% -*.*f", NUM_DIGITS_BEFORE_DEC,
          NUM_DIGITS_AFTER_DEC).c_str());
  output += tmpstring;
  if (!isInfinity(std::abs(boundInfo.gmic_obj))) {
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char),
        "%-*.*s%s (%d cuts)\n", NAME_WIDTH, NAME_WIDTH, "GMICs: ",
        stringValue(boundInfo.gmic_obj, "% -*.*f", NUM_DIGITS_BEFORE_DEC,
            NUM_DIGITS_AFTER_DEC).c_str(), boundInfo.num_gmic);
    output += tmpstring;
  }
  if (!isInfinity(std::abs(boundInfo.vpc_obj))) {
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char),
        "%-*.*s%s (%d cuts)\n", NAME_WIDTH, NAME_WIDTH, "VPCs: ",
        stringValue(boundInfo.vpc_obj, "% -*.*f", NUM_DIGITS_BEFORE_DEC,
            NUM_DIGITS_AFTER_DEC).c_str(), boundInfo.num_vpc);
    output += tmpstring;
  }
  if (!isInfinity(std::abs(boundInfo.best_disj_obj))) {
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%-*.*s%s\n",
        NAME_WIDTH, NAME_WIDTH, "Disjunctive lb: ",
        stringValue(boundInfo.best_disj_obj, "% -*.*f", NUM_DIGITS_BEFORE_DEC,
            NUM_DIGITS_AFTER_DEC).c_str());
    output += tmpstring;
  }
  if (!isInfinity(std::abs(boundInfo.worst_disj_obj))) {
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%-*.*s%s\n",
        NAME_WIDTH, NAME_WIDTH, "Disjunctive ub: ",
        stringValue(boundInfo.worst_disj_obj, "% -*.*f", NUM_DIGITS_BEFORE_DEC,
            NUM_DIGITS_AFTER_DEC).c_str());
    output += tmpstring;
  }
  if (!isInfinity(std::abs(boundInfo.ip_obj))) {
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%-*.*s%s\n",
        NAME_WIDTH, NAME_WIDTH, "IP: ",
        stringValue(boundInfo.ip_obj, "% -*.*f", NUM_DIGITS_BEFORE_DEC,
            NUM_DIGITS_AFTER_DEC).c_str());
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
