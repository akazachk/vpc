/**
 * @file analysis.cpp
 * @author A. M. Kazachkov
 * @date 2019-03-04
 */
#include "analysis.hpp"

// COIN-OR
#include <OsiSolverInterface.hpp>
#include <OsiCuts.hpp>
#include <CglGMI.hpp>

// Project files
#include "BBHelper.hpp"
#include "CglVPC.hpp"
#include "CutHelper.hpp"
#include "Disjunction.hpp"
#include "PartialBBDisjunction.hpp"
#include "PRLP.hpp"
#include "SolverHelper.hpp"
#include "VPCParameters.hpp"
using namespace VPCParametersNamespace;
#include "utility.hpp" // isInfinity, stringValue

const int countBoundInfoEntries = 11;
const int countGapInfoEntries = 4;
const int countSummaryBBInfoEntries = 4 * 2;
const int countFullBBInfoEntries = static_cast<int>(BB_INFO_CONTENTS.size()) * 4 * 2;
const int countOrigProbEntries = 13;
const int countPostCutProbEntries = 10;
const int countDisjInfoEntries = 12;
const int countCutInfoEntries = 10;
const int countObjInfoEntries = 1 + 4 * static_cast<int>(CglVPC::ObjectiveType::NUM_OBJECTIVE_TYPES);
const int countFailInfoEntries = 1 + static_cast<int>(CglVPC::FailureType::NUM_FAILURE_TYPES);
const int countParamInfoEntries = intParam::NUM_INT_PARAMS + doubleParam::NUM_DOUBLE_PARAMS + stringParam::NUM_STRING_PARAMS;
int countTimeInfoEntries = 0; // set in printHeader
const int countVersionInfoEntries = 5;
const int countExtraInfoEntries = 4;

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
  fprintf(logfile, "%s", "PARAM INFO");
  tmpstring.assign(countParamInfoEntries, SEP);
  fprintf(logfile, "%s", tmpstring.c_str());
  fprintf(logfile, "%s", "BOUND INFO");
  tmpstring.assign(countBoundInfoEntries, SEP);
  fprintf(logfile, "%s", tmpstring.c_str());
  fprintf(logfile, "%s", "GAP INFO");
  tmpstring.assign(countGapInfoEntries, SEP);
  fprintf(logfile, "%s", tmpstring.c_str());
  fprintf(logfile, "%s", "BB INFO");
  tmpstring.assign(countSummaryBBInfoEntries, SEP);
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
  fprintf(logfile, "%s", "FULL BB INFO");
  tmpstring.assign(countFullBBInfoEntries, SEP);
  fprintf(logfile, "%s", tmpstring.c_str());
  fprintf(logfile, "%s", "TIME INFO");
  tmpstring.assign(countTimeInfoEntries, SEP);
  fprintf(logfile, "%s", tmpstring.c_str());
  fprintf(logfile, "%s", "VERSION INFO");
  tmpstring.assign(countVersionInfoEntries, SEP);
  fprintf(logfile, "%s", tmpstring.c_str());
  fprintf(logfile, "%s", "WRAP UP INFO");
  tmpstring.assign(countExtraInfoEntries, SEP);
  fprintf(logfile, "%s", tmpstring.c_str());
  fprintf(logfile, "%s", "END");
  fprintf(logfile, "\n");

  fprintf(logfile, "%s%c", "INSTANCE", SEP);
  { // PARAM INFO
    //printParams(params, logfile, 1); // names for int/double params
    printParams(params, logfile, 7); // names for int/double/string params
  } // PARAM INFO
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
    fprintf(logfile, "%s%c", "VPC+GMIC OBJ", SEP); count++; // 11
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
    std::vector<std::string> nameVec = {"NODES", "TIME"};
    for (auto name : nameVec) {
      fprintf(logfile, "%s%c", ("FIRST REF " + name).c_str(), SEP); count++;
      fprintf(logfile, "%s%c", ("FIRST REF+V " + name).c_str(), SEP); count++;
      fprintf(logfile, "%s%c", ("BEST REF " + name).c_str(), SEP); count++;
      fprintf(logfile, "%s%c", ("BEST REF+V " + name).c_str(), SEP); count++;
    }
    assert(count == countSummaryBBInfoEntries);
  } // BB INFO
  { // ORIG PROB
    int count = 0;
    fprintf(logfile, "%s%c", "ROWS", SEP); count++;
    fprintf(logfile, "%s%c", "COLS", SEP); count++;
    fprintf(logfile, "%s%c", "NUM FRAC", SEP); count++;
    fprintf(logfile, "%s%c", "MIN FRACTIONALITY", SEP); count++;
    fprintf(logfile, "%s%c", "MAX FRACTIONALITY", SEP); count++;
    fprintf(logfile, "%s%c", "EQ ROWS", SEP); count++;
    fprintf(logfile, "%s%c", "BOUND ROWS", SEP); count++;
    fprintf(logfile, "%s%c", "ASSIGN ROWS", SEP); count++;
    fprintf(logfile, "%s%c", "FIXED COLS", SEP); count++;
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
    fprintf(logfile, "%s%c", "ACTIVE GMIC (gmics)", SEP); count++;
    fprintf(logfile, "%s%c", "ACTIVE VPC (gmics)", SEP); count++;
    fprintf(logfile, "%s%c", "ACTIVE GMIC (vpcs)", SEP); count++;
    fprintf(logfile, "%s%c", "ACTIVE VPC (vpcs)", SEP); count++;
    fprintf(logfile, "%s%c", "ACTIVE GMIC (all cuts)", SEP); count++;
    fprintf(logfile, "%s%c", "ACTIVE VPC (all cuts)", SEP); count++;
    assert(count == countPostCutProbEntries);
  } // POST-CUT PROB
  { // DISJ INFO
    int count = 0;
    fprintf(logfile, "%s%c", "NUM DISJ TERMS", SEP); count++;
    fprintf(logfile, "%s%c", "NUM INTEGER SOL", SEP); count++;
    fprintf(logfile, "%s%c", "NUM DISJ", SEP); count++;
//    fprintf(logfile, "%s%c", "MIN DENSITY PRLP", SEP); count++;
//    fprintf(logfile, "%s%c", "MAX DENSITY PRLP", SEP); count++;
    fprintf(logfile, "%s%c", "AVG DENSITY PRLP", SEP); count++;
//    fprintf(logfile, "%s%c", "MIN ROWS PRLP", SEP); count++;
//    fprintf(logfile, "%s%c", "MAX ROWS PRLP", SEP); count++;
    fprintf(logfile, "%s%c", "AVG ROWS PRLP", SEP); count++;
//    fprintf(logfile, "%s%c", "MIN COLS PRLP", SEP); count++;
//    fprintf(logfile, "%s%c", "MAX COLS PRLP", SEP); count++;
    fprintf(logfile, "%s%c", "AVG COLS PRLP", SEP); count++;
//    fprintf(logfile, "%s%c", "MIN POINTS PRLP", SEP); count++;
//    fprintf(logfile, "%s%c", "MAX POINTS PRLP", SEP); count++;
    fprintf(logfile, "%s%c", "AVG POINTS PRLP", SEP); count++;
//    fprintf(logfile, "%s%c", "MIN RAYS PRLP", SEP); count++;
//    fprintf(logfile, "%s%c", "MAX RAYS PRLP", SEP); count++;
    fprintf(logfile, "%s%c", "AVG RAYS PRLP", SEP); count++;
    fprintf(logfile, "%s%c", "AVG PARTIAL BB EXPLORED NODES", SEP); count++;
    fprintf(logfile, "%s%c", "AVG PARTIAL BB PRUNED NODES", SEP); count++;
    fprintf(logfile, "%s%c", "AVG PARTIAL BB MIN DEPTH", SEP); count++;
    fprintf(logfile, "%s%c", "AVG PARTIAL BB MAX DEPTH", SEP); count++;
    assert(count == countDisjInfoEntries);
  } // DISJ INFO
  { // CUT INFO
    int count = 0;
    fprintf(logfile, "%s%c", "NUM ROUNDS", SEP); count++;
    fprintf(logfile, "%s%c", "NUM CUTS", SEP); count++; // repeat, but it's ok
    fprintf(logfile, "%s%c", "NUM ONE SIDED CUTS", SEP); count++;
    fprintf(logfile, "%s%c", "NUM OPTIMALITY CUTS", SEP); count++;
    fprintf(logfile, "%s%c", "MIN SUPPORT VPC", SEP); count++;
    fprintf(logfile, "%s%c", "MAX SUPPORT VPC", SEP); count++;
    fprintf(logfile, "%s%c", "AVG SUPPORT VPC", SEP); count++;
    fprintf(logfile, "%s%c", "MIN SUPPORT GOMORY", SEP); count++;
    fprintf(logfile, "%s%c", "MAX SUPPORT GOMORY", SEP); count++;
    fprintf(logfile, "%s%c", "AVG SUPPORT GOMORY", SEP); count++;
    assert(count == countCutInfoEntries);
  } // CUT INFO
  { // OBJ INFO
    // For each objective: num obj, num fails, num active
    int count = 0;
    fprintf(logfile, "%s%c", "NUM OBJ", SEP); count++;
    for (int obj_ind = 0; obj_ind < static_cast<int>(CglVPC::ObjectiveType::NUM_OBJECTIVE_TYPES); obj_ind++) {
      fprintf(logfile, "NUM OBJ %s%c", CglVPC::ObjectiveTypeName[obj_ind].c_str(), SEP); count++;
      fprintf(logfile, "NUM CUTS %s%c", CglVPC::ObjectiveTypeName[obj_ind].c_str(), SEP); count++;
      fprintf(logfile, "NUM FAILS %s%c", CglVPC::ObjectiveTypeName[obj_ind].c_str(), SEP); count++;
      fprintf(logfile, "NUM ACTIVE %s%c", CglVPC::ObjectiveTypeName[obj_ind].c_str(), SEP); count++;
    }
    assert(count == countObjInfoEntries);
  } // OBJ INFO
  { // FAIL INFO
    int count = 0;
    fprintf(logfile, "%s%c", "NUM FAILS", SEP); count++;
    for (int fail_ind = 0; fail_ind < static_cast<int>(CglVPC::FailureType::NUM_FAILURE_TYPES); fail_ind++) {
      fprintf(logfile, "%s%c", CglVPC::FailureTypeName[fail_ind].c_str(), SEP); count++;
    }
    assert(count == countFailInfoEntries);
  } // FAIL INFO
  { // FULL BB INFO
    int count = 0;
    for (std::string name : BB_INFO_CONTENTS) {
      fprintf(logfile, "%s%c", ("FIRST REF " + name).c_str(), SEP); count++;
      fprintf(logfile, "%s%c", ("FIRST REF+V " + name).c_str(), SEP); count++;
      fprintf(logfile, "%s%c", ("BEST REF " + name).c_str(), SEP); count++;
      fprintf(logfile, "%s%c", ("BEST REF+V " + name).c_str(), SEP); count++;
      fprintf(logfile, "%s%c", ("AVG REF " + name).c_str(), SEP); count++;
      fprintf(logfile, "%s%c", ("AVG REF+V " + name).c_str(), SEP); count++;
    }
    for (std::string name : BB_INFO_CONTENTS) {
      fprintf(logfile, "%s%c", ("ALL REF " + name).c_str(), SEP); count++;
    }
    for (std::string name : BB_INFO_CONTENTS) {
      fprintf(logfile, "%s%c", ("ALL REF+V " + name).c_str(), SEP); count++;
    }
    assert(count == countFullBBInfoEntries);
  } // FULL BB INFO
  { // TIME INFO
    int count = 0;
    for (int t = 0; t < (int) time_name.size(); t++) {
      fprintf(logfile, "%s%c", time_name[t].c_str(), SEP); count++;
    }
    assert(count == countTimeInfoEntries);
  } // TIME INFO
  { // VERSION INFO
    fprintf(logfile, "%s%c", "vpc_version", SEP);
    fprintf(logfile, "%s%c", "cbc_version", SEP);
    fprintf(logfile, "%s%c", "clp_version", SEP);
    fprintf(logfile, "%s%c", "gurobi_version", SEP);
    fprintf(logfile, "%s%c", "cplex_version", SEP);
  } // VERSION INFO
  { // WRAP UP INFO
    fprintf(logfile, "%s%c", "ExitReason", SEP);
    fprintf(logfile, "%s%c", "end_time_string", SEP);
    fprintf(logfile, "%s%c", "time elapsed", SEP);
    fprintf(logfile, "%s%c", "instname", SEP);
  } // WRAP UP INFO
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
  fflush(logfile);
} /* printBoundAndGapInfo */

void printSummaryBBInfo(const std::vector<SummaryBBInfo>& info_vec, FILE* logfile,
    const bool print_blanks, const char SEP) {
  if (!logfile)
      return;

  int count = 0;
  for (auto info : info_vec) {
    if (!print_blanks)
      fprintf(logfile, "%ld%c", info.first_bb_info.nodes, SEP);
    else
      fprintf(logfile, "%c", SEP);
    count++;
  }
  for (auto info : info_vec) {
    if (!print_blanks)
      fprintf(logfile, "%ld%c", info.best_bb_info.nodes, SEP);
    else
      fprintf(logfile, "%c", SEP);
    count++;
  }
  for (auto info : info_vec) {
    if (!print_blanks)
      fprintf(logfile, "%2.3f%c", info.first_bb_info.time, SEP);
    else
      fprintf(logfile, "%c", SEP);
    count++;
  }
  for (auto info : info_vec) {
    if (!print_blanks)
      fprintf(logfile, "%2.3f%c", info.best_bb_info.time, SEP);
    else
      fprintf(logfile, "%c", SEP);
    count++;
  }
  fflush(logfile);
  assert(count == countSummaryBBInfoEntries);
} /* printSummaryBBInfo */

void printFullBBInfo(const std::vector<SummaryBBInfo>& info_vec, FILE* logfile,
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
      fprintf(logfile, "%ld%c", info.first_bb_info.root_iters, SEP);
      count++;
    }
    for (auto info : info_vec) {
      fprintf(logfile, "%ld%c", info.best_bb_info.root_iters, SEP);
      count++;
    }
    for (auto info : info_vec) {
      fprintf(logfile, "%ld%c", info.avg_bb_info.root_iters, SEP);
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
    for (auto info : info_vec) {
      std::vector<std::string> vec_str;
      createStringFromBBInfoVec(info.vec_bb_info, vec_str);
      for (unsigned i = 0; i < vec_str.size(); i++) {
        fprintf(logfile, "%s%c", vec_str[i].c_str(), SEP);
        count++;
      }
    }
  } else {
    for (unsigned i = 0; i < BB_INFO_CONTENTS.size() * info_vec.size() * 4; i++) {
      fprintf(logfile, "%c", SEP); count++;
    }
  }
  fflush(logfile);
  assert(count == countFullBBInfoEntries);
} /* printFullBBInfo */

void printOrigProbInfo(const OsiSolverInterface* const solver, FILE* logfile,
    const char SEP) {
  if (!logfile)
    return;

  const int num_rows = solver->getNumRows();
  const int num_cols = solver->getNumCols();
  // Get row stats
  int num_eq_rows = 0, num_bound_rows = 0, num_assign_rows = 0;
  const CoinPackedMatrix* mat = solver->getMatrixByRow();
  for (int row = 0; row < num_rows; row++) {
    const double row_lb = solver->getRowLower()[row];
    const double row_ub = solver->getRowUpper()[row];
    if (isVal(row_lb, row_ub))
      num_eq_rows++;
    if (mat->getVectorSize(row) == 1) {
      if (isVal(row_lb, row_ub))
        num_assign_rows++;
      else
        num_bound_rows++;
    }
  }
  // Calculate fractionality
  int num_frac = 0;
  int num_fixed = 0, num_gen_int = 0, num_bin = 0, num_cont = 0;
  double min_frac = 1., max_frac = 0.;
  for (int col = 0; col < num_cols; col++) {
    const double col_lb = solver->getColLower()[col];
    const double col_ub = solver->getColUpper()[col];
    if (isVal(col_lb, col_ub))
      num_fixed++;
    if (!solver->isInteger(col)) {
      num_cont++;
      continue;
    }
    if (solver->isBinary(col))
      num_bin++;
    else
      num_gen_int++;
    const double val = solver->getColSolution()[col];
    const double floorxk = std::floor(val);
    const double ceilxk = std::ceil(val);
    const double frac = CoinMin(val - floorxk, ceilxk - val);
    if (!isVal(frac, 0., 1e-5)) {
      num_frac++;
      if (frac < min_frac)
        min_frac = frac;
      if (frac > max_frac)
        max_frac = frac;
    }
  }

  int count = 0;
  fprintf(logfile, "%s%c", stringValue(num_rows).c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue(num_cols).c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue(num_frac).c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue(min_frac, "%.5f").c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue(max_frac, "%.5f").c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue(num_eq_rows).c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue(num_bound_rows).c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue(num_assign_rows).c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue(num_fixed).c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue(num_gen_int).c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue(num_bin).c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue(num_cont).c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue((double) mat->getNumElements() / (num_rows * num_cols)).c_str(), SEP); count++;
  fflush(logfile);
  assert(count == countOrigProbEntries);
} /* printOrigProbInfo */

/**
 * @details Assumed that solver is already with cuts added
 */
void printPostCutProbInfo(const OsiSolverInterface* const solver,
    const SummaryCutInfo& cutInfoGMICs, const SummaryCutInfo& cutInfoVPCs,
    FILE* logfile, const char SEP) {
  if (!logfile)
    return;

  const int num_rows = solver->getNumRows();
  const int num_cols = solver->getNumCols();

  // Calculate fractionality
  int num_frac = 0;
  double min_frac = 1., max_frac = 0.;
  for (int col = 0; col < num_cols; col++) {
    if (!solver->isInteger(col)) {
      continue;
    }
    const double val = solver->getColSolution()[col];
    const double floorxk = std::floor(val);
    const double ceilxk = std::ceil(val);
    const double frac = CoinMin(val - floorxk, ceilxk - val);
    if (!isVal(frac, 0., 1e-5)) {
      num_frac++;
      if (frac < min_frac)
        min_frac = frac;
      if (frac > max_frac)
        max_frac = frac;
    }
  }

  int count = 0;
  fprintf(logfile, "%s%c", stringValue(num_frac).c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue(min_frac, "%.5f").c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue(max_frac, "%.5f").c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue((double) solver->getMatrixByCol()->getNumElements() / (num_rows * num_cols)).c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue(cutInfoGMICs.num_active_gmic).c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue(cutInfoVPCs.num_active_gmic).c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue(cutInfoGMICs.num_active_vpc).c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue(cutInfoVPCs.num_active_vpc).c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue(cutInfoGMICs.num_active_all).c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue(cutInfoVPCs.num_active_all).c_str(), SEP); count++;
  fflush(logfile);
  assert(count == countPostCutProbEntries);
} /* printPostCutProbInfo */

void printDisjInfo(const SummaryDisjunctionInfo& disjInfo, FILE* logfile,
    const char SEP) {
  if (!logfile)
    return;

  int count = 0;
  fprintf(logfile, "%s%c", stringValue(disjInfo.avg_num_terms, "%g").c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue(disjInfo.num_integer_sol, "%d").c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue(disjInfo.num_disj, "%d").c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue(disjInfo.avg_density_prlp, "%g").c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue(disjInfo.avg_num_rows_prlp, "%g").c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue(disjInfo.avg_num_cols_prlp, "%g").c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue(disjInfo.avg_num_points_prlp, "%g").c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue(disjInfo.avg_num_rays_prlp, "%g").c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue(disjInfo.avg_explored_nodes, "%g").c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue(disjInfo.avg_pruned_nodes, "%g").c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue(disjInfo.avg_min_depth, "%g").c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue(disjInfo.avg_max_depth, "%g").c_str(), SEP); count++;
  fflush(logfile);
  assert(count == countDisjInfoEntries);
} /* printDisjInfo */

void printCutInfo(const SummaryCutInfo& cutInfoGMICs,
    const SummaryCutInfo& cutInfo, FILE* logfile, const char SEP) {
  if (!logfile)
    return;

  { // CUT INFO
    int count = 0;
    fprintf(logfile, "%s%c", stringValue(cutInfo.num_rounds).c_str(), SEP); count++;
    fprintf(logfile, "%s%c", stringValue(cutInfo.num_cuts).c_str(), SEP); count++;
    fprintf(logfile, "%s%c", stringValue(cutInfo.numCutsOfType[static_cast<int>(CglVPC::CutType::ONE_SIDED_CUT)]).c_str(), SEP); count++;
    fprintf(logfile, "%s%c", stringValue(cutInfo.numCutsOfType[static_cast<int>(CglVPC::CutType::OPTIMALITY_CUT)]).c_str(), SEP); count++;
    if (cutInfo.num_cuts > 0) {
      fprintf(logfile, "%s%c", stringValue(cutInfo.min_support).c_str(), SEP); count++;
      fprintf(logfile, "%s%c", stringValue(cutInfo.max_support).c_str(), SEP); count++;
      fprintf(logfile, "%s%c", stringValue(cutInfo.avg_support, "%.3f").c_str(), SEP); count++;
    } else {
      fprintf(logfile, "%s%c", stringValue(0).c_str(), SEP); count++;
      fprintf(logfile, "%s%c", stringValue(0).c_str(), SEP); count++;
      fprintf(logfile, "%s%c", stringValue(0).c_str(), SEP); count++;
    }
    if (cutInfoGMICs.num_cuts > 0) {
      fprintf(logfile, "%s%c", stringValue(cutInfoGMICs.min_support).c_str(), SEP); count++;
      fprintf(logfile, "%s%c", stringValue(cutInfoGMICs.max_support).c_str(), SEP); count++;
      fprintf(logfile, "%s%c", stringValue(cutInfoGMICs.avg_support, "%.3f").c_str(), SEP); count++;
    } else {
      fprintf(logfile, "%s%c", stringValue(0).c_str(), SEP); count++;
      fprintf(logfile, "%s%c", stringValue(0).c_str(), SEP); count++;
      fprintf(logfile, "%s%c", stringValue(0).c_str(), SEP); count++;
    }
    assert(count == countCutInfoEntries);
  } // CUT INFO
  { // OBJ INFO
    int count = 0;
    fprintf(logfile, "%s%c", stringValue(cutInfo.num_obj_tried).c_str(), SEP); count++;
    for (int obj_ind = 0; obj_ind < static_cast<int>(CglVPC::ObjectiveType::NUM_OBJECTIVE_TYPES); obj_ind++) {
      fprintf(logfile, "%s%c", stringValue(cutInfo.numObjFromHeur[obj_ind]).c_str(), SEP); count++;
      fprintf(logfile, "%s%c", stringValue(cutInfo.numCutsFromHeur[obj_ind]).c_str(), SEP); count++;
      fprintf(logfile, "%s%c", stringValue(cutInfo.numFailsFromHeur[obj_ind]).c_str(), SEP); count++;
      fprintf(logfile, "%s%c", stringValue(cutInfo.numActiveFromHeur[obj_ind]).c_str(), SEP); count++;
    }
    assert(count == countObjInfoEntries);
  } // OBJ INFO
  { // FAIL INFO
    int count = 0;
    fprintf(logfile, "%s%c", stringValue(cutInfo.num_failures).c_str(), SEP); count++;
    for (int fail_ind = 0; fail_ind < static_cast<int>(CglVPC::FailureType::NUM_FAILURE_TYPES); fail_ind++) {
      fprintf(logfile, "%s%c", stringValue(cutInfo.numFails[fail_ind]).c_str(), SEP); count++;
    }
    assert(count == countFailInfoEntries);
  } // FAIL INFO
  fflush(logfile);
} /* printCutInfo */

/// @details Gets cut support size and updates min/max component of \p cutInfo
int checkCutDensity(
    /// [in,out] Where to save min and max support
    SummaryCutInfo& cutInfo,
    /// [in] Row that we want to check
    const OsiRowCut* const cut,
    /// [in] What counts as a zero coefficient
    const double EPS) {
  int num_elem = cut->row().getNumElements();
  const double* el = cut->row().getElements();
  for (int i = 0; i < cut->row().getNumElements(); i++) {
    if (isZero(el[i], EPS)) {
      num_elem--;
    }
  }
  if (num_elem < cutInfo.min_support)
    cutInfo.min_support = num_elem;
  if (num_elem > cutInfo.max_support)
    cutInfo.max_support = num_elem;
  return num_elem;
} // checkCutDensity

/// @brief Find active cuts, and also report density of cuts
bool checkCutActivity(
  const OsiSolverInterface* const solver,
  const OsiRowCut* const cut) {
  if (solver && solver->isProvenOptimal()) {
    const double activity = dotProduct(cut->row(), solver->getColSolution());
    return isVal(activity, cut->rhs());
  } else {
    return false;
  }
} /* checkCutActivity */

/**
 * @details The cut properties we want to look at are:
 * 1. Gap closed
 * 2. Activity (after adding cuts)
 * 3. Density
 */
void analyzeStrength(
    const VPCParameters& params, 
    const OsiSolverInterface* const solver_gmic,
    const OsiSolverInterface* const solver_vpc,
    const OsiSolverInterface* const solver_all,
    SummaryCutInfo& cutInfoGMICs, SummaryCutInfo& cutInfoVPCs,
    const OsiCuts* const gmics, const OsiCuts* const vpcs,
    const SummaryBoundInfo& boundInfo, std::string& output) {
  cutInfoGMICs.num_active_gmic = 0;
  cutInfoGMICs.num_active_vpc = 0;
  cutInfoGMICs.num_active_all = 0;
  cutInfoVPCs.num_active_gmic = 0;
  cutInfoVPCs.num_active_vpc = 0;
  cutInfoVPCs.num_active_all = 0;
  cutInfoVPCs.numActiveFromHeur.resize(static_cast<int>(CglVPC::ObjectiveType::NUM_OBJECTIVE_TYPES), 0);
  if (vpcs) {
    const int num_vpcs = vpcs->sizeCuts();
    int total_support = 0;
    for (int cut_ind = 0; cut_ind < num_vpcs; cut_ind++) {
      const OsiRowCut* const cut = vpcs->rowCutPtr(cut_ind);
      if (checkCutActivity(solver_gmic, cut)) {
        cutInfoVPCs.num_active_gmic++;
      }
      if (checkCutActivity(solver_vpc, cut)) {
        cutInfoVPCs.num_active_vpc++;
        cutInfoVPCs.numActiveFromHeur[static_cast<int>(cutInfoVPCs.objType[cut_ind])]++;
      }
      if (checkCutActivity(solver_all, cut)) {
        cutInfoVPCs.num_active_all++;
      }
      total_support += checkCutDensity(cutInfoVPCs, cut, params.get(EPS) / 2.);
    }
    cutInfoVPCs.avg_support = (double) total_support / num_vpcs;
  }
  if (gmics) {
    const int num_gmics = gmics->sizeCuts();
    cutInfoGMICs.num_cuts = num_gmics;
    int total_support = 0;
    for (int cut_ind = 0; cut_ind < num_gmics; cut_ind++) {
      const OsiRowCut* const cut = gmics->rowCutPtr(cut_ind);
      if (checkCutActivity(solver_gmic, cut)) {
        cutInfoGMICs.num_active_gmic++;
      }
      if (checkCutActivity(solver_vpc, cut)) {
        cutInfoGMICs.num_active_vpc++;
      }
      if (checkCutActivity(solver_all, cut)) {
        cutInfoGMICs.num_active_all++;
      }
      total_support += checkCutDensity(cutInfoGMICs, cut, params.get(EPS) / 2.);
    }
    cutInfoGMICs.avg_support = (double) total_support / num_gmics;
  }

  // Print results from adding cuts
  int NAME_WIDTH = 25;
  int NUM_DIGITS_BEFORE_DEC = 7;
  int NUM_DIGITS_AFTER_DEC = 7;
  const double INF = std::numeric_limits<double>::max();
  char tmpstring[300];

  snprintf(tmpstring, sizeof(tmpstring) / sizeof(char),
      "\n## Results from adding cuts ##\n");
  output += tmpstring;
  snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%-*.*s%s\n",
      NAME_WIDTH, NAME_WIDTH, "LP: ",
      stringValue(boundInfo.lp_obj, "% -*.*g",
        INF,
        NUM_DIGITS_BEFORE_DEC,
        NUM_DIGITS_AFTER_DEC).c_str());
  output += tmpstring;
  if (!isInfinity(std::abs(boundInfo.gmic_obj))) {
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char),
        "%-*.*s%s (%d cuts", NAME_WIDTH, NAME_WIDTH, "GMICs: ",
        stringValue(boundInfo.gmic_obj, "% -*.*g",
          INF,
          NUM_DIGITS_BEFORE_DEC,
          NUM_DIGITS_AFTER_DEC).c_str(),
        boundInfo.num_gmic);
    output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char),
        ", %d active GMICs", cutInfoGMICs.num_active_gmic);
    output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char),
        ", %d active VPCs", cutInfoVPCs.num_active_gmic);
    output += tmpstring;
    output += ")\n";
  }
  if (!isInfinity(std::abs(boundInfo.vpc_obj))) {
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char),
        "%-*.*s%s (%d cuts", NAME_WIDTH, NAME_WIDTH, "VPCs: ",
        stringValue(boundInfo.vpc_obj, "% -*.*g",
          INF,
          NUM_DIGITS_BEFORE_DEC,
          NUM_DIGITS_AFTER_DEC).c_str(),
        boundInfo.num_vpc);
    output += tmpstring;
    if (gmics && gmics->sizeCuts() > 0) {
      snprintf(tmpstring, sizeof(tmpstring) / sizeof(char),
          ", %d active GMICs", cutInfoGMICs.num_active_vpc);
      output += tmpstring;
    }
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char),
        ", %d active VPCs", cutInfoVPCs.num_active_vpc);
    output += tmpstring;
    output += ")\n";
  }
  if (boundInfo.num_gmic + boundInfo.num_lpc > 0) {
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char),
        "%-*.*s%s (%d cuts", NAME_WIDTH, NAME_WIDTH, "All: ",
        stringValue(boundInfo.all_cuts_obj, "% -*.*g", INF, NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC).c_str(), 
        boundInfo.num_gmic + boundInfo.num_lpc + boundInfo.num_vpc);
    output += tmpstring;
    if (gmics && gmics->sizeCuts() > 0) {
      snprintf(tmpstring, sizeof(tmpstring) / sizeof(char),
          ", %d active GMICs", cutInfoGMICs.num_active_all);
      output += tmpstring;
    }
    if (vpcs && vpcs->sizeCuts() > 0) {
      snprintf(tmpstring, sizeof(tmpstring) / sizeof(char),
          ", %d active VPCs", cutInfoVPCs.num_active_all);
      output += tmpstring;
    }
    output += ")\n";
  }
  if (!isInfinity(std::abs(boundInfo.best_disj_obj))) {
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%-*.*s%s\n",
        NAME_WIDTH, NAME_WIDTH, "Disjunctive lb: ",
        stringValue(boundInfo.best_disj_obj, "% -*.*g",
          INF,
          NUM_DIGITS_BEFORE_DEC,
          NUM_DIGITS_AFTER_DEC).c_str());
    output += tmpstring;
  }
  if (!isInfinity(std::abs(boundInfo.worst_disj_obj))) {
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%-*.*s%s\n",
        NAME_WIDTH, NAME_WIDTH, "Disjunctive ub: ",
        stringValue(boundInfo.worst_disj_obj, "% -*.*g",
          INF,
          NUM_DIGITS_BEFORE_DEC,
          NUM_DIGITS_AFTER_DEC).c_str());
    output += tmpstring;
  }
  if (!isInfinity(std::abs(boundInfo.ip_obj))) {
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%-*.*s%s\n",
        NAME_WIDTH, NAME_WIDTH, "IP: ",
        stringValue(boundInfo.ip_obj, "% -*.*g",
          INF,
          NUM_DIGITS_BEFORE_DEC,
          NUM_DIGITS_AFTER_DEC).c_str());
    output += tmpstring;
  }
} /* analyzeStrength */

/** @details Branch-and-bound itself has already been performed */
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
  const int NAME_WIDTH = 10; //25
  const int NUM_DIGITS_BEFORE_DEC = 15; //10
  const int NUM_DIGITS_AFTER_DEC = 2; //2
  const double INF = 1e50; //params.get(doubleParam::INF);
  const bool use_gurobi = use_bb_option(params.get(intParam::BB_STRATEGY), BB_Strategy_Options::gurobi);
  const bool use_cplex = use_bb_option(params.get(intParam::BB_STRATEGY), BB_Strategy_Options::cplex);
  const bool use_cbc = use_bb_option(params.get(intParam::BB_STRATEGY), BB_Strategy_Options::cbc);
  const std::string SOLVER =
    use_gurobi ? "Gur" :
    (use_cplex ? "Cpx" :
     (use_cbc ? "Cbc" : "")
    );
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
    const std::string CURR_NAME = SOLVER;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%-*s", NAME_WIDTH, CURR_NAME.c_str()); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%s", stringValue(info_nocuts.avg_bb_info.obj, "%-*.*f", INF, NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC).c_str()); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%s", stringValue(info_nocuts.avg_bb_info.bound, "%-*.*f", INF, NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC).c_str()); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%s", stringValue(info_nocuts.avg_bb_info.iters, "%-*ld", std::numeric_limits<long double>::max(), NUM_DIGITS_BEFORE_DEC).c_str()); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%s", stringValue(info_nocuts.avg_bb_info.nodes, "%-*ld", std::numeric_limits<long double>::max(), NUM_DIGITS_BEFORE_DEC).c_str()); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%s", stringValue(info_nocuts.avg_bb_info.root_passes, "%-*ld", std::numeric_limits<long double>::max(), NUM_DIGITS_BEFORE_DEC).c_str()); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%s", stringValue(info_nocuts.avg_bb_info.first_cut_pass, "%-*.*f", INF, NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC).c_str()); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%s", stringValue(info_nocuts.avg_bb_info.last_cut_pass, "%-*.*f", INF, NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC).c_str()); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%s", stringValue(info_nocuts.avg_bb_info.root_time, "%-*.*f", INF, NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC).c_str()); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%s", stringValue(info_nocuts.avg_bb_info.last_sol_time, "%-*.*f", INF, NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC).c_str()); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%s", stringValue(info_nocuts.avg_bb_info.time, "%-*.*f", INF, NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC).c_str()); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "\n"); output += tmpstring;
  } // no cuts
  if (branch_with_vpcs) {
    const std::string CURR_NAME = SOLVER + "V";
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%-*s", NAME_WIDTH, CURR_NAME.c_str()); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%s", stringValue(info_mycuts.avg_bb_info.obj, "%-*.*f", INF, NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC).c_str()); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%s", stringValue(info_mycuts.avg_bb_info.bound, "%-*.*f", INF, NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC).c_str()); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%s", stringValue(info_mycuts.avg_bb_info.iters, "%-*ld", std::numeric_limits<long double>::max(), NUM_DIGITS_BEFORE_DEC).c_str()); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%s", stringValue(info_mycuts.avg_bb_info.nodes, "%-*ld", std::numeric_limits<long double>::max(), NUM_DIGITS_BEFORE_DEC).c_str()); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%s", stringValue(info_mycuts.avg_bb_info.root_passes, "%-*ld", std::numeric_limits<long double>::max(), NUM_DIGITS_BEFORE_DEC).c_str()); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%s", stringValue(info_mycuts.avg_bb_info.first_cut_pass, "%-*.*f", INF, NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC).c_str()); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%s", stringValue(info_mycuts.avg_bb_info.last_cut_pass, "%-*.*f", INF, NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC).c_str()); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%s", stringValue(info_mycuts.avg_bb_info.root_time, "%-*.*f", INF, NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC).c_str()); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%s", stringValue(info_mycuts.avg_bb_info.last_sol_time, "%-*.*f", INF, NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC).c_str()); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%s", stringValue(info_mycuts.avg_bb_info.time, "%-*.*f", INF, NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC).c_str()); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "\n"); output += tmpstring;
  } // vpcs
  if (branch_with_gmics) {
    const std::string CURR_NAME = SOLVER + "V+G";
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%-*s", NAME_WIDTH, CURR_NAME.c_str()); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%s", stringValue(info_allcuts.avg_bb_info.obj, "%-*.*f", INF, NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC).c_str()); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%s", stringValue(info_allcuts.avg_bb_info.bound, "%-*.*f", INF, NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC).c_str()); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%s", stringValue(info_allcuts.avg_bb_info.iters, "%-*ld", std::numeric_limits<long double>::max(), NUM_DIGITS_BEFORE_DEC).c_str()); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%s", stringValue(info_allcuts.avg_bb_info.nodes, "%-*ld", std::numeric_limits<long double>::max(), NUM_DIGITS_BEFORE_DEC).c_str()); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%s", stringValue(info_allcuts.avg_bb_info.root_passes, "%-*ld", std::numeric_limits<long double>::max(), NUM_DIGITS_BEFORE_DEC).c_str()); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%s", stringValue(info_allcuts.avg_bb_info.first_cut_pass, "%-*.*f", INF, NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC).c_str()); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%s", stringValue(info_allcuts.avg_bb_info.last_cut_pass, "%-*.*f", INF, NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC).c_str()); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%s", stringValue(info_allcuts.avg_bb_info.root_time, "%-*.*f", INF, NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC).c_str()); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%s", stringValue(info_allcuts.avg_bb_info.last_sol_time, "%-*.*f", INF, NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC).c_str()); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%s", stringValue(info_allcuts.avg_bb_info.time, "%-*.*f", INF, NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC).c_str()); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "\n"); output += tmpstring;
  } // gmics
} /* analyzeBB */

double getNumGomoryRounds(const VPCParameters& params,
    const OsiSolverInterface* const origSolver,
    const OsiSolverInterface* const postCutSolver) {
#ifdef TRACE
  printf("\nGetting number rounds of Gomory cuts req'd to get bound.\n");
#endif
  const int num_cuts = postCutSolver->getNumRows() - origSolver->getNumRows();
  const double post_cut_opt = postCutSolver->getObjValue();
  const int min_sic_rounds = (params.get(STRENGTHEN) == 2) ? 2 : 0;
  int max_rounds = 1000;

  int total_num_sics = 0;
  int num_sic_rounds = 0;
  double curr_sic_opt = 0.;
  std::vector<int> numCutsByRoundSIC;
  std::vector<double> boundByRoundSIC;
  OsiSolverInterface* copySolver = origSolver->clone();
  if (!copySolver->isProvenOptimal()) {
    copySolver->initialSolve();
    checkSolverOptimality(copySolver, false);
  }
  while (num_sic_rounds < min_sic_rounds
      || (lessThanVal(curr_sic_opt, post_cut_opt) && total_num_sics < num_cuts)) {
    OsiCuts GMICs;
    CglGMI gen;
    gen.generateCuts(*copySolver, GMICs);
    const int curr_num_cuts = GMICs.sizeCuts();
    if (curr_num_cuts == 0)
      break;

    num_sic_rounds++;
    total_num_sics += curr_num_cuts;
    numCutsByRoundSIC.push_back(curr_num_cuts);
    curr_sic_opt = applyCutsCustom(copySolver, GMICs, params.logfile);
    boundByRoundSIC.push_back(curr_sic_opt);

    // Other stopping conditions:
    // Bound does not improve at all after one round
    if (num_sic_rounds >= 2
        && !greaterThanVal(curr_sic_opt, boundByRoundSIC[num_sic_rounds - 2])) {
      break;
    }
    // Bound does not significantly improve after five rounds
    if (num_sic_rounds > 4) {
      const double delta = curr_sic_opt - boundByRoundSIC[num_sic_rounds - 4];
      if (!greaterThanVal(delta, 1e-3)) {
        break;
      }
    }
  } // do rounds of Gomory cuts
  if (max_rounds < num_sic_rounds) {
    max_rounds = boundByRoundSIC.size();
  }
  const double final_sic_bound = copySolver->getObjValue();
  return final_sic_bound;
} /* getNumGomoryRounds */

void updateDisjInfo(SummaryDisjunctionInfo& disjInfo, const int num_disj, const CglVPC& gen) {
  if (num_disj <= 0)
    return;
  const Disjunction* const disj = gen.getDisjunction();
  const PRLP* const prlp = gen.getPRLP();
  if (!prlp)
    return;
  disjInfo.num_disj = num_disj;
  disjInfo.num_integer_sol += !(disj->integer_sol.empty());
  disjInfo.avg_num_terms = (disjInfo.avg_num_terms * (num_disj - 1) + disj->num_terms) / num_disj;
  disjInfo.avg_density_prlp = (disjInfo.avg_density_prlp * (num_disj - 1) + prlp->density) / num_disj;
  disjInfo.avg_num_rows_prlp += (disjInfo.avg_num_rows_prlp * (num_disj - 1) + prlp->getNumRows()) / num_disj;
  disjInfo.avg_num_cols_prlp += (disjInfo.avg_num_cols_prlp * (num_disj - 1) + prlp->getNumCols()) / num_disj;
  disjInfo.avg_num_points_prlp += (disjInfo.avg_num_points_prlp * (num_disj - 1) + prlp->numPoints) / num_disj;
  disjInfo.avg_num_rays_prlp += (disjInfo.avg_num_rays_prlp * (num_disj - 1) + prlp->numRays) / num_disj;
  try {
    const PartialBBDisjunction* const partialDisj =
        dynamic_cast<const PartialBBDisjunction* const >(disj);
    disjInfo.avg_explored_nodes += (disjInfo.avg_explored_nodes * (num_disj - 1) + partialDisj->data.num_nodes_on_tree) / num_disj;
    disjInfo.avg_pruned_nodes += (disjInfo.avg_pruned_nodes * (num_disj - 1) + partialDisj->data.num_pruned_nodes) / num_disj;
    disjInfo.avg_min_depth += (disjInfo.avg_min_depth * (num_disj - 1) + partialDisj->data.min_node_depth) / num_disj;
    disjInfo.avg_max_depth += (disjInfo.avg_max_depth * (num_disj - 1) + partialDisj->data.max_node_depth) / num_disj;
  } catch (std::exception& e) {

  }
} /* updateDisjInfo */

/**
 * @details Use this to add to cutInfo (but within one round,
 * because the cutType and objType vectors are cleared in gen in each round
 * (so tracking that based on isSetupForRepeatedUse does not work,
 * and the old cutType and objType stored in cutInfo would be overwritten)
 */
void updateCutInfo(SummaryCutInfo& cutInfo, const CglVPC& gen) {
  cutInfo.num_cuts += gen.num_cuts;
  cutInfo.num_obj_tried += gen.num_obj_tried;
  cutInfo.num_failures += gen.num_failures;

  // For cutType and objType, what we do depends on whether the generator is setup for repeated use or not
  if (gen.isSetupForRepeatedUse) {
    cutInfo.cutType = gen.cutType;
    cutInfo.objType = gen.objType;
  } else {
    cutInfo.cutType.insert(cutInfo.cutType.end(), gen.cutType.begin(), gen.cutType.end());
    cutInfo.objType.insert(cutInfo.objType.end(), gen.objType.begin(), gen.objType.end());
  }

  if (cutInfo.numCutsOfType.size() > 0) {
    for (int i = 0; i < static_cast<int>(CglVPC::CutType::NUM_CUT_TYPES); i++) {
      cutInfo.numCutsOfType[i] += gen.numCutsOfType[i];
    }
  } else {
    cutInfo.numCutsOfType = gen.numCutsOfType;
  }

  if (cutInfo.numCutsFromHeur.size() > 0) {
    for (int i = 0; i < static_cast<int>(CglVPC::ObjectiveType::NUM_OBJECTIVE_TYPES); i++) {
      cutInfo.numCutsFromHeur[i] += gen.numCutsFromHeur[i];
      cutInfo.numObjFromHeur[i] += gen.numObjFromHeur[i];
      cutInfo.numFailsFromHeur[i] += gen.numFailsFromHeur[i];
    }
  } else {
    cutInfo.numCutsFromHeur = gen.numCutsFromHeur;
    cutInfo.numObjFromHeur = gen.numObjFromHeur;
    cutInfo.numFailsFromHeur = gen.numFailsFromHeur;
  }

  if (cutInfo.numFails.size() > 0) {
    for (int i = 0; i < static_cast<int>(CglVPC::FailureType::NUM_FAILURE_TYPES); i++) {
      cutInfo.numFails[i] += gen.numFails[i];
    }
  } else {
    cutInfo.numFails = gen.numFails;
  }
} /* updateCutInfo (within one round) */

/**
 * @details Compute total number of cuts / objectives / failures of various types, as well as total activity
 */
void setCutInfo(SummaryCutInfo& cutInfo, const int num_rounds,
    const SummaryCutInfo* const oldCutInfos) {
  const int numCutTypes = static_cast<int>(CglVPC::CutType::NUM_CUT_TYPES);
  const int numObjTypes = static_cast<int>(CglVPC::ObjectiveType::NUM_OBJECTIVE_TYPES);
  const int numFailTypes = static_cast<int>(CglVPC::FailureType::NUM_FAILURE_TYPES);

  cutInfo.num_cuts = 0;
  cutInfo.num_active_gmic = 0;
  cutInfo.num_active_vpc = 0;
  cutInfo.num_active_all = 0;
  cutInfo.num_obj_tried = 0;
  cutInfo.num_failures = 0;
  cutInfo.num_rounds = num_rounds;
  cutInfo.cutType.resize(0);
  cutInfo.objType.resize(0);
  cutInfo.numCutsOfType.clear();
  cutInfo.numCutsOfType.resize(numCutTypes, 0);
  cutInfo.numCutsFromHeur.clear();
  cutInfo.numCutsFromHeur.resize(numObjTypes, 0);
  cutInfo.numObjFromHeur.clear();
  cutInfo.numObjFromHeur.resize(numObjTypes, 0);
  cutInfo.numFailsFromHeur.clear();
  cutInfo.numFailsFromHeur.resize(numObjTypes, 0);
  cutInfo.numActiveFromHeur.clear();
  cutInfo.numActiveFromHeur.resize(numObjTypes, 0);
  cutInfo.numFails.clear();
  cutInfo.numFails.resize(numFailTypes, 0);

  for (int round = 0; round < num_rounds; round++) {
    cutInfo.num_cuts += oldCutInfos[round].num_cuts;
    cutInfo.num_active_gmic += oldCutInfos[round].num_active_gmic;
    cutInfo.num_active_vpc += oldCutInfos[round].num_active_vpc;
    cutInfo.num_active_all += oldCutInfos[round].num_active_all;
    cutInfo.num_obj_tried += oldCutInfos[round].num_obj_tried;
    cutInfo.num_failures += oldCutInfos[round].num_failures;

    for (int i = 0; i < numCutTypes; i++) {
      cutInfo.numCutsOfType[i] += oldCutInfos[round].numCutsOfType[i];
    }
    for (int i = 0; i < numObjTypes; i++) {
      cutInfo.numCutsFromHeur[i] += oldCutInfos[round].numCutsFromHeur[i];
      cutInfo.numObjFromHeur[i] += oldCutInfos[round].numObjFromHeur[i];
      cutInfo.numFailsFromHeur[i] += oldCutInfos[round].numFailsFromHeur[i];
    }
    if (oldCutInfos[round].numActiveFromHeur.size() > 0) {
      for (int i = 0; i < numObjTypes; i++) {
        cutInfo.numActiveFromHeur[i] += oldCutInfos[round].numActiveFromHeur[i];
      }
    }
    for (int i = 0; i < numFailTypes; i++) {
      cutInfo.numFails[i] += oldCutInfos[round].numFails[i];
    }
  }

  cutInfo.cutType.resize(cutInfo.num_cuts);
  cutInfo.objType.resize(cutInfo.num_cuts);
  int cut_ind = 0;
  for (int round = 0; round < num_rounds; round++) {
    for (int i = 0; i < oldCutInfos[round].num_cuts; i++) {
      cutInfo.cutType[cut_ind] = oldCutInfos[round].cutType[i];
      cutInfo.objType[cut_ind] = oldCutInfos[round].objType[i];
      cut_ind++;
    }
  }
} /* setCutInfo (merge from multiple rounds) */
