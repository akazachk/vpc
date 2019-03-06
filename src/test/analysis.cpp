// Name:     analysis.cpp
// Author:   A. M. Kazachkov
// Date:     2019-Mar-04
//-----------------------------------------------------------------------------
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
#include "utility.hpp" // isInfinity, stringValue

const int countBoundInfoEntries = 11;
const int countGapInfoEntries = 4;
const int countBBInfoEntries = static_cast<int>(BB_INFO_CONTENTS.size()) * 4 * 2;
const int countOrigProbEntries = 13;
const int countPostCutProbEntries = 6;
const int countDisjInfoEntries = 10;
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
    fprintf(logfile, "%s%c", "VPC+GMIC BOUND", SEP); count++; // 11
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
    fprintf(logfile, "%s%c", "ACTIVE GMIC", SEP); count++;
    fprintf(logfile, "%s%c", "ACTIVE VPC", SEP); count++;
    assert(count == countPostCutProbEntries);
  } // POST-CUT PROB
  { // DISJ INFO
    int count = 0;
    fprintf(logfile, "%s%c", "NUM DISJ TERMS", SEP); count++;
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
  fflush(logfile);
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
  fflush(logfile);
  assert(count == countBBInfoEntries);
} /* printBBInfo */

void printOrigProbInfo(const OsiSolverInterface* const solver, FILE* logfile,
    const char SEP) {
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
 * Assumed that solver is already with cuts added
 */
void printPostCutProbInfo(const OsiSolverInterface* const solver,
    const OsiCuts* const vpcs, const OsiCuts* const gmics, FILE* logfile,
    const char SEP) {
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

  int num_active_vpcs = 0, num_active_gmics = 0;
  if (vpcs) {
    const int num_vpcs = vpcs->sizeCuts();
    for (int cut_ind = 0; cut_ind < num_vpcs; cut_ind++) {
      const double activity = dotProduct(vpcs->rowCutPtr(cut_ind)->row(),
          solver->getColSolution());
      if (isVal(activity, vpcs->rowCutPtr(cut_ind)->rhs())) {
        num_active_vpcs++;
      }
    }
  }
  if (gmics) {
    const int num_gmics = gmics->sizeCuts();
    for (int cut_ind = 0; cut_ind < num_gmics; cut_ind++) {
      const double activity = dotProduct(gmics->rowCutPtr(cut_ind)->row(),
          solver->getColSolution());
      if (isVal(activity, gmics->rowCutPtr(cut_ind)->rhs())) {
        num_active_gmics++;
      }
    }
  }

  int count = 0;
  fprintf(logfile, "%s%c", stringValue(num_frac).c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue(min_frac, "%.5f").c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue(max_frac, "%.5f").c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue((double) solver->getMatrixByCol()->getNumElements() / (num_rows * num_cols)).c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue(num_active_gmics).c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue(num_active_vpcs).c_str(), SEP); count++;
  fflush(logfile);
  assert(count == countPostCutProbEntries);
} /* printPostCutProbInfo */

void printDisjInfo(const SummaryDisjunctionInfo& disjInfo, FILE* logfile,
    const char SEP) {
  if (!logfile)
    return;

  int count = 0;
  fprintf(logfile, "%s%c", stringValue(disjInfo.avg_num_terms, "%.0f").c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue(disjInfo.avg_density_prlp, "%.0f").c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue(disjInfo.avg_num_rows_prlp, "%.0f").c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue(disjInfo.avg_num_cols_prlp, "%.0f").c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue(disjInfo.avg_num_points_prlp, "%.0f").c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue(disjInfo.avg_num_rays_prlp, "%.0f").c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue(disjInfo.avg_explored_nodes, "%.0f").c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue(disjInfo.avg_pruned_nodes, "%.0f").c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue(disjInfo.avg_min_depth, "%.0f").c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue(disjInfo.avg_max_depth, "%.0f").c_str(), SEP); count++;
  fflush(logfile);
  assert(count == countDisjInfoEntries);
}

/**
 * The cut properties we want to look at are:
 * 1. Gap closed
 * 2. Activity (after adding cuts)
 * 3. Density
 */
void analyzeStrength(const VPCParameters& params,
//    const OsiSolverInterface* const origSolver, const OsiCuts* const vpc,
    const SummaryBoundInfo& boundInfo,
    std::string& output) {

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

void getNumGomoryRounds(const VPCParameters& params,
    const OsiSolverInterface* const origSolver,
    const OsiSolverInterface* const postCutSolver) {
  // Get number rounds of SICs needed to meet bound from GICs+SICs
#ifdef TRACE
  printf("\nGetting number rounds of Gomory cuts req'd to get bound.\n");
#endif
  const int num_cuts = postCutSolver->getNumRows() - origSolver->getNumRows();
  const double post_cut_opt = postCutSolver->getObjValue();
  const int min_sic_rounds = (params.get(STRENGTHEN) == 2) ? 2 : 0;
  int max_rounds = 1000;

  int total_num_sics = 0;
  double final_sic_bound = 0.;
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
      return;

    num_sic_rounds++;
    total_num_sics += curr_num_cuts;
    numCutsByRoundSIC.push_back(curr_num_cuts);
    curr_sic_opt = applyCutsCustom(copySolver, GMICs);
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
  final_sic_bound = copySolver->getObjValue();
  if (max_rounds < num_sic_rounds) {
    max_rounds = boundByRoundSIC.size();
  }
} /* getNumGomoryRounds */

void updateDisjInfo(SummaryDisjunctionInfo& disjInfo, const int num_disj, const CglVPC& gen) {
  if (num_disj <= 0)
    return;
  const Disjunction* const disj = gen.getDisjunction();
  const PRLP* const prlp = gen.getPRLP();
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
