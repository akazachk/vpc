#pragma once
#include <map>
#include <string>
#include <vector>
#include <cstdio> // fprintf

enum intParam {
  CUTLIMIT, // max number of cuts generated
  NUM_DISJ_TERMS,
  // Total used to decide the choose:
  // variable decision => hundreds digit: 0: default, 1: default+second criterion, 2: max min change+second (max max change), 3: second-best default, 4: second-best max-min change, -x: -1 * (1+x);
  // branch decision => tens digit: 0: default, 1: dynamic, 2: strong, 3: none;
  // node comparison decision => ones digit: 0: default: 1: bfs, 2: depth, 3: estimate, 4: objective
  PARTIAL_BB_STRATEGY,
  PARTIAL_BB_NUM_STRONG, // -1: num cols, -2: sqrt(num cols), >= 0: that many
  PRLP_BETA, // rhs in nb space; 1: cut away LP opt,  -1: do not cut away LP opt, 0: both
  STRENGTHEN, // 0: no, 1: yes, when possible, 2: same as 1 plus add GMICs to strengthen each disjunctive term
  TEMP, // useful for various temporary parameter changes
  NUM_INT_PARAMS
}; /* intParam */
const std::vector<std::string> intParamName {
  "CUTLIMIT",
  "NUM_DISJ_TERMS",
  "PARTIAL_BB_STRATEGY",
  "PARTIAL_BB_NUM_STRONG",
  "STRENGTHEN",
  "TEMP"
};

enum doubleParam {
  AWAY,
  DIFFEPS,
  EPS,
  RAYEPS,
  PARTIAL_BB_TIMELIMIT,
  PRLP_TIMELIMIT, // -1 = unlimited time for initial solve, and then half that time subsequently; -2 = unlimited time always
  TIMELIMIT,
  NUM_DOUBLE_PARAMS
}; /* doubleParam */
const std::vector<std::string> doubleParamName {
  "AWAY",
  "DIFFEPS",
  "EPS",
  "RAYEPS",
  "PARTIAL_BB_TIMELIMIT",
  "PRLP_TIMELIMIT",
  "TIMELIMIT"
};

enum stringParam {
  FILENAME,
  LOGFILE,
  OPTFILE,
  NUM_STRING_PARAMS
}; /* stringParam */
const std::vector<std::string> stringParamName {
  "FILENAME",
  "LOGFILE",
  "OPTFILE"
};

struct VPCParameters {
  std::map<intParam, int> intParamValues {
    {intParam::CUTLIMIT, 0}, // no limit
    {intParam::NUM_DISJ_TERMS, 0}, // no disjunction (=> no cuts)
    {intParam::PARTIAL_BB_STRATEGY, 4}, // 004 => default var & branch decisions, and choose next node by min objective
    {intParam::PARTIAL_BB_NUM_STRONG, 5}, // consider 5 strong branching candidates
    {intParam::PRLP_BETA, 1}, // cut away the LP optimum
    {intParam::STRENGTHEN, 1}, // strengthen GMICs but not VPCs
    {intParam::TEMP, 0} // do not enable any temporary options
  };
  int get(intParam param) const { return intParamValues.find(param)->second; }
  void set(intParam param, int value) { intParamValues[param] = value; }

  std::map<doubleParam, double> doubleParamValues {
    {doubleParam::AWAY, 1e-3},
    {doubleParam::DIFFEPS, 1e-3}, // to check whether something is different enough to throw an error
    {doubleParam::EPS, 1e-7},
    {doubleParam::RAYEPS, 1e-7},
    {doubleParam::PARTIAL_BB_TIMELIMIT, 3600},
    {doubleParam::PRLP_TIMELIMIT, -1},
    {doubleParam::TIMELIMIT, 60}
  };
  double get(doubleParam param) const { return doubleParamValues.find(param)->second; }
  void set(doubleParam param, double value) {
    doubleParamValues[param] = value;
  }

  std::map<stringParam, std::string> stringParamValues {
    {stringParam::FILENAME, ""},
    {stringParam::LOGFILE, ""},
    {stringParam::OPTFILE, ""}
  };
  std::string get(stringParam param) const { return stringParamValues.find(param)->second; }
  void set(stringParam param, std::string value) {
    stringParamValues[param] = value;
  }

  FILE* logfile = NULL;
}; /* struct VPCParameters */

inline void printParams(VPCParameters param, FILE* outfile = stdout) {
  for (int i = 0; i < intParam::NUM_INT_PARAMS; i++) {
    fprintf(outfile, "%s,%d\n", intParamName[i].c_str(), param.get(static_cast<intParam>(i)));
  }
  for (int i = 0; i < doubleParam::NUM_DOUBLE_PARAMS; i++) {
    fprintf(outfile, "%s,%e\n", doubleParamName[i].c_str(), param.get(static_cast<doubleParam>(i)));
  }
  for (int i = 0; i < stringParam::NUM_STRING_PARAMS; i++) {
    fprintf(outfile, "%s,%s\n", stringParamName[i].c_str(), param.get(static_cast<stringParam>(i)).c_str());
  }
  fflush(outfile);
} /* printParams */
