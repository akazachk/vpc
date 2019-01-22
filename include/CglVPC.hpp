// Name:     CglVPC.hpp
// Author:   A. M. Kazachkov
// Date:     2018-Dec-24
//-----------------------------------------------------------------------------

#pragma once

//#include <chrono>
#include <vector>
#include <string>

// COIN-OR files
#include <CglCutGenerator.hpp>
#include <OsiClpSolverInterface.hpp>

// Project files
#include "params.hpp"
#include "timestats.hpp"

class CglVPC : public CglCutGenerator {
public:
	// Enums and static variables
  enum ExitReason {
    SUCCESS_EXIT = 0,
    CUT_LIMIT_EXIT,
    FAIL_LIMIT_EXIT,
    PARTIAL_BB_INTEGER_SOLUTION_FOUND_EXIT,
    TIME_LIMIT_EXIT,
    UNKNOWN,
    NUM_EXIT_REASONS
  }; /* ExitReason */
  const std::vector<std::string> ExitReasonName {
    "SUCCESS",
    "CUT_LIMIT",
    "FAIL_LIMIT"
    "PARTIAL_BB_INTEGER_SOLUTION_FOUND",
    "TIME_LIMIT",
    "UNKNOWN"
  }; /* ExitReasonName */

  enum VPCTimeStats {
    TOTAL_TIME,
    INIT_SOLVE_TIME,
    DISJ_SETUP_TIME,
    GEN_TREE_TIME,
    PRLP_SETUP_TIME,
    PRLP_SOLVE_TIME,
    GEN_CUTS_TIME,
    NUM_TIME_STATS
  }; /* VPCTimeStats */
  const std::vector<std::string> VPCTimeStatsName {
    "TOTAL_TIME",
    "INIT_SOLVE_TIME",
    "DISJ_SETUP_TIME",
    "GEN_TREE_TIME",
    "PRLP_SETUP_TIME",
    "PRLP_SOLVE_TIME",
    "GEN_CUTS_TIME"
  }; /* VPCTimeStatsName */

	enum CutType {
		ONE_SIDED_CUT, VPC, NUM_CUT_TYPE_STATS
	}; /* CutType */
	const std::vector<std::string> CutTypeName {
		"ONE_SIDED_CUT", "VPC"
	}; /* CutTypeName */

	enum FailureType {
		NUMERICAL_ISSUES_WARNING, NUM_FAILURES
	}; /* FailureType */

	// Class variables
  VPCParameters params;
  ExitReason exitReason;
  TimeStats timer;

  std::vector<CutType> cutType;
  std::vector<int> numCutsOfType;
  std::vector<int> numFails;

  double branching_lb, branching_ub, min_nb_obj_val, ip_opt;
  int num_cgs, num_cuts;
  int num_disj_terms, num_partial_bb_nodes, num_pruned_nodes, num_obj_tried;

  /** Default constructor */
  CglVPC();

  /** Param constructor */
  CglVPC(const VPCParameters& param);

  /** Copy constructor */
  CglVPC(const CglVPC& source);

  /** Destructor */
  ~CglVPC();

  /** Assignment operator */
  CglVPC& operator=(const CglVPC& source);

  /** Clone */
  virtual CglCutGenerator* clone() const;

  /** setParams based on CglVPCParams */
  void setParams(const VPCParameters& param);

  /** generateCuts */
  virtual void generateCuts(const OsiSolverInterface&, OsiCuts&, const CglTreeInfo = CglTreeInfo());

protected:
	struct ProblemData {
		double EPS;
//		std::vector<int> cstat, rstat;
	  std::vector<int> NBVarIndex, rowOfVar, varBasicInRow;
	  std::vector<int> fractionalCore;
	  std::vector<double> NBReducedCost;
	} probData;

	void getProblemData(OsiSolverInterface* const solver, ProblemData& probData);

  ExitReason prepareDisjunction(OsiClpSolverInterface* solver, OsiCuts& cuts);
  void setupPRLP(const OsiSolverInterface& si);
  void tryObjectives(const OsiSolverInterface& PRLP);
  void addCut(const OsiRowCut& cut, const CutType& type, OsiCuts& cuts);

  void initialize(const CglVPC* const source = NULL, const VPCParameters* const param = NULL);

  int getCutLimit() const;
  inline bool reachedCutLimit(const int num_cuts) const {
    const bool reached_limit = (num_cuts >= getCutLimit());
//    if (reached_limit)
//      exitReason = ExitReason::CUT_LIMIT_EXIT;
    return reached_limit;
  } /* reachedCutLimit */
  inline bool reachedTimeLimit(const std::string& timeName, const double max_time) const {
    const bool reached_limit = (timer.get_total_time(timeName) > max_time);
//    if (reached_limit)
//      exitReason = ExitReason::TIME_LIMIT_EXIT;
    return reached_limit;
  } /* reachedTimeLimit */
//  inline bool reachedTimelimit(const std::chrono::time_point<std::chrono::high_resolution_clock>& start_chrono) const {
//    std::chrono::duration<double> elapsed_seconds = std::chrono::high_resolution_clock::now() - start_chrono;
//    const bool reached_limit = (elapsed_seconds.count() > this->param.get(TIMELIMIT));
////    if (reached_limit)
////      exitReason = ExitReason::TIME_LIMIT_EXIT;
//    return reached_limit;
//  }
  bool reachedFailureLimit(const int num_cuts, const int num_fails,
      //const double time,
      const double few_cuts_fail_threshold = 0.95,
      const double many_cuts_fail_threshold = 0.90,
      const double many_obj_fail_threshold = 0.80,
      const double time_fail_threshold = 0.66) const;

  inline void finish(ExitReason exitReason = ExitReason::UNKNOWN) {
    this->exitReason = exitReason;
    this->timer.end_all();
  }

  inline void setCgsName(std::string& cgsName, const std::string& disjTermName) const {
    if (disjTermName.empty()) {
      return;
    }
    if (!cgsName.empty()) {
      cgsName += " V ";
    }
    cgsName += "(";
    cgsName += disjTermName;
    cgsName += ")";
  } /* setCgsName (given disj term name) */

  inline void setCgsName(std::string& cgsName, const int num_ineq_per_term,
      const std::vector<std::vector<int> >& termIndices,
      const std::vector<std::vector<double> >& termCoeff,
      const std::vector<double>& termRHS, const bool append) const {
    if (num_ineq_per_term == 0) {
      return;
    }
    if (!cgsName.empty()) {
      if (!append) {
        cgsName += " V ";
      } else {
        cgsName.resize(cgsName.size() - 1);
        cgsName += "; ";
      }
    }
    cgsName += append ? "" : "(";
    for (int i = 0; i < num_ineq_per_term; i++) {
      if (i > 0) {
        cgsName += "; ";
      }
      for (int coeff_ind = 0; coeff_ind < (int) termIndices[i].size();
          coeff_ind++) {
        cgsName += (termCoeff[i][coeff_ind] > 0) ? "+" : "-";
        cgsName += "x";
        cgsName += std::to_string(termIndices[i][coeff_ind]);
      }
      cgsName += " >= ";
      cgsName += std::to_string((int) termRHS[i]);
    }
    cgsName += ")";
  } /* setCgsName */
}; /* CglVPC */
