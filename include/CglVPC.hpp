// Name:     CglVPC.hpp
// Author:   A. M. Kazachkov
// Date:     2018-Dec-24
//-----------------------------------------------------------------------------

#pragma once

#include <vector>
#include <string>

// COIN-OR files
#include <CglCutGenerator.hpp>

// Project files
#include "VPCEventHandler.hpp"
#include "VPCParameters.hpp" // SolverInterface and VPCParameters
#include "TimeStats.hpp"

class Disjunction; // include is in source file
class PartialBBDisjunction;

enum class ExitReason {
  SUCCESS_EXIT = 0,
  CUT_LIMIT_EXIT,
  FAIL_LIMIT_EXIT,
  PARTIAL_BB_OPTIMAL_SOLUTION_FOUND_EXIT,
  TIME_LIMIT_EXIT,
  NO_DISJUNCTION_EXIT,
  UNKNOWN,
  NUM_EXIT_REASONS
};
/* ExitReason */
const std::vector<std::string> ExitReasonName {
  "SUCCESS",
  "CUT_LIMIT",
  "FAIL_LIMIT",
  "PARTIAL_BB_INTEGER_SOLUTION_FOUND",
  "TIME_LIMIT",
  "NO_DISJUNCTION",
  "UNKNOWN"
}; /* ExitReasonName */

class CglVPC : public CglCutGenerator {
public:
  friend class PRLP;

  // Enums
  // If adding to any of these enums, do not forget to add the appropriate name (in the right position) in *Name in the source file
  enum class VPCMode {
    PARTIAL_BB,
    SPLITS,
    CROSSES,
    NUM_VPC_MODES
  };

  enum class VPCTimeStats {
    TOTAL_TIME,
    INIT_SOLVE_TIME,
    DISJ_SETUP_TIME,
    DISJ_GEN_TIME,
    PRLP_SETUP_TIME,
    PRLP_SOLVE_TIME,
    GEN_CUTS_TIME,
    NUM_TIME_STATS
  }; /* VPCTimeStats */

  enum class CutType {
    ONE_SIDED_CUT,
    OPTIMALITY_CUT,
    VPC,
    NUM_CUT_TYPES
  }; /* CutType */

  enum class CutHeuristics {
    DUMMY_OBJ,
    ALL_ONES,
    CUT_VERTICES,
    ITER_BILINEAR,
    UNIT_VECTORS,
    DISJ_LB,
    TIGHT_POINTS,
    TIGHT_RAYS,
    TIGHT_POINTS2,
    TIGHT_RAYS2,
    ONE_SIDED,
    NUM_CUT_HEUR
  };

  enum class FailureType {
    ABANDONED,
    BAD_DYNAMISM,
    BAD_SUPPORT,
    BAD_VIOLATION,
    CUT_LIMIT,
    DUAL_INFEASIBLE,
    DUPLICATE_SIC,
    DUPLICATE_VPC,
    ITERATION_LIMIT,
    ORTHOGONALITY_SIC,
    ORTHOGONALITY_VPC,
    PRIMAL_INFEASIBLE,
    TIME_LIMIT,
    NUMERICAL_ISSUES_WARNING,
    PRIMAL_INFEASIBLE_NO_OBJ,
    NUMERICAL_ISSUES_NO_OBJ,
    UNKNOWN,
    NUM_FAILURES
  }; /* FailureType */

  // Static variables
  static const std::vector<std::string> VPCTimeStatsName;
  static const std::vector<std::string> CutTypeName;
  static const std::vector<std::string> CutHeuristicsName;
  static const std::vector<std::string> FailureTypeName;
  static const std::string time_T1; // = "TIME_TYPE1_";
  static const std::string time_T2; // = "TIME_TYPE2_";

  // Class variables
  VPCParameters params;
  VPCMode mode;
//  Disjunction* disj;
  PartialBBDisjunction* disj;
  ExitReason exitReason;
  TimeStats timer;

  std::vector<CutType> cutType; // one entry per cut
  std::vector<CutHeuristics> cutHeurVec; // one entry per cut

  std::vector<int> numCutsOfType;
  std::vector<int> numObjFromHeur, numCutsFromHeur;
  std::vector<int> numFails;

  double ip_opt;
  int num_cgs, num_cgs_actually_used, num_cgs_leading_to_cuts;
  int num_cuts;
  int num_obj_tried;

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

  /** setParams based on VPCParameters */
  void setParams(const VPCParameters& param);

  /** generateCuts */
  virtual void generateCuts(const OsiSolverInterface&, OsiCuts&, const CglTreeInfo = CglTreeInfo());

protected:
  struct ProblemData {
    int num_cols;
    double lp_opt;
    double minAbsCoeff, maxAbsCoeff; // useful for dynamism calculation
    double EPS;
    std::vector<int> NBVarIndex, varBasicInRow;
    std::vector<int> rowOfVar; // row in which var is basic; if var is nonbasic, then it is (-(1 + nonbasic index))
    std::vector<int> rowOfOrigNBVar; // row in which "original" nonbasic vars are basic; if var is still nonbasic, then it is (-(1 + new nonbasic index))
    std::vector<double> NBReducedCost;

    /** Get index of variable in NBVarIndex, and -1 if it is basic */
    int getVarNBIndex(const int var) const {
      if (rowOfVar[var] <= -1) {
        return -1 * (1 + rowOfVar[var]);
      } else {
        return -1;
      }
    }
  } probData;

  struct PRLPData {
    std::vector<CoinPackedVector> constraints;
    std::vector<double> rhs;
    std::vector<int> term;
    std::vector<double> objViolation;
    void addConstraint(const CoinPackedVector& vec, const double rhs, const int term, const double viol) {
      this->constraints.push_back(vec);
      this->rhs.push_back(rhs);
      this->term.push_back(term);
      this->objViolation.push_back(viol);
    }
  } prlpData;

  void initialize(const CglVPC* const source = NULL, const VPCParameters* const param = NULL);
  void getProblemData(SolverInterface* const solver, ProblemData& probData,
      const ProblemData* const origProbData = NULL,
      const bool enable_factorization = true);

//  ExitReason prepareDisjunction(SolverInterface* const solver, OsiCuts& cuts);

  ExitReason setupConstraints(const SolverInterface* const si, OsiCuts& cuts);
  bool setupDisjunctiveTerm(const int term_ind, const int branching_index,
      const int branching_variable, const int branching_way,
      const double branching_value, std::vector<std::vector<int> >& termIndices,
      std::vector<std::vector<double> >& termCoeff,
      std::vector<double>& termRHS, const SolverInterface* const vpcsolver,
      const SolverInterface* const tmpSolverBase);
  void genDepth1PRCollection(const SolverInterface* const vpcsolver,
      const SolverInterface* const tmpSolver, const ProblemData& origProbData,
      const ProblemData& tmpProbData, const int term_ind);

  ExitReason tryObjectives(OsiCuts& cuts,
      const OsiSolverInterface* const origSolver, const OsiCuts* const structSICs,
      const std::string& timeName);
  void addCut(const OsiRowCut& cut, const CutType& type, OsiCuts& cuts);

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
//      const double time,
      const double few_cuts_fail_threshold = 0.95,
      const double many_cuts_fail_threshold = 0.90,
      const double many_obj_fail_threshold = 0.80,
      const double time_fail_threshold = 0.66) const;

  inline void finish(ExitReason exitReason = ExitReason::UNKNOWN) {
    this->exitReason = exitReason;
    this->timer.end_all();
#ifdef TRACE
    printf("CglVPC: Finishing with exit reason: %s. Number cuts: %d.\n", ExitReasonName[static_cast<int>(exitReason)].c_str(), num_cuts);
#endif
  }

  /** Set/update name of cut generating set (disjunction) */
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
      const std::vector<double>& termRHS, const bool append = false) const {
    if (num_ineq_per_term == 0) {
      return;
    }
    if (!cgsName.empty()) {
      if (!append) {
        cgsName += " V ";
      } else {
        cgsName.resize(cgsName.size() - 1); // remove last ")"
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
