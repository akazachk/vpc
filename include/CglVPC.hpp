/**
 * @file CglVPC.hpp
 * @brief File containing CglVPC class and related methods
 *
 * @author A. M. Kazachkov
 * @date 2018-12-24
 */
#pragma once

#include <vector>
#include <string>

// COIN-OR files
#include <CglCutGenerator.hpp>

// Project files
#include "VPCParameters.hpp" // VPCParameters has to be declared here because of the class variable params
#include "TimeStats.hpp"

class Disjunction; // include is in source file
class PRLP;
enum class DisjExitReason;

/// @brief CglVPC Class, implemented as a CglCutGenerator
/// @details Takes a Disjunction as an input and generates cuts from it based on the relaxed V-polyhedral strategy.
///
/// Stops when exit criterion has been reached:
/// either we have exhausted our objectives,
/// or we have reachedTimeLimit(), reachedCutLimit(), or reachedFailureLimit(),
/// or we may be exiting due to numerical issues;
/// the ultimate exiting reason is stored in #exitReason.
class CglVPC : public CglCutGenerator {
public:
  friend class PRLP;

  ///@{
  /// @name Class enums and structs
  /// If adding to enum XXX, do not forget to add the appropriate name (in the right position) in the vector<string> XXXName in the source file

  /// @brief Documenting the exit status after / whether cuts are generated
  enum class ExitReason {
    SUCCESS_EXIT = 0,            ///< everything seems okay
    CUT_LIMIT_EXIT,              ///< reachedCutLimit() is true
    FAIL_LIMIT_EXIT,             ///< reachedFailureLimit() is true
    OPTIMAL_SOLUTION_FOUND_EXIT, ///< found integer-optimal solution
    PRLP_INFEASIBLE_EXIT,        ///< PRLP is primal infeasible
    PRLP_NUMERICAL_ISSUES_EXIT,  ///< encountered numerical issues trying to solve PRLP (TODO expand)
    PRLP_TIME_LIMIT_EXIT,        ///< could not solve PRLP within time limit
    TIME_LIMIT_EXIT,             ///< reachedTimeLimit() is true
    TOO_FEW_TERMS_EXIT,          ///< e.g., only one disjunctive term is feasible
    NO_CUTS_LIKELY_EXIT,         ///< e.g., when all leaf nodes have same value as LP solution
    NO_DISJUNCTION_EXIT,         ///< inital solution is integer-feasible
    UNKNOWN,                     ///< TODO expand
    NUM_EXIT_REASONS             ///< number of possible exit reasons
  }; /* ExitReason */

  /// @brief Mode in which VPCs will be generated (disjunction type, in fact)
  enum class VPCMode {
    PARTIAL_BB,   ///< generate partial branch-and-bound tree #PartialBBDisjunction
    SPLITS,       ///< set of split disjunctions
    CROSSES,      ///< set of cross disjunctions
    CUSTOM,       ///< user-defined disjunction
    NUM_VPC_MODES ///< number of vpc modes
  }; /* VPCMode */

  /// @brief Timing
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

  /// @brief Type of cuts generated
  enum class CutType {
    ONE_SIDED_CUT,  ///< cuts from disjunctions in which one side is empty
    OPTIMALITY_CUT, ///< objective cut added when optimal solution is found
    VPC,            ///< generic vpcs
    NUM_CUT_TYPES   ///< number of cut types
  }; /* CutType */

  /// @brief Types of objectives used to generate the cuts
  enum class ObjectiveType {
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
    USER,
    OBJ_CUT,
    ONE_SIDED,
    NUM_OBJECTIVE_TYPES ///< number of objective types
  }; /* ObjectiveType */

  /// @brief Types of failures that occurred during cut generation
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
    DLB_EQUALS_DUB_NO_OBJ,
    DLB_EQUALS_LPOPT_NO_OBJ,
    PRIMAL_INFEASIBLE_NO_OBJ,
    NUMERICAL_ISSUES_NO_OBJ,
    UNKNOWN,
    NUM_FAILURE_TYPES ///< number of failure types
  }; /* FailureType */

  /// Data used as input into PRLP construction
  struct PRLPData {
    std::vector<CoinPackedVector> constraints; ///< PRLP constraints
    std::vector<double> rhs; ///< vector of constant sides (always assumed >=)
    std::vector<int> term; ///< which term the point/ray belongs to
    std::vector<double> objViolation; ///< how much the corresponding point/ray improves on the original objective value

    /// Add constraint, rhs, index of disjunctive term, and how much this point violates the objective inequality wrt the LP optimum
    void addConstraint(const CoinPackedVector& vec, const double rhs, const int term, const double viol) {
      this->constraints.push_back(vec);
      this->rhs.push_back(rhs);
      this->term.push_back(term);
      this->objViolation.push_back(viol);
    }
  }; /* PRLPData */

  /// @brief Quick reference to information from current instance
  /// @details This information is key for converting to / from nonbasic space
  struct ProblemData {
    int num_cols; ///< number of variables
    double lp_opt; ///< value of LP optimal solution
    double minAbsCoeff; ///< useful for dynamism calculation
    double maxAbsCoeff; ///< useful for dynamism calculation
    double EPS; ///< epsilon computed from this particular instance data
    std::vector<int> NBVarIndex; /// indices of nonbasic vectors at LP optimum
    std::vector<int> varBasicInRow; /// for the optimal basis, what are the variables basic in each row
    std::vector<int> rowOfVar; ///< row in which var is basic; if var is nonbasic, then it is (-(1 + nonbasic index))
    std::vector<int> rowOfOrigNBVar; ///< row in which "original" nonbasic vars are basic; if var is still nonbasic, then it is (-(1 + new nonbasic index))
    std::vector<double> NBReducedCost; ///< objective coefficients in the nonbasic space

    /// Get index of variable in NBVarIndex, and -1 if it is basic
    int getVarNBIndex(const int var) const {
      if (rowOfVar[var] <= -1) {
        return -1 * (1 + rowOfVar[var]);
      } else {
        return -1;
      }
    }
  } probData; ///< Quick reference to key information from current instance

  ///@} // class enums and structs

  ///@{
  /// @name Static variables/functions
  static const std::vector<std::string> ExitReasonName;    ///< #ExitReason strings
  static const std::vector<std::string> VPCModeName;       ///< #VPCMode strings
  static const std::vector<std::string> VPCTimeStatsName;  ///< #VPCTimeStats strings
  static const std::vector<std::string> CutTypeName;       ///< #CutType strings
  static const std::vector<std::string> ObjectiveTypeName; ///< #ObjectiveType strings
  static const std::vector<std::string> FailureTypeName;   ///< #FailureType strings
  static const std::string time_T1; ///< = "TIME_TYPE1_"
  static const std::string time_T2; ///< = "TIME_TYPE2_"

  /// @brief Universal way to check whether we reached the limit for the number of cuts for each split
  static int getCutLimit(const int CUTLIMIT, const int numFracVar);
  ///@} // static vars/funcs

  ///@{
  /// @name Class variables
  VPCParametersNamespace::VPCParameters params; ///< store parameters for this specific run
  VPCMode mode;                                 ///< type of Disjunction to use
  ExitReason exitReason;                        ///< track why CglVPC instance exits
  TimeStats timer;                              ///< holds time statistics for all timers enumerated in #VPCTimeStats
  PRLP* prlp = NULL;                            ///< pointer to PRLP instance constructed from PRLPData

  std::vector<CutType> cutType;       ///< one entry per cut
  std::vector<ObjectiveType> objType; ///< one entry per cut

  std::vector<int> numCutsOfType;     ///< one entry per cut type
  std::vector<int> numCutsFromHeur;   ///< one entry per objective type
  std::vector<int> numObjFromHeur;    ///< one entry per objective type
  std::vector<int> numFailsFromHeur;  ///< one entry per objective type
  std::vector<int> numFails;          ///< one entry per failure type

  int init_num_cuts; ///< useful for keeping track of # of old cuts when #isSetupForRepeatedUse == true
  int num_rounds;    ///< which round of cuts is being generated (can potentially inform parameter choices)
  int num_cuts;      ///< this is always the number of cuts generated by the current call to generateCuts(const OsiSolverInterface&, OsiCuts&, const CglTreeInfo)
  int num_obj_tried; ///< number of times PRLP was solved (except initial solve with no objective) = #num_cuts + #num_failures
  int num_failures;  ///< number of times PRLP was solved and no new cut was added

  /// Determines whether user or class deletes disjunction; default is false
  bool ownsDisjunction;

  /// @brief When setup for repeated use, information about previous cuts is not deleted
  /// @details This is useful, e.g., when generating cuts from multiple split disjunctions
  /// (the user does not want to regenerate cuts from previous splits,
  /// and, moreover, would want to replace old dominated cuts from different splits)
  bool isSetupForRepeatedUse; // default is false

  /// @brief When we find a VPC that is stronger than an existing cut from before this run, can we replace it? Default is false 
  /// @details It can be the case that the instance is setup for repeated use, but we do not wish to replace any given cuts (merely prevent duplicating them).
  /// On the other hand, if the instance is not setup for repeated use, we require that canReplaceGivenCuts is false.
  /// (Though the user can still give a set of old cuts to prevent duplicating them)
  bool canReplaceGivenCuts;
  ///@} // class vars

  /// Default constructor
  CglVPC();

  /// Param constructor
  CglVPC(const VPCParametersNamespace::VPCParameters& param, const int round_ind = 0);

  /// Copy constructor
  CglVPC(const CglVPC& source);

  /// Destructor
  ~CglVPC();

  /// Assignment operator
  CglVPC& operator=(const CglVPC& source);

  /// Clone
  virtual CglCutGenerator* clone() const;

  /// @brief Prepare for repeated use (e.g., multiple disjunctions in one round)
  /// @details When setup for repeated use, information about previous cuts is not deleted
  /// This is useful, e.g., when generating cuts from multiple split disjunctions
  /// (the user does not want to regenerate cuts from previous splits,
  /// and, moreover, would want to replace old dominated cuts from different splits)
  inline void setupRepeatedUse(bool useRepeatedly = true) {
    isSetupForRepeatedUse = useRepeatedly;
  }

  ///@{
  /// @name Get/set methods

  /// Set #params based on VPCParameters
  void setParams(const VPCParametersNamespace::VPCParameters& param);

  /// Get the cut limit for this generator
  int getCutLimit() const;

  /// Disjunction being used by the class
  inline Disjunction* const disj() const { return this->disjunction; }
  /// alias for disj()
  inline Disjunction* const getDisjunction() const { return this->disj(); } // alias
  /// Specify a #disjunction to use and whether the CglVPC class #ownsDisjunction (and is responsible for freeing memory)
  void setDisjunction(Disjunction* const sourceDisj, int ownIt = -1);

  /// Return #prlpData 
  inline const PRLPData& getPRLPData() const { return this->prlpData; }
  /// Return #prlp
  inline const PRLP* const getPRLP() const { return this->prlp; }

  /// Return #user_objectives
  inline std::vector<std::vector<double> > getUserObjectives() const { return this->user_objectives; }
  /// @brief User can provide objectives for PRLP to try, saved in #user_objectives
  void setUserObjectives(const std::vector<std::vector<double> >& obj);

  /// Return #user_tight_points
  inline std::vector<std::vector<double> > getUserTightPoints() const { return this->user_tight_points; }
  /// @brief User can provide set of points for PRLP to try and find cuts that are tight on those points, saved in #user_objectives
  void setUserTightPoints(const std::vector<std::vector<double> >& tight_points);
  ///@} // get/set methods

  /// @brief Generate VPCs from a disjunction (e.g., arising from a partial branch-and-bound tree)
  virtual void generateCuts(const OsiSolverInterface&, OsiCuts&, const CglTreeInfo = CglTreeInfo());

  /// @brief Any time a cut is added, it should go through this method
  bool addCut(const OsiRowCut& cut, OsiCuts& cuts, const CutType& type, const ObjectiveType& cutHeur, 
      const OsiSolverInterface* const origSolver = NULL,
      const bool check_violation = false);

  ///@{
  /// @name Print methods

  /// @brief For each cut type in \link CglVPC::CutType::NUM_CUT_TYPES NUM_CUT_TYPES \endlink, print #CutTypeName, #numCutsOfType
  inline void printCutsOfType(FILE* logfile = stdout) const {
    for (int i = 0; i < static_cast<int>(CutType::NUM_CUT_TYPES); i++) {
      fprintf(logfile, "%s,", CutTypeName[i].c_str());
      fprintf(logfile, "%d,", numCutsOfType[i]);
      fprintf(logfile, "\n");
    }
    fflush(logfile);
  } /* printCutsOfType */

  /// @brief For each objective in \link CglVPC::ObjectiveType::NUM_OBJECTIVE_TYPES NUM_OBJECTIVE_TYPES \endlink, print #ObjectiveTypeName, #numObjFromHeur
  inline void printObjFromHeur(FILE* logfile = stdout) const {
      for (int i = 0; i < static_cast<int>(ObjectiveType::NUM_OBJECTIVE_TYPES); i++) {
        fprintf(logfile, "OBJ_%s,", ObjectiveTypeName[i].c_str());
        fprintf(logfile, "%d,", numObjFromHeur[i]);
        fprintf(logfile, "\n");
      }
      fflush(logfile);
    } /* printObjFromHeur */

  /// @brief For each objective in \link CglVPC::ObjectiveType::NUM_OBJECTIVE_TYPES NUM_OBJECTIVE_TYPES \endlink, print #ObjectiveTypeName, #numCutsFromHeur
  inline void printCutsFromHeur(FILE* logfile = stdout) const {
      for (int i = 0; i < static_cast<int>(ObjectiveType::NUM_OBJECTIVE_TYPES); i++) {
        fprintf(logfile, "CUTS_%s,", ObjectiveTypeName[i].c_str());
        fprintf(logfile, "%d,", numCutsFromHeur[i]);
        fprintf(logfile, "\n");
      }
      fflush(logfile);
    } /* printCutsFromHeur */

  /// @brief For each objective in \link CglVPC::ObjectiveType::NUM_OBJECTIVE_TYPES NUM_OBJECTIVE_TYPES \endlink, print #ObjectiveTypeName, #numFailsFromHeur
  inline void printFailsFromHeur(FILE* logfile = stdout) const {
    for (int i = 0; i < static_cast<int>(ObjectiveType::NUM_OBJECTIVE_TYPES); i++) {
      fprintf(logfile, "FAILS_%s,", ObjectiveTypeName[i].c_str());
      fprintf(logfile, "%d,", numFailsFromHeur[i]);
      fprintf(logfile, "\n");
    }
    fflush(logfile);
  } /* printFailsFromHeur */

  /// @brief For each objective in \link CglVPC::FailureType::NUM_FAILURE_TYPES NUM_FAILURE_TYPES \endlink, print #FailureTypeName, #numFails
  inline void printFailures(FILE* logfile = stdout) const {
    for (int i = 0; i < static_cast<int>(FailureType::NUM_FAILURE_TYPES); i++) {
      fprintf(logfile, "%s,", FailureTypeName[i].c_str());
      fprintf(logfile, "%d,", numFails[i]);
      fprintf(logfile, "\n");
    }
    fflush(logfile);
  } /* printFailures */
  ///@}

protected:
  std::vector<std::vector<double> > user_objectives; ///< User can provide objectives to try for PRLP, in original (structural) space
  std::vector<std::vector<double> > user_tight_points; ///< User can provide points that PRLP will attempt to find cuts with small distance to, where the points are in the original (structural) space
  Disjunction* disjunction = NULL; ///< Pointer to Disjunction used for this round of cuts
  PRLPData prlpData; ///< PRLPData for this round of cuts

  /// @brief Clear old information before another round of cuts
  void setupAsNew();

  /// @brief Initialize class members
  void initialize(const CglVPC* const source = NULL, const VPCParametersNamespace::VPCParameters* const param = NULL);

  /// @brief Calculate values from given instance held in \p solver, and store this in \p probData
  void getProblemData(OsiSolverInterface* const solver, ProblemData& probData,
      const ProblemData* const origProbData = NULL,
      const bool enable_factorization = true);

  /// @brief Do everything necessary to setup PRLP constraints
  CglVPC::ExitReason setupConstraints(OsiSolverInterface* const vpcsolver, OsiCuts& cuts);

  /// @brief Get basis cone for each disjunctive term
  void genDepth1PRCollection(const OsiSolverInterface* const vpcsolver,
      const OsiSolverInterface* const tmpSolver, const ProblemData& origProbData,
      const ProblemData& tmpProbData, const int term_ind);

  /// @brief Try objectives in nonbasic space or structural space, as desired
  CglVPC::ExitReason tryObjectives(OsiCuts& cuts,
      const OsiSolverInterface* const origSolver, const OsiCuts* const structSICs);

  /// @brief Return #num_cuts >= getCutLimit()
  inline bool reachedCutLimit() const {
    return (num_cuts >= getCutLimit()); 
  } /* reachedCutLimit */

  /// @brief Call reachedTimeLimit(const std::string&,const double) const
  inline bool reachedTimeLimit(const VPCTimeStats& timeName, const double max_time) const {
    return timer.reachedTimeLimit(VPCTimeStatsName[static_cast<int>(timeName)], max_time);
  } /* reachedTimeLimit */

  /// @brief Return \link TimeStats::get_total_time() timer.get_total_time(timeName) \endlink > \p max_time
  inline bool reachedTimeLimit(const std::string& timeName, const double max_time) const {
    //return (timer.get_total_time(timeName) > max_time);
    return timer.reachedTimeLimit(timeName, max_time);
  } /* reachedTimeLimit */

//  inline bool reachedTimelimit(const std::chrono::time_point<std::chrono::high_resolution_clock>& start_chrono) const {
//    std::chrono::duration<double> elapsed_seconds = std::chrono::high_resolution_clock::now() - start_chrono;
//    const bool reached_limit = (elapsed_seconds.count() > this->param.get(TIMELIMIT));
////    if (reached_limit)
////      exitReason = ExitReason::TIME_LIMIT_EXIT;
//    return reached_limit;
//  }

  /// @brief Check whether too many unsuccessful objective attempts have been made
  bool reachedFailureLimit(const int num_cuts, const int num_fails,
      const double few_cuts_fail_threshold = 0.95,
      const double many_cuts_fail_threshold = 0.90,
      const double many_obj_fail_threshold = 0.80,
      const double time_fail_threshold = 0.66) const;

  /// @brief Set #exitReason to \p exitReason, end all running clocks in #timer
  void finish(CglVPC::ExitReason exitReason = CglVPC::ExitReason::UNKNOWN);
}; /* Class CglVPC */

/** @brief Match status from generating disjunction to status of CglVPC exit */
CglVPC::ExitReason matchStatus(const DisjExitReason status);

