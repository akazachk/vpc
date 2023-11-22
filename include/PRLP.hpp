/**
 * @file PRLP.hpp
 * @brief PRLP class, adapting OsiSolverInterface to VPC setting
 *
 * Specific adaptation of an OsiSolverInterface to iteratively generate V-polyhedral cuts
 * by selecting objectives for a given disjunction and point-ray collection derived from that disjunction
 *
 * Though it is designed as an OsiClpSolverInterface, it can be adapted for other OsiSolverInterfaces
 *
 * @author A. M. Kazachkov
 * @date 2019-02-14
 */
#pragma once
#include "CglVPC.hpp"
#include "utility.hpp"
#include "SolverInterface.hpp"

/**
 * @brief Specific adaptation of an OsiSolverInterface to iteratively generate V-polyhedral cuts
 *
 * Selects objectives for a given disjunction and point-ray collection derived from that disjunction
 * Though it is designed as an OsiClpSolverInterface, it can be adapted (with some work) for other OsiSolverInterfaces
 */
class PRLP : public SolverInterface {
public:
  CglVPC* owner; ///< CglVPC instance; not owned by PRLP
  const Disjunction* disj; ///< disjunction corresponding to this PRLP; not owned by PRLP
  int disjID; ///< ID of disjunction corresponding to this PRLP
  const CglVPC::PRLPData* prlpData; ///< data for this PRLP; not owned by PRLP
  std::vector<int> nonZeroColIndex;
  std::vector<double> ortho; ///< orthogonality of each row
  int numPoints; ///< number of rows with nonzero constant side
  int numRays; ///< number of rows with zero constant side
  double density; ///< average number of nonzeros in the constraint matrix
  int num_obj_tried; ///< number of times PRLP was solved
  int num_cuts; ///< number of times a cut was generated
  int num_failures; ///< number of times the solve was unsuccessful (e.g., duplicate cut)

  /// @brief Default constructor
  PRLP();

  /// @brief Owner constructor
  PRLP(CglVPC* owner, const Disjunction* const disj, const int disjID, const CglVPC::PRLPData* const prlpData);

  /// @brief Copy constructor
  PRLP(const PRLP& source);

  /// @brief Destructor
  ~PRLP();

  /// @brief Assignment operator 
  PRLP& operator=(const PRLP& source);

  /// @brief Clone
  virtual OsiSolverInterface* clone(bool copyData = true) const;

  /// @brief Setup PRLP
  CglVPC::ExitReason setup(const double scale);

  /** 
   * @brief Loop for generating cuts
   * @return \# cuts generated 
   */
  int targetStrongAndDifferentCuts(const double beta, OsiCuts& cuts,
      const OsiSolverInterface* const origSolver,
      const OsiCuts* const structSICs, const std::string& timeName);

  /**
   * @brief Solve with a structural-space objective
   * @return \# cuts generated
   */
  int solve(OsiCuts& cuts, const std::vector<double>& obj,
      const OsiSolverInterface* const origSolver, const double beta,
      const OsiCuts* const structSICs);

protected:
  /// @brief Keep a row index and the activity on that row (to sort by)
  struct rowAndActivity {
    int row;
    double activity;
  };

  /// @brief Define sorting of rowAndActivity by activity, descending
  struct by_activity_desc {
    bool operator()(rowAndActivity const &a, rowAndActivity const &b) {
      return a.activity > b.activity;
    }
  };

  /// @brief Define sorting of rowAndActivity by activity, ascending
  struct by_activity_asc {
    bool operator()(rowAndActivity const &a, rowAndActivity const &b) {
      return a.activity < b.activity;
    }
  };

  /// @brief Define a way to sort rays by orthogonality
  struct compareRayInfo {
    int index; ///< index of the ray
    double min_ortho; ///< minimum angle with any other ray
    double avg_ortho; ///< average ange with other rays

    /// @brief Sort by #min_ortho, or if equal, by #avg_ortho
    bool operator<(const compareRayInfo& rhs) const {
      return (this->min_ortho > rhs.min_ortho
          || (isVal(this->min_ortho, rhs.min_ortho, 0.1)
              && this->avg_ortho > rhs.avg_ortho));
    }
  };

  /// @brief Setup PRLP members
  void initialize(const PRLP* const source = NULL);

  /// @brief Convert structural point to nonbasic space and install as objective
  void setObjectiveFromStructuralPoint(const double* const pointVals,
      const double* const pointSlack, const OsiSolverInterface* const origSolver,
      const bool inNBSpace);

  /// @brief Try PRLP on current objective
  int tryOneObjective(std::vector<int>& numTimesTightRow,
      std::vector<int>& numTimesTightColLB,
      std::vector<int>& numTimesTightColUB, OsiCuts& cuts,
      const OsiSolverInterface* const origSolver, const double beta,
      const OsiCuts* const structSICs, const bool inNBSpace,
      const CglVPC::ObjectiveType cutHeur, const bool tryExtraHard = false);

  /// @brief Resolve and then reset limits
  int resolvePRLP(const bool tryExtraHard = false);

  /// @brief Generate one cut
  int genCut(OsiCuts& cuts,
      const OsiSolverInterface* const origSolver, const double beta,
      const OsiCuts* const structSICs, const bool inNBSpace,
      const CglVPC::ObjectiveType cutHeur, const bool tryExtraHard = false);

  /// @brief Helper function to check a cut
  int genCutHelper(OsiCuts & cuts, const OsiSolverInterface* const origSolver,
      const double beta, const OsiCuts* const structSICs, const bool inNBSpace,
      const CglVPC::ObjectiveType cutHeur);

  /// @brief Process whether a cut was generated, and if not, why not
  int exitGenCut(const int num_cuts_generated, const CglVPC::ObjectiveType cutHeur);

  /// @brief Sort points and rays
  void setupForTargetedCutGeneration(std::vector<rowAndActivity>& pointIndex,
      std::vector<rowAndActivity>& rayIndex);

  /// @brief Update points and rays in sorted order
  int updateStepForTargetedCutGeneration(std::vector<int>& numTimesTightRow,
      std::vector<int>& numTimesTightColLB, std::vector<int>& numTimesTightColUB) const;

  /// @brief Update sorting of rays
  void updateRaySetForTargetedCutGeneration(std::set<compareRayInfo>& sortedRays,
      const std::vector<rowAndActivity>& rayIndex, const int init_num_cuts, int& num_old_cuts,
      const OsiCuts& cuts);

  /// @brief Target cuts that have small slack wrt to a given point
  int findCutsTightOnPoint(std::vector<int>& numTimesTightRow,
      std::vector<int>& numTimesTightColLB,
      std::vector<int>& numTimesTightColUB, const int point_row_ind,
      const CglVPC::ObjectiveType& cutHeur,
      const double beta, OsiCuts& cuts,
      const OsiSolverInterface* const origSolver,
      const OsiCuts* const structSICs, const std::string& timeName,
      const bool inNBSpace, const int MAX_NUM_OBJ_PER_POINT);

  /// @brief Add Gomory cuts, then try to separate the resulting point
  int iterateDeepestCutPostGomory(OsiCuts & cuts,
      const OsiSolverInterface* const origSolver, const double beta,
      const OsiCuts* structSICs, const bool inNBSpace);

  /// @brief Track why we exit from iterateDeepestCutPostGomory()
  int exitFromIterateDeepestCutPostGomory(const int num_cuts_gen,
      OsiSolverInterface* const MSICSolver);
}; /* class PRLP */
