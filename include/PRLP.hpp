// Name:     PRLP.hpp
// Author:   A. M. Kazachkov
// Date:     2019-Feb-14
//-----------------------------------------------------------------------------
#pragma once
#include "CglVPC.hpp"
#include "utility.hpp"
#include "VPCParameters.hpp"

/**********************************************************************************************************
 * PRLP Class
 * Specific adaptation of an OsiSolverInterface to iteratively generate V-polyhedral cuts
 * by selecting objectives for a given disjunction and point-ray collection derived from that disjunction
 *
 * Though it is designed as an OsiClpSolverInterface, it can be adapted for other OsiSolverInterfaces
 **********************************************************************************************************/

class PRLP : public SolverInterface {
public:
  CglVPC* owner; // not owned by PRLP
  std::vector<int> nonZeroColIndex;
  std::vector<double> ortho;
  int numPoints, numRays;
  double density;
  int num_obj_tried, num_cuts, num_failures;

  /** Default constructor */
  PRLP();

  /** Owner constructor */
  PRLP(CglVPC* owner);

  /** Copy constructor */
  PRLP(const PRLP& source);

  /** Destructor */
  ~PRLP();

  /** Assignment operator */
  PRLP& operator=(const PRLP& source);

  /** Clone */
  virtual OsiSolverInterface* clone(bool copyData = true) const;

  /** Setup PRLP */
  bool setup(const double scale);

  /** Loop for generating cuts; return # cuts generated */
  int targetStrongAndDifferentCuts(const double beta, OsiCuts& cuts,
      const OsiSolverInterface* const origSolver,
      const OsiCuts* const structSICs, const std::string& timeName);

protected:
  struct rowAndActivity {
    int row;
    double activity;
  };

  struct by_activity_desc {
    bool operator()(rowAndActivity const &a, rowAndActivity const &b) {
      return a.activity > b.activity;
    }
  };

  struct by_activity_asc {
    bool operator()(rowAndActivity const &a, rowAndActivity const &b) {
      return a.activity < b.activity;
    }
  };

  struct compareRayInfo {
    int index;
    double min_ortho;
    double avg_ortho;

    bool operator<(const compareRayInfo& rhs) const {
      return (this->min_ortho > rhs.min_ortho
          || (isVal(this->min_ortho, rhs.min_ortho, 0.1)
              && this->avg_ortho > rhs.avg_ortho));
    }
  };

  void initialize(const PRLP* const source = NULL);

  void setObjectiveFromStructuralPoint(const double* const pointVals,
      const double* const pointSlack, const OsiSolverInterface* const origSolver,
      const bool inNBSpace);
  int tryOneObjective(std::vector<int>& numTimesTightRow,
      std::vector<int>& numTimesTightColLB,
      std::vector<int>& numTimesTightColUB, OsiCuts& cuts,
      const OsiSolverInterface* const origSolver, const double beta,
      const OsiCuts* const structSICs, const bool inNBSpace,
      const CglVPC::ObjectiveType cutHeur, const bool tryExtraHard = false);

  int resolvePRLP(const bool tryExtraHard = false);
  int genCut(OsiCuts& cuts,
      const OsiSolverInterface* const origSolver, const double beta,
      const OsiCuts* const structSICs, const bool inNBSpace,
      const CglVPC::ObjectiveType cutHeur, const bool tryExtraHard = false);
  int genCutHelper(OsiCuts & cuts, const OsiSolverInterface* const origSolver,
      const double beta, const OsiCuts* const structSICs, const bool inNBSpace,
      const CglVPC::ObjectiveType cutHeur);
  int exitGenCut(const int num_cuts_generated, const CglVPC::ObjectiveType cutHeur);

  void setupForTargetedCutGeneration(std::vector<rowAndActivity>& pointIndex,
      std::vector<rowAndActivity>& rayIndex);
  int updateStepForTargetedCutGeneration(std::vector<int>& numTimesTightRow,
      std::vector<int>& numTimesTightColLB, std::vector<int>& numTimesTightColUB) const;
  void updateRaySetForTargetedCutGeneration(std::set<compareRayInfo>& sortedRays,
      const std::vector<rowAndActivity>& rayIndex, const int init_num_cuts, int& num_old_cuts,
      const OsiCuts& cuts);
  int findCutsTightOnPoint(std::vector<int>& numTimesTightRow,
      std::vector<int>& numTimesTightColLB,
      std::vector<int>& numTimesTightColUB, const int point_row_ind,
      const CglVPC::ObjectiveType& cutHeur,
      const double beta, OsiCuts& cuts,
      const OsiSolverInterface* const origSolver,
      const OsiCuts* const structSICs, const std::string& timeName,
      const bool inNBSpace, const int MAX_NUM_OBJ_PER_POINT);

  int iterateDeepestCutPostGomory(OsiCuts & cuts,
      const OsiSolverInterface* const origSolver, const double beta,
      const OsiCuts* structSICs, const bool inNBSpace);
  int exitFromIterateDeepestCutPostGomory(const int num_cuts_gen,
      OsiSolverInterface* const MSICSolver);
}; /* class PRLP */
