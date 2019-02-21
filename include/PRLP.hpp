// Name:     PRLP.hpp
// Author:   A. M. Kazachkov
// Date:     2019-Feb-14
//-----------------------------------------------------------------------------
#pragma once
#include "CglVPC.hpp"
#include "utility.hpp"
#include <OsiClpSolverInterface.hpp>

/**
 * PRLP Class
 * Specific adaptation of OsiClpSolverInterface to iteratively generate V-polyhedral cuts
 * by selecting objectives for a given disjunction and point-ray collection derived from that disjunction
 *
 * This class can be changed to work with other OsiSolverInterfaces
 */
class PRLP : public OsiClpSolverInterface {
public:

  CglVPC* owner;
  std::vector<int> nonZeroColIndex;
  std::vector<double> ortho;
  int numPoints, numRays;
  double density;

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
  bool setup(const double scale, const bool fixFirstPoint = false);

  /** Loop for generating cuts */
  int targetStrongAndDifferentCuts(const double beta, OsiCuts& cuts,
      int& num_cuts_total, int& num_obj_tried,
      const OsiSolverInterface* const origSolver,
      const OsiCuts* const structSICs, const std::string& timeName,
      const bool inNBSpace);

protected:
  void initialize(const PRLP* const source = NULL);
  bool setupHelper();
  int genCut(OsiCuts& cuts,
      const OsiSolverInterface* const origSolver, const double beta,
      const OsiCuts* const structSICs, int& num_cuts_total, const bool inNBSpace,
      const CglVPC::CutHeuristics cutHeur, const bool tryExtraHard);
  int genCutHelper(OsiCuts & cuts, const OsiSolverInterface* const origSolver,
      const double beta, const OsiCuts* const structSICs, int& num_cuts_total,
      const bool inNBSpace, const CglVPC::CutHeuristics cutHeur);
  int exitGenCut(const int num_cuts_generated);

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

  void setupForTargetedCutGeneration(std::vector<rowAndActivity>& pointIndex,
      std::vector<rowAndActivity>& rayIndex);
  int updateStepForTargetedCutGeneration(std::vector<int>& numTimesTightRow,
      std::vector<int>& numTimesTightColLB, std::vector<int>& numTimesTightColUB);
  void updateRaySetForTargetedCutGeneration(std::set<compareRayInfo>& sortedRays,
      const std::vector<rowAndActivity>& rayIndex, const int init_num_cuts, int& num_old_cuts,
      const OsiCuts& cuts);
  int tryOneObjective(int& num_failures, std::vector<int>& numTimesTightRow,
      std::vector<int>& numTimesTightColLB, std::vector<int>& numTimesTightColUB,
      const CglVPC::CutHeuristics& cutHeur, const double beta, OsiCuts& cuts,
      int& num_cuts_total, int& num_obj_tried,
      const OsiSolverInterface* const origSolver, const OsiCuts* const structSICs,
      const bool inNBSpace, const bool tryExtraHard = false);
  int findCutsTightOnPoint(int& num_failures,
      std::vector<int>& numTimesTightRow, std::vector<int>& numTimesTightColLB,
      std::vector<int>& numTimesTightColUB, const int point_row_ind,
      const int init_num_cuts, const CglVPC::CutHeuristics& cutHeur, const double beta,
      OsiCuts& cuts, int& num_cuts_total, int& num_obj_tried,
      const OsiSolverInterface* const origSolver, const OsiCuts* const structSICs,
      const std::string& timeName, const bool inNBSpace,
      const int MAX_NUM_OBJ_PER_POINT);
}; /* class PRLP */
