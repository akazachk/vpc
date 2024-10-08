/**
 * @file OsiChooseStrongCustom.hpp
 *
 * @brief This code has been adapted (by Aleksandr M. Kazachkov) from OsiChooseVariable.hpp in the Osi project of COIN-OR.
 * 
 * The key details are located in the description of \link OsiChooseStrongCustom::method_ \endlink.
 *
 * @author A. M. Kazachkov
 * @date 2019-01-22
 */

// Old copyright info below:
// Copyright (C) 2006, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#pragma once

#include <string>
#include <vector>

#include "CoinWarmStartBasis.hpp"
#include "OsiBranchingObject.hpp"
#include "OsiChooseVariable.hpp"

/** @brief This class chooses a variable to branch on 
    
    By default, it will do strong branching, but also can be set to use two candidates, or pick in the reverse order.
    This can be useful if we wish to test a tree for cut generation that is different from a tree that would be used for branching.

    This chooses the variable and direction with reversereliability strong branching.

    The flow is :
    a) initialize the process.  This decides on strong branching list
       and stores indices of all infeasible objects  
    b) do strong branching on list.  If list is empty then just
       choose one candidate and return without strong branching.  If not empty then
       go through list and return best.  However we may find that the node is infeasible
       or that we can fix a variable.  If so we return and it is up to user to call
       again (after fixing a variable).
*/

class OsiChooseStrongCustom  : public OsiChooseVariable {
 
public:
    
  /// Default Constructor 
  OsiChooseStrongCustom ();

  /// Constructor from solver (so we can set up arrays etc)
  OsiChooseStrongCustom (const OsiSolverInterface * solver);

  /// Copy constructor 
  OsiChooseStrongCustom (const OsiChooseStrongCustom &);
   
  /// Assignment operator 
  OsiChooseStrongCustom & operator= (const OsiChooseStrongCustom& rhs);

  /// Clone
  virtual OsiChooseVariable * clone() const;

  virtual void copyOurStuff(const OsiChooseStrongCustom* const rhs);

  /// Destructor 
  virtual ~OsiChooseStrongCustom ();

  /** Sets up strong list and clears all if initialize is true.
      Returns number of infeasibilities. 
      If returns -1 then has worked out node is infeasible!
  */
  virtual int setupList ( OsiBranchingInformation *info, bool initialize);

  virtual double calculateValue(const int chosen_method,
    const int variable_index, const double downEstimate,
    const double upEstimate, const int whichCriterion);
  virtual int processObject(const int chosen_method, OsiSolverInterface* solver, 
    const int iObject, const double downEstimate, const double upEstimate,
    double& bestTrusted, double& bestTrustedSecondCriterion, 
    double& secondBestTrusted, double& secondBestTrustedSecondCriterion,
    const int givenReturnCode);

  /** Choose a variable
      Returns - 
     -1 Node is infeasible
     0  Normal termination - we have a candidate
     1  All looks satisfied - no candidate
     2  We can change the bound on a variable - but we also have a strong branching candidate
     3  We can change the bound on a variable - but we have a non-strong branching candidate
     4  We can change the bound on a variable - no other candidates
     We can pick up branch from bestObjectIndex() and bestWhichWay()
     We can pick up a forced branch (can change bound) from firstForcedObjectIndex() and firstForcedWhichWay()
     If we have a solution then we can pick up from goodObjectiveValue() and goodSolution()
     If fixVariables is true then 2,3,4 are all really same as problem changed
  */
  virtual int chooseVariable( OsiSolverInterface * solver, OsiBranchingInformation *info, bool fixVariables);

  /** Pseudo Shadow Price mode
      0 - off
      1 - use if no strong info
      2 - use if strong not trusted
      3 - use even if trusted
  */
  inline int shadowPriceMode() const
  { return shadowPriceMode_;}
  /// Set Shadow price mode
  inline void setShadowPriceMode(int value)
  { shadowPriceMode_ = value;}

  /** Accessor method to pseudo cost object*/
  const OsiPseudoCosts& pseudoCosts() const
  { return pseudoCosts_; }

  /** Accessor method to pseudo cost object*/
  OsiPseudoCosts& pseudoCosts()
  { return pseudoCosts_; }

  /** A few pass-through methods to access members of pseudoCosts_ as if they
      were members of OsiChooseStrongCustom object */
  inline int numberBeforeTrusted() const {
    return pseudoCosts_.numberBeforeTrusted(); }
  inline void setNumberBeforeTrusted(int value) {
    pseudoCosts_.setNumberBeforeTrusted(value); }
  inline int numberObjects() const {
    return pseudoCosts_.numberObjects(); }

  /// @brief return #method_
  inline int getMethod() const {
    return method_;
  }
  /// @brief set #method_
  inline void setMethod(int method) {
    method_ = method;
  }

protected:
  /**  This is a utility function which does strong branching on
       a list of objects and stores the results in OsiHotInfo.objects.
       On entry the object sequence is stored in the OsiHotInfo object
       and maybe more.
       It returns -
       -1 - one branch was infeasible both ways
       0 - all inspected - nothing can be fixed
       1 - all inspected - some can be fixed (returnCriterion==0)
       2 - may be returning early - one can be fixed (last one done) (returnCriterion==1) 
       3 - returning because max time
       
  */
  int doStrongBranching( OsiSolverInterface * solver, 
       OsiBranchingInformation *info,
       int numberToDo, int returnCriterion);

  /** Clear out the results array */
  void resetResults(int num);

  /// @brief Method for choosing the variable
  ///
  /// When #method_ != 0, calculate second value (flip weight on maxmin_coeff) for each object to break ties in processObject()
  /// when firstValue = bestTrusted and secondValue > bestTrustedSecondCriterion.
  /// When #method_ is 3 or 4, the first criterion is ignored, ando only the second is used
  ///
  /// \li <0: when value is x<0, proceed by reverse order of -1*(x+1)
  /// \li 0/default: what is in OsiChooseStrong (0.85 * min change + 0.15 * max change)
  /// \li 1: default, but with tie-breaking
  /// \li 2: choose var that maximizes min change, then select by maximizes max change; maxmin_coeff = 1
  /// \li 3: choose var that is second-best in terms of default, with tie-breaking
  /// \li 4: choose var that is second-best in terms of maximizes min change, then select by maximizes max change; maxmin_coeff = 1
  /// \li 5: maxmin_coeff = 0
  /// \li 6: reverse order, #bestWhichWay_ always 0
  /// \li 7: choose first not-fixed integer variable
  int method_;
  int secondBestObjectIndex_;
  int secondBestWhichWay_;
  /** Pseudo Shadow Price mode
      0 - off
      1 - use and multiply by strong info
      2 - use 
  */
  int shadowPriceMode_;

  /** The pseudo costs for the chooser */
  OsiPseudoCosts pseudoCosts_;

  /** The results of the strong branching done on the candidates where the
      pseudocosts were not sufficient */
  OsiHotInfo* results_;
  /** The number of OsiHotInfo objetcs that contain information */
  int numResults_;

  static constexpr int DEFAULT_METHOD_ = 0; ///< the default method is what is in the original OsiChooseStrong class
  static constexpr double SAVED_MAXMIN_CRITERION_ = 0.85; ///< factor for weighing max and min obj of children
};
