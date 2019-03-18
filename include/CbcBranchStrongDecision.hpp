// Original license is:
// $Id: CbcBranchStrongDecision.hpp Last modified: 2017-03-21 by Aleksandr M. Kazachkov $
// The original code is licensed under the terms of the Eclipse Public License (EPL).

// Edwin 11/10/2009-- carved out of CbcBranchActual
//
// Aleksandr M. Kazachkov 03/21/2017 Set method to 1

#ifndef CbcBranchStrongDecision_H
#define CbcBranchStrongDecision_H

#include "CbcBranchDecision.hpp"
/** Branching decision default class (modified by Aleksandr M. Kazachkov)

  This class implements a simple default algorithm
  (betterBranch()) for choosing a branching variable.
*/

class CbcBranchStrongDecision : public CbcBranchDecision {
public:
    // Default Constructor
    CbcBranchStrongDecision ();

    // Copy constructor
    CbcBranchStrongDecision ( const CbcBranchStrongDecision &);

    virtual ~CbcBranchStrongDecision();

    /// Clone
    virtual CbcBranchDecision * clone() const;

    /// Initialize, <i>e.g.</i> before the start of branch selection at a node
    virtual void initialize(CbcModel * model);

    /* Says whether this method can handle both branching methods - 1 betterBranch, 2 bestBranch, 3 both. */
    inline virtual int whichMethod() {
      return 3;
    }

    OsiChooseVariable* chooseMethod() const;

    /** \brief Compare two branching objects. Return nonzero if \p thisOne is
           better than \p bestSoFar.

      The routine compares branches using the values supplied in \p numInfUp and
      \p numInfDn until a solution is found by search, after which it uses the
      values supplied in \p changeUp and \p changeDn. The best branching object
      seen so far and the associated parameter values are remembered in the
      \c CbcBranchStrongDecision object. The nonzero return value is +1 if the
      up branch is preferred, -1 if the down branch is preferred.

      As the names imply, the assumption is that the values supplied for
      \p numInfUp and \p numInfDn will be the number of infeasibilities reported
      by the branching object, and \p changeUp and \p changeDn will be the
      estimated change in objective. Other measures can be used if desired.

      Because an \c CbcBranchStrongDecision object remembers the current best
      branching candidate (#bestObject_) as well as the values used in the
      comparison, the parameter \p bestSoFar is redundant, hence unused.
    */
    virtual int betterBranch(CbcBranchingObject * thisOne,
                             CbcBranchingObject * bestSoFar,
                             double changeUp, int numInfUp,
                             double changeDn, int numInfDn);
    /** Sets or gets best criterion so far */
    virtual void setBestCriterion(double value);
    virtual double getBestCriterion() const;

    /** \brief Compare N branching objects. Return index of best
        and sets way of branching in chosen object.

      This routine is used only after strong branching.
    */

    virtual int
    bestBranch (CbcBranchingObject ** objects, int numberObjects, int numberUnsatisfied,
                double * changeUp, int * numberInfeasibilitiesUp,
                double * changeDown, int * numberInfeasibilitiesDown,
                double objectiveValue) ;
private:

    /// Illegal Assignment operator
    CbcBranchStrongDecision & operator=(const CbcBranchStrongDecision& rhs);

    /// data
    /// OsiChooseStrong
//    OsiChooseStrong* chooseMethod_;

    /// "best" so far
    double bestCriterion_;

    /// Change up for best
    double bestChangeUp_;

    /// Number of infeasibilities for up
    int bestNumberUp_;

    /// Change down for best
    double bestChangeDown_;

    /// Pointer to best branching object
    CbcBranchingObject * bestObject_;

    /// Number of infeasibilities for down
    int bestNumberDown_;

};

#endif
