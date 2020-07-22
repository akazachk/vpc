// This code has been adpated (by Aleksandr M. Kazachkov) from OsiChooseVariable.cpp in the Osi project of COIN-OR
//
// The original license notice is below:
// Copyright (C) 2006, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include <string>
#include <cassert>
#include <cfloat>
#include <cmath>
#include "CoinPragma.hpp"
#include "OsiSolverInterface.hpp"
#include "OsiAuxInfo.hpp"
#include "OsiSolverBranch.hpp"
#include "CoinWarmStartBasis.hpp"
#include "CoinPackedMatrix.hpp"
#include "CoinTime.hpp"
#include "CoinSort.hpp"
#include "CoinFinite.hpp"
#include "OsiChooseStrongCustom.hpp"

/*  This is a utility function which does strong branching on
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
int OsiChooseStrongCustom::doStrongBranching(OsiSolverInterface * solver,
    OsiBranchingInformation *info, int numberToDo, int returnCriterion) {

  // Might be faster to extend branch() to return bounds changed
  double * saveLower = NULL;
  double * saveUpper = NULL;
  int numberColumns = solver->getNumCols();
  solver->markHotStart();
  const double * lower = info->lower_;
  const double * upper = info->upper_;
  saveLower = CoinCopyOfArray(info->lower_, numberColumns);
  saveUpper = CoinCopyOfArray(info->upper_, numberColumns);
  numResults_ = 0;
  int returnCode = 0;
  double timeStart = CoinCpuTime();
  for (int iDo = 0; iDo < numberToDo; iDo++) {
    OsiHotInfo * result = results_ + iDo;
    // For now just 2 way
    OsiBranchingObject * branch = result->branchingObject();
    assert(branch->numberBranches() == 2);
    /*
     Try the first direction.  Each subsequent call to branch() performs the
     specified branch and advances the branch object state to the next branch
     alternative.)
     */
    OsiSolverInterface * thisSolver = solver;
    if (branch->boundBranch()) {
      // ordinary
      branch->branch(solver);
      // maybe we should check bounds for stupidities here?
      solver->solveFromHotStart();
    } else {
      // adding cuts or something 
      thisSolver = solver->clone();
      branch->branch(thisSolver);
      // set hot start iterations
      int limit;
      thisSolver->getIntParam(OsiMaxNumIterationHotStart, limit);
      thisSolver->setIntParam(OsiMaxNumIteration, limit);
      thisSolver->resolve();
    }
    // can check if we got solution
    // status is 0 finished, 1 infeasible and 2 unfinished and 3 is solution
    int status0 = result->updateInformation(thisSolver, info, this);
    numberStrongIterations_ += thisSolver->getIterationCount();
    if (status0 == 3) {
      // new solution already saved
      if (trustStrongForSolution_) {
        info->cutoff_ = goodObjectiveValue_;
        status0 = 0;
      }
    }
    if (solver != thisSolver)
      delete thisSolver;
    // Restore bounds
    for (int j = 0; j < numberColumns; j++) {
      if (saveLower[j] != lower[j])
        solver->setColLower(j, saveLower[j]);
      if (saveUpper[j] != upper[j])
        solver->setColUpper(j, saveUpper[j]);
    }
    /*
     Try the next direction
     */
    thisSolver = solver;
    if (branch->boundBranch()) {
      // ordinary
      branch->branch(solver);
      // maybe we should check bounds for stupidities here?
      solver->solveFromHotStart();
    } else {
      // adding cuts or something 
      thisSolver = solver->clone();
      branch->branch(thisSolver);
      // set hot start iterations
      int limit;
      thisSolver->getIntParam(OsiMaxNumIterationHotStart, limit);
      thisSolver->setIntParam(OsiMaxNumIteration, limit);
      thisSolver->resolve();
    }
    // can check if we got solution
    // status is 0 finished, 1 infeasible and 2 unfinished and 3 is solution
    int status1 = result->updateInformation(thisSolver, info, this);
    numberStrongDone_++;
    numberStrongIterations_ += thisSolver->getIterationCount();
    if (status1 == 3) {
      // new solution already saved
      if (trustStrongForSolution_) {
        info->cutoff_ = goodObjectiveValue_;
        status1 = 0;
      }
    }
    if (solver != thisSolver)
      delete thisSolver;
    // Restore bounds
    for (int j = 0; j < numberColumns; j++) {
      if (saveLower[j] != lower[j])
        solver->setColLower(j, saveLower[j]);
      if (saveUpper[j] != upper[j])
        solver->setColUpper(j, saveUpper[j]);
    }
    /*
     End of evaluation for this candidate variable. Possibilities are:
     * Both sides below cutoff; this variable is a candidate for branching.
     * Both sides infeasible or above the objective cutoff: no further action
     here. Break from the evaluation loop and assume the node will be purged
     by the caller.
     * One side below cutoff: Install the branch (i.e., fix the variable). Possibly break
     from the evaluation loop and assume the node will be reoptimised by the
     caller.
     */
    numResults_++;
    if (status0 == 1 && status1 == 1) {
      // infeasible or both sides have objective value above cutoff
      returnCode = -1;
      break; // exit loop
    } else if (status0 == 1 || status1 == 1) {
      numberStrongFixed_++;
      if (!returnCriterion) {
        returnCode = 1;
      } else {
        returnCode = 2;
        break;
      }
    }
    bool hitMaxTime = (CoinCpuTime() - timeStart > info->timeRemaining_);
    if (hitMaxTime) {
      returnCode = 3;
      break;
    }
  }
  delete[] saveLower;
  delete[] saveUpper;
  // Delete the snapshot
  solver->unmarkHotStart();
  return returnCode;
} /* doStrongBranching */

//##############################################################################

OsiChooseStrongCustom::OsiChooseStrongCustom() :
    OsiChooseVariable(), shadowPriceMode_(0), pseudoCosts_(), results_(NULL), numResults_(
        0) {
  copyOurStuff(NULL);
}

OsiChooseStrongCustom::OsiChooseStrongCustom(const OsiSolverInterface * solver) :
    OsiChooseVariable(solver), shadowPriceMode_(0), pseudoCosts_(), results_(
        NULL), numResults_(0) {
  // create useful arrays
  pseudoCosts_.initialize(solver_->numberObjects());
  copyOurStuff(NULL);
}

OsiChooseStrongCustom::OsiChooseStrongCustom(const OsiChooseStrongCustom & rhs) :
    OsiChooseVariable(rhs), shadowPriceMode_(rhs.shadowPriceMode_), pseudoCosts_(
        rhs.pseudoCosts_), results_(NULL), numResults_(0) {
  copyOurStuff(&rhs);
}

OsiChooseStrongCustom &
OsiChooseStrongCustom::operator=(const OsiChooseStrongCustom & rhs) {
  if (this != &rhs) {
    OsiChooseVariable::operator=(rhs);
    copyOurStuff(&rhs);
    shadowPriceMode_ = rhs.shadowPriceMode_;
    pseudoCosts_ = rhs.pseudoCosts_;
    delete[] results_;
    results_ = NULL;
    numResults_ = 0;
  }
  return *this;
}

OsiChooseStrongCustom::~OsiChooseStrongCustom() {
  delete[] results_;
}

// Clone
OsiChooseVariable *
OsiChooseStrongCustom::clone() const {
  return new OsiChooseStrongCustom(*this);
}

// Save our information
void OsiChooseStrongCustom::copyOurStuff(
    const OsiChooseStrongCustom* const rhs) {
  if (rhs) {
    this->method_ = rhs->method_;
    this->secondBestObjectIndex_ = rhs->secondBestObjectIndex_;
    this->secondBestWhichWay_ = rhs->secondBestWhichWay_;
  } else {
    this->method_ = DEFAULT_METHOD_; // see description in the processObject method
    this->secondBestObjectIndex_ = -1;
    this->secondBestWhichWay_ = 0;
  }
} /* copyOurStuff */

// Initialize
int OsiChooseStrongCustom::setupList(OsiBranchingInformation *info,
    bool initialize) {
  if (initialize) {
    status_ = -2;
    delete[] goodSolution_;
    bestObjectIndex_ = -1;
    numberStrongDone_ = 0;
    numberStrongIterations_ = 0;
    numberStrongFixed_ = 0;
    goodSolution_ = NULL;
    goodObjectiveValue_ = COIN_DBL_MAX;
  }
  numberOnList_ = 0;
  numberUnsatisfied_ = 0;
  int numberObjects = solver_->numberObjects();
  if (numberObjects > pseudoCosts_.numberObjects()) {
    // redo useful arrays
    pseudoCosts_.initialize(numberObjects);
  }
  double check = -COIN_DBL_MAX;
  int checkIndex = 0;
  int bestPriority = COIN_INT_MAX;
  int maximumStrong = CoinMin(numberStrong_, numberObjects);
  int putOther = numberObjects;
  int i;
  for (i = 0; i < numberObjects; i++) {
    list_[i] = -1;
    useful_[i] = 0.0;
  }
  OsiObject ** object = info->solver_->objects();
  // Get average pseudo costs and see if pseudo shadow prices possible
  int shadowPossible = shadowPriceMode_;
  if (shadowPossible) {
    for (i = 0; i < numberObjects; i++) {
      if (!object[i]->canHandleShadowPrices()) {
        shadowPossible = 0;
        break;
      }
    }
    if (shadowPossible) {
      int numberRows = solver_->getNumRows();
      const double * pi = info->pi_;
      double sumPi = 0.0;
      for (i = 0; i < numberRows; i++)
        sumPi += fabs(pi[i]);
      sumPi /= static_cast<double>(numberRows);
      // and scale back
      sumPi *= 0.01;
      info->defaultDual_ = sumPi; // switch on
      int numberColumns = solver_->getNumCols();
      int size = CoinMax(numberColumns, 2 * numberRows);
      info->usefulRegion_ = new double[size];
      CoinZeroN(info->usefulRegion_, size);
      info->indexRegion_ = new int[size];
    }
  }
  double sumUp = 0.0;
  double numberUp = 0.0;
  double sumDown = 0.0;
  double numberDown = 0.0;
  const double* upTotalChange = pseudoCosts_.upTotalChange();
  const double* downTotalChange = pseudoCosts_.downTotalChange();
  const int* upNumber = pseudoCosts_.upNumber();
  const int* downNumber = pseudoCosts_.downNumber();
  const int numberBeforeTrusted = pseudoCosts_.numberBeforeTrusted();
  for (i = 0; i < numberObjects; i++) {
    sumUp += upTotalChange[i];
    numberUp += upNumber[i];
    sumDown += downTotalChange[i];
    numberDown += downNumber[i];
  }
  double upMultiplier = (1.0 + sumUp) / (1.0 + numberUp);
  double downMultiplier = (1.0 + sumDown) / (1.0 + numberDown);
  // Say feasible
  bool feasible = true;
#if 0
  int pri[]= {10,1000,10000};
  int priCount[]= {0,0,0};
#endif
  for (i = 0; i < numberObjects; i++) {
    int way;
    double value = object[i]->infeasibility(info, way);
    if (value > 0.0) {
      numberUnsatisfied_++;
      if (value == COIN_DBL_MAX) {
        // infeasible
        feasible = false;
        break;
      }
      int priorityLevel = object[i]->priority();
#if 0
      for (int k=0;k<3;k++) {
        if (priorityLevel==pri[k])
        priCount[k]++;
      }
#endif
      // Better priority? Flush choices.
      if (priorityLevel < bestPriority) {
        for (int j = maximumStrong - 1; j >= 0; j--) {
          if (list_[j] >= 0) {
            int iObject = list_[j];
            list_[j] = -1;
            useful_[j] = 0.0;
            list_[--putOther] = iObject;
          }
        }
        maximumStrong = CoinMin(maximumStrong, putOther);
        bestPriority = priorityLevel;
        check = -COIN_DBL_MAX;
        checkIndex = 0;
      }
      if (priorityLevel == bestPriority) {
        // Modify value
        sumUp = upTotalChange[i] + 1.0e-30;
        numberUp = upNumber[i];
        sumDown = downTotalChange[i] + 1.0e-30;
        numberDown = downNumber[i];
        double upEstimate = object[i]->upEstimate();
        double downEstimate = object[i]->downEstimate();
        if (shadowPossible < 2) {
          upEstimate =
              numberUp ?
                  ((upEstimate * sumUp) / numberUp) :
                  (upEstimate * upMultiplier);
          if (numberUp < numberBeforeTrusted)
            upEstimate *= (numberBeforeTrusted + 1.0) / (numberUp + 1.0);
          downEstimate =
              numberDown ?
                  ((downEstimate * sumDown) / numberDown) :
                  (downEstimate * downMultiplier);
          if (numberDown < numberBeforeTrusted)
            downEstimate *= (numberBeforeTrusted + 1.0) / (numberDown + 1.0);
        } else {
          // use shadow prices always
        }
        value = calculateValue(method_, object[i]->columnNumber(), downEstimate, upEstimate, 0);
        if (value > check) { // this is probably not what we want in this new code
          //add to list
          int iObject = list_[checkIndex];
          if (iObject >= 0) {
            assert(list_[putOther - 1] < 0);
            list_[--putOther] = iObject;  // to end
          }
          list_[checkIndex] = i;
          assert(checkIndex < putOther);
          useful_[checkIndex] = value;
          // find worst
          check = COIN_DBL_MAX;
          maximumStrong = CoinMin(maximumStrong, putOther);
          for (int j = 0; j < maximumStrong; j++) {
            if (list_[j] >= 0) {
              if (useful_[j] < check) {
                check = useful_[j];
                checkIndex = j;
              }
            } else {
              check = 0.0;
              checkIndex = j;
              break;
            }
          }
        } else {
          // to end
          assert(list_[putOther - 1] < 0);
          list_[--putOther] = i;
          maximumStrong = CoinMin(maximumStrong, putOther);
        }
      } else {
        // worse priority
        // to end
        assert(list_[putOther - 1] < 0);
        list_[--putOther] = i;
        maximumStrong = CoinMin(maximumStrong, putOther);
      }
    }
  }
#if 0
  printf("%d at %d, %d at %d and %d at %d\n",priCount[0],pri[0],
      priCount[1],pri[1],priCount[2],pri[2]);
#endif
  // Get list
  numberOnList_ = 0;
  if (feasible) {
    for (i = 0; i < CoinMin(maximumStrong, putOther); i++) {
      if (list_[i] >= 0) {
        list_[numberOnList_] = list_[i];
        useful_[numberOnList_++] = -useful_[i];
      }
    }
    if (numberOnList_) {
      // Sort 
      CoinSort_2(useful_, useful_ + numberOnList_, list_);
      // move others
      i = numberOnList_;
      for (; putOther < numberObjects; putOther++)
        list_[i++] = list_[putOther];
      assert(i == numberUnsatisfied_);
      if (!numberStrong_)
        numberOnList_ = 0;
    }
  } else {
    // not feasible
    numberUnsatisfied_ = -1;
  }
  // Get rid of any shadow prices info
  info->defaultDual_ = -1.0; // switch off
  delete[] info->usefulRegion_;
  delete[] info->indexRegion_;
  return numberUnsatisfied_;
}

void OsiChooseStrongCustom::resetResults(int num) {
  delete[] results_;
  numResults_ = 0;
  results_ = new OsiHotInfo[num];
}

/**
 * AMK changed here 2018-08-14
 * method:
 * 0/default: what is in OsiChooseStrong (0.85 * min change + 0.15 * max change)
 * 1: default, but with tie-breaking
 * 2: choose var that maximizes min change, then select by maximizes max change
 * 3: choose var that is second-best in terms of default, with tie-breaking
 * 4: choose var that is second-best in terms of maximizes min change, then select by maximizes max change
 * -x: -1 * (x+1)
 */
double OsiChooseStrongCustom::calculateValue(const int chosen_method,
    const int variable_index, const double downEstimate,
    const double upEstimate, const int whichCriterion) {
  if (chosen_method == 6) {
    return -1 * variable_index;
  }
  int method = chosen_method;
  // Process the method
  double mult = 1.;
  if (method < 0) {
    mult = -1.;
    method = -1 * (method + 1);
  }
  double maxmin_coeff = SAVED_MAXMIN_CRITERION_;
  switch (method) {
    case 2:
    case 4: {
      maxmin_coeff = 1.;
      break;
    }
    case 5: {
      maxmin_coeff = 0.;
      break;
    }
  }
  if (whichCriterion) {
    maxmin_coeff = 1. - maxmin_coeff;
  }
  return mult
      * (maxmin_coeff * CoinMin(upEstimate, downEstimate)
          + (1.0 - maxmin_coeff) * CoinMax(upEstimate, downEstimate));
} /* calculateValue */

int OsiChooseStrongCustom::processObject(const int chosen_method,
    OsiSolverInterface* solver, const int iObject, const double downEstimate,
    const double upEstimate, double& bestTrusted,
    double& bestTrustedSecondCriterion, double& secondBestTrusted,
    double& secondBestTrustedSecondCriterion, const int givenReturnCode) {
  int returnCode = givenReturnCode;
  const double firstValue = calculateValue(chosen_method, solver->object(iObject)->columnNumber(), downEstimate,
      upEstimate, 0);
  const double secondValue = calculateValue(chosen_method, solver->object(iObject)->columnNumber(), downEstimate,
      upEstimate, 1);
  if (firstValue > bestTrusted) {
    secondBestTrusted = bestTrusted;
    secondBestTrustedSecondCriterion = bestTrustedSecondCriterion;
    secondBestObjectIndex_ = bestObjectIndex_;
    secondBestWhichWay_ = bestWhichWay_;

    bestTrusted = firstValue;
    bestTrustedSecondCriterion = secondValue;
    bestObjectIndex_ = iObject;
    bestWhichWay_ = upEstimate > downEstimate ? 0 : 1;
    // but override if there is a preferred way
    const OsiObject * obj = solver->object(iObject);
    if (obj->preferredWay() >= 0 && obj->infeasibility())
      bestWhichWay_ = obj->preferredWay();
    if (returnCode)
      returnCode = 2;
  } else if (chosen_method != 0 && firstValue == bestTrusted
      && secondValue > bestTrustedSecondCriterion) {
    secondBestTrusted = bestTrusted;
    secondBestTrustedSecondCriterion = bestTrustedSecondCriterion;
    secondBestObjectIndex_ = bestObjectIndex_;
    secondBestWhichWay_ = bestWhichWay_;

    bestTrusted = firstValue;
    bestTrustedSecondCriterion = secondValue;
    bestObjectIndex_ = iObject;
    bestWhichWay_ = upEstimate > downEstimate ? 0 : 1;
    // but override if there is a preferred way
    const OsiObject * obj = solver->object(iObject);
    if (obj->preferredWay() >= 0 && obj->infeasibility())
      bestWhichWay_ = obj->preferredWay();
    if (returnCode)
      returnCode = 2;
  } else if (firstValue > secondBestTrusted) {
    secondBestTrusted = firstValue;
    secondBestTrustedSecondCriterion = secondValue;
    secondBestObjectIndex_ = iObject;
    secondBestWhichWay_ = upEstimate > downEstimate ? 0 : 1;
    // but override if there is a preferred way
    const OsiObject * obj = solver->object(iObject);
    if (obj->preferredWay() >= 0 && obj->infeasibility())
      secondBestWhichWay_ = obj->preferredWay();
    if (returnCode)
      returnCode = 2;
  } else if (chosen_method != 0 && firstValue == secondBestTrusted
      && secondValue > secondBestTrustedSecondCriterion) {
    secondBestTrusted = firstValue;
    secondBestTrustedSecondCriterion = secondValue;
    secondBestObjectIndex_ = iObject;
    secondBestWhichWay_ = upEstimate > downEstimate ? 0 : 1;
    // but override if there is a preferred way
    const OsiObject * obj = solver->object(iObject);
    if (obj->preferredWay() >= 0 && obj->infeasibility())
      secondBestWhichWay_ = obj->preferredWay();
    if (returnCode)
      returnCode = 2;
  }

  return returnCode;
} /* processObject */

/* Choose a variable
 Returns -
 -1 Node is infeasible, or pruned for other reasons (e.g., rlp2_presolved: both sides are integer-feasible for some variable)
 0  Normal termination - we have a candidate
 1  All looks satisfied - no candidate
 2  We can change the bound on a variable - but we also have a strong branching candidate
 3  We can change the bound on a variable - but we have a non-strong branching candidate
 4  We can change the bound on a variable - no other candidates
 We can pick up branch from whichObject() and whichWay()
 We can pick up a forced branch (can change bound) from whichForcedObject() and whichForcedWay()
 If we have a solution then we can pick up from goodObjectiveValue() and goodSolution()
 */
int OsiChooseStrongCustom::chooseVariable(OsiSolverInterface * solver,
    OsiBranchingInformation *info, bool fixVariables) {
  if (method_ == 7) {
    bestObjectIndex_ = -1;
    bestWhichWay_ = -1;

    // Check if it is integer-feasible or infeasible; otherwise try to find a branching variable
    OsiSolverInterface* resolver = solver->clone();
    resolver->initialSolve();
    int returnCode = (!resolver->isProvenOptimal()) ? -1 : 0;

    if (this->feasibleSolution(info, solver->getColSolution(),
        solver->numberObjects(),
        const_cast<const OsiObject **>(solver->objects()))) {
      saveSolution(solver);
      returnCode = 1;
    } /* found feasible solution */
    else if (returnCode == 0 && info->depth_ < solver->numberObjects()) {
      double* sol = const_cast<double*>(info->solution_);
      for (int obj_ind = info->depth_; obj_ind < solver->numberObjects(); obj_ind++) {
        // Change solution to be fractional temporarily, just in case it is not
        // This will not affect anything on the Osi side of things, because chooseOsiBranch replaces info's solution
        // It is important to do this so that when createBranch is called, we get the right down_ and up_ values
        const int var = solver->objects()[obj_ind]->columnNumber();
        double val = sol[var];
        const double lb = info->lower_[var]; //solver->getColLower()[var];
        const double ub = info->upper_[var]; //solver->getColUpper()[var];
        if (ub > lb) { // need to check this in case some were fixed previously during tightenBounds call in CbcModel
          if (val >= ub - 1e-3) {
            val -= (ub - lb) / 2.;
          } else if (val <= lb + 1e-3) {
            val += (ub - lb) / 2.;
          }
          sol[var] = val;
          info->integerTolerance_ = -1;
          bestObjectIndex_ = obj_ind;
          bestWhichWay_ = 0;
          solver->objects()[obj_ind]->setWhichWay(bestWhichWay_);
          break;
        }
      }
      assert(bestObjectIndex_ >= 0);
    } /* solver feasible but not integer-feasible */
    else { // should be infeasible
      assert(resolver->isProvenPrimalInfeasible());
    } /* solver not feasible */
    if (resolver) {
      delete resolver;
    }
    return returnCode;
  } /* method == 7 */
//  else if (method_ == 6) {
//    // Check if it is integer-feasible or infeasible; otherwise try to find a branching variable
//    OsiSolverInterface* resolver = solver->clone();
//    resolver->initialSolve();
//    int returnCode = (!resolver->isProvenOptimal()) ? -1 : 0;
//    if (returnCode == 0) {
//      bestObjectIndex_ = solver->getNumCols();
//      bestWhichWay_ = -1;
//      int numberLeft = CoinMin(numberStrong_ - numberStrongDone_,
//          numberUnsatisfied_);
//      for (int i = 0; i < numberLeft; i++) {
//        int iObject = list_[i];
//        if (iObject < bestObjectIndex_) {
//          bestObjectIndex_ = iObject;
//        }
//      }
//      if (bestObjectIndex_ == solver->getNumCols()) {
//        bestObjectIndex_ = -1;
//        returnCode = 1;
//      } else {
//        bestWhichWay_ = 0;
//        solver->objects()[bestObjectIndex_]->setWhichWay(bestWhichWay_);
//      }
//    } else {
//      assert(resolver->isProvenPrimalInfeasible());
//    }
//    if (resolver) {
//      delete resolver;
//    }
//    return returnCode;
//  } /* method == 6 */
  if (numberUnsatisfied_) {
    const double* upTotalChange = pseudoCosts_.upTotalChange();
    const double* downTotalChange = pseudoCosts_.downTotalChange();
    const int* upNumber = pseudoCosts_.upNumber();
    const int* downNumber = pseudoCosts_.downNumber();
    int numberBeforeTrusted = pseudoCosts_.numberBeforeTrusted();
    // Somehow we can get here with it 0 !
    if (!numberBeforeTrusted) {
      numberBeforeTrusted = 5;
      pseudoCosts_.setNumberBeforeTrusted(numberBeforeTrusted);
    }

    int numberLeft = CoinMin(numberStrong_ - numberStrongDone_,
        numberUnsatisfied_);
    int numberToDo = 0;
    resetResults(numberLeft);
    int returnCode = 0;
    bestObjectIndex_ = -1;
    bestWhichWay_ = -1;
    firstForcedObjectIndex_ = -1;
    firstForcedWhichWay_ = -1;
    double bestTrusted = -COIN_DBL_MAX;
    double bestTrustedSecondCriterion = -COIN_DBL_MAX;
    double secondBestTrusted = -COIN_DBL_MAX;
    double secondBestTrustedSecondCriterion = -COIN_DBL_MAX;
    for (int i = 0; i < numberLeft; i++) {
      int iObject = list_[i];
      if (upNumber[iObject] < numberBeforeTrusted
          || downNumber[iObject] < numberBeforeTrusted) {
        results_[numberToDo++] = OsiHotInfo(solver, info, solver->objects(),
            iObject);
      } else {
        const OsiObject * obj = solver->object(iObject);
        double upEstimate = (upTotalChange[iObject] * obj->upEstimate())
            / upNumber[iObject];
        double downEstimate = (downTotalChange[iObject] * obj->downEstimate())
            / downNumber[iObject];
        processObject(method_, solver, iObject, downEstimate, upEstimate,
            bestTrusted, bestTrustedSecondCriterion, secondBestTrusted,
            secondBestTrustedSecondCriterion, 0);
      }
    }
    int numberFixed = 0;
    if (numberToDo) {
      returnCode = doStrongBranching(solver, info, numberToDo, 1);
      if (returnCode >= 0 && returnCode <= 2) {
        if (returnCode) {
          returnCode = 4;
          if (bestObjectIndex_ >= 0)
            returnCode = 3;
        }
        for (int i = 0; i < numResults_; i++) {
          int iObject = results_[i].whichObject();
          double upEstimate;
          if (results_[i].upStatus() != 1) {
            assert(results_[i].upStatus() >= 0);
            upEstimate = results_[i].upChange();
          } else {
            // infeasible - just say expensive
            if (info->cutoff_ < 1.0e50)
              upEstimate = 2.0 * (info->cutoff_ - info->objectiveValue_);
            else
              upEstimate = 2.0 * fabs(info->objectiveValue_);
            if (firstForcedObjectIndex_ < 0) {
              firstForcedObjectIndex_ = iObject;
              firstForcedWhichWay_ = 0;
            }
            numberFixed++;
            if (fixVariables) {
              const OsiObject * obj = solver->object(iObject);
              OsiBranchingObject * branch = obj->createBranch(solver, info, 0);
              branch->branch(solver);
              delete branch;
            }
          }
          double downEstimate;
          if (results_[i].downStatus() != 1) {
            assert(results_[i].downStatus() >= 0);
            downEstimate = results_[i].downChange();
          } else {
            // infeasible - just say expensive
            if (info->cutoff_ < 1.0e50)
              downEstimate = 2.0 * (info->cutoff_ - info->objectiveValue_);
            else
              downEstimate = 2.0 * fabs(info->objectiveValue_);
            if (firstForcedObjectIndex_ < 0) {
              firstForcedObjectIndex_ = iObject;
              firstForcedWhichWay_ = 1;
            }
            numberFixed++;
            if (fixVariables) {
              const OsiObject * obj = solver->object(iObject);
              OsiBranchingObject * branch = obj->createBranch(solver, info, 1);
              branch->branch(solver);
              delete branch;
            }
          }
          // AMK changed here 2018-08-14
          returnCode = processObject(method_, solver, iObject, downEstimate,
              upEstimate, bestTrusted, bestTrustedSecondCriterion,
              secondBestTrusted, secondBestTrustedSecondCriterion, returnCode);
        }
      } else if (returnCode == 3) {
        // max time - just choose one
        bestObjectIndex_ = list_[0];
        bestWhichWay_ = 0;
        if (numberOnList_ > 1) {
          secondBestObjectIndex_ = list_[1];
          secondBestWhichWay_ = 0;
        }
        returnCode = 0;
      }
    } else {
      bestObjectIndex_ = list_[0];
      bestWhichWay_ = 0;
      if (numberOnList_ > 1) {
        secondBestObjectIndex_ = list_[1];
        secondBestWhichWay_ = 0;
      }
    }
    const int method = (method_ < 0) ? -1 * (method_ + 1) : method_;
    const bool useSecondBest = (method == 3 || method == 4);
    if (useSecondBest && secondBestObjectIndex_ >= 0) {
      bestObjectIndex_ = secondBestObjectIndex_;
      bestWhichWay_ = secondBestWhichWay_;
    }
    if (method_ == 6) {
      bestWhichWay_ = 0;
    }
    if (bestObjectIndex_ >= 0) {
      OsiObject * obj = solver->objects()[bestObjectIndex_];
      obj->setWhichWay(bestWhichWay_);
    }
    if (numberFixed == numberUnsatisfied_ && numberFixed)
      returnCode = 4;
    return returnCode;
  } else {
    return 1;
  }
} /* chooseVariable */
