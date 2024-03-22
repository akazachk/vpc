/** 
 * @file nbspace.cpp
 * @author A. M. Kazachkov
 * @date 2019-Feb-20
 */
#include "nbspace.hpp"

#include <cmath> // abs
#include <CoinPackedVector.hpp>
#include <stdexcept> // Include the necessary header for exceptions

#include "CutHelper.hpp" // setOsiRowCut
#include "SolverHelper.hpp"
#include "VPCParameters.hpp"
using namespace VPCParametersNamespace;
#include "utility.hpp" // dotProduct

/**
 * Assumed to be on the same set of variables
 * We also need that the right-hand side of the original solver remains unchanged (except the bounds)
 */
void setCompNBCoorPoint(CoinPackedVector& vec, double& objViolation,
    const VPCParameters& params, const OsiSolverInterface* const tmpSolver,
    const OsiSolverInterface* const origSolver,
    const std::vector<int>& nonBasicVarIndex,
    const std::vector<double>& nonBasicReducedCost, const int deletedVar) {
  // If a variable was deleted to make tmpSolver, then we will need to adjust the index
  // It is -1 because we need the index in the *new* space
  const int deletedVarAdjustment = (deletedVar >= 0) ? -1 : 0;

  const int numCols = origSolver->getNumCols();
  const int numNB = numCols; // again that nb fixed vars are treated as nb

  // For keeping the packed coordinates to input as a new point
  std::vector<int> packedIndex;
  std::vector<double> packedElem;

  packedIndex.reserve(numNB);
  packedElem.reserve(numNB);

  objViolation = tmpSolver->getObjValue() - origSolver->getObjValue();

  // Sometimes points have tiny values
  // These may be "necessary" for validity
  // E.g., if the point is (p1,p2), and c is the nb obj, we need c1 * p1 + c2 * p2 >= 1 (actually 1 - 1e-7)
  // It may be that p1 is tiny, but c1 is very large, and c2 * p2 < 0 by that amount, roughly
  // Other times, if the non-tiny values satisfy the objective cut on their own, let's ignore the tiny ones
  // We will not add the tiny indices + values until the end, when we check if they are necessary
  std::vector<int> tinyIndex;
  std::vector<double> tinyElem;
  // double tinyObjOffset = 0.; // how much c . p changes when tiny values are ignored
  double nonTinyObj = 0.;

  // All coefficients are going to be nonnegative, since we work in the complemented NB space
  for (int i = 0; i < numNB; i++) {
    const int currVar = nonBasicVarIndex[i];
    const int currVarTmpIndex = (currVar < deletedVar) ? currVar : currVar + deletedVarAdjustment;
    double tmpVal = 0., origVal = 0.;
    if (currVar < numCols) {
      tmpVal = tmpSolver->getColSolution()[currVarTmpIndex];
      origVal = origSolver->getColSolution()[currVar];
    } else {
      const int currRow = currVar - numCols;
      tmpVal = tmpSolver->getRightHandSide()[currRow] - tmpSolver->getRowActivity()[currRow];
      origVal = 0.0;
    }
    // ***
    // For *struct* var:
    // If the variable is lower-bounded, then we need to shift it up for the comp space
    // That is, xComp = x - lb
    // If the variable is upper-bounded, we need to get the value in the complemented space
    // That is, xComp = ub - x
    // What to do if the variable is fixed NB or free NB?
    // In the former case, thinking of it as a lower-bounded variable, it will simply be zero
    // in any solution, including in the original one at v, in which case we don't need to
    // worry about it for the packed ray.
    // NB: The latter case we could have handled by preprocessing the free variables out
    // ****
    // For *slack* var:
    // In the original NB space at v, all slacks were at zero.
    // We simply look at b[row] - activity[row]
    // For >= rows, slack is non-positive
    // For <= rows, slack is non-negative
    // We take absolute value to deal with this
    const double newVal = std::abs(tmpVal - origVal);

    if (!isZero(newVal, params.get(EPS))) {
      packedIndex.push_back(i);
      packedElem.push_back(newVal);
      nonTinyObj += newVal * nonBasicReducedCost[i];
    } // non-tiny
    else if (!isZero(newVal * nonBasicReducedCost[i], params.get(EPS))) {
      tinyIndex.push_back(i);
      tinyElem.push_back(newVal);
      // tinyObjOffset += newVal * nonBasicReducedCost[i];
    } // tiny
  } // end iterating over nonbasic elements from original basis

  // Check whether tiny elements are needed
  bool useTinyElements = //lessThanVal(nonTinyObj, 0., param.getEPS())
          !isVal(nonTinyObj, objViolation, params.get(EPS));
  const int numNonTiny = packedIndex.size();
  if (useTinyElements) {
    // We may not need all of the tiny elements
    // Sort the tiny elements in decreasing order and keep adding them until the criterion is reached
    const int numTinyTotal = tinyIndex.size();
    CoinPackedVector tinyVec(numTinyTotal, tinyIndex.data(), tinyElem.data());
    tinyVec.sortDecrElement();
    int numTiny = 0;
    for (int i = 0; i < numTinyTotal && useTinyElements; i++) {
      const double val = tinyVec.getElements()[i];
      numTiny++;
      nonTinyObj += val * nonBasicReducedCost[tinyVec.getIndices()[i]];
      useTinyElements = //lessThanVal(nonTinyObj, 0., param.getEPS())
          !isVal(nonTinyObj, objViolation, params.get(EPS));
    }

    std::vector<int> newPackedIndex(numNonTiny + numTiny, 0);
    std::vector<double> newPackedElem(numNonTiny + numTiny, 0);
    int ind = 0, ind1 = 0, ind2 = 0;
    const double scale = 1.; // perhaps in a future version we will scale things
    while (ind1 < numNonTiny) {
      while (ind2 < numTiny && tinyIndex[ind2] < packedIndex[ind1]) {
        newPackedIndex[ind] = tinyVec.getIndices()[ind2];
        newPackedElem[ind] = tinyVec.getElements()[ind2] * scale;
        ind++;
        ind2++;
      }
      newPackedIndex[ind] = packedIndex[ind1];
      newPackedElem[ind] = packedElem[ind1] * scale;
      ind++;
      ind1++;
    }
    while (ind2 < numTiny) {
      newPackedIndex[ind] = tinyVec.getIndices()[ind2];
      newPackedElem[ind] = tinyVec.getElements()[ind2] * scale;
      ind++;
      ind2++;
    }

    vec.setVector((int) newPackedIndex.size(),
        newPackedIndex.data(), newPackedElem.data(), false);
  } else {
    vec.setVector((int) packedIndex.size(), packedIndex.data(),
        packedElem.data(), false);
  }

  if (!isVal(nonTinyObj, objViolation, params.get(EPS))) {
    warning_msg(warnstring,
                "Point: Calculated dot product with obj differs from solver's. Obj viol from solver: %.8f. Calculated: %.8f.\n",
                objViolation, nonTinyObj);
    if (!isVal(nonTinyObj, objViolation, params.get(doubleConst::DIFFEPS))) {
      throw std::runtime_error("Point: Calculated dot product with obj differs from solver's");
    }
  }
} /* setCompNBCoorPoint */

/**
 * Note that we may have deleted a variable to get to tmpSolver
 * In that case, tmpNBVar is in the space of tmpSolver
 */
void setCompNBCoorRay(CoinPackedVector& vec, const double* ray, double& objViolation, double& scale,
    const VPCParameters& params, const OsiSolverInterface* const tmpSolver,
    const OsiSolverInterface* const origSolver,
    const std::vector<int>& rowOfOrigNBVar,
    const std::vector<int>& nonBasicVarIndex,
    const std::vector<double>& nonBasicReducedCost,
    const int tmpNBVar, const int deletedVar,
    const bool rayNeedsCalculation) {
  const int numCols = origSolver->getNumCols();
  const int numNB = numCols; // TODO again assuming that nb fixed vars are treated as ub vars

  std::vector<int> NBRayIndex;
  std::vector<double> NBRayElem;

  NBRayIndex.reserve(numNB);
  NBRayElem.reserve(numNB);

  // Adjustments for the deletedVar (if there is one)
  const int deletedVarAdjustment = (deletedVar >= 0) ? 1 : 0;
  const int tmpNBVarOrigIndex = (tmpNBVar < deletedVar) ? tmpNBVar : tmpNBVar + deletedVarAdjustment;

  // If the variable corresponding to this ray is upper-bounded, its ray direction is negative
  const int tmpRayMult = (isNonBasicUBVar(tmpSolver, tmpNBVar)) ? -1 : 1;

  // Get objective dot product
  if (tmpNBVar < tmpSolver->getNumCols()) {
    objViolation = (double) tmpRayMult * tmpSolver->getReducedCost()[tmpNBVar];
  } else {
    const int tempIndex = tmpNBVar - tmpSolver->getNumCols();
    objViolation = (double) -1 * tmpRayMult * tmpSolver->getRowPrice()[tempIndex];
  }

  // Perhaps we will scale the rays to get better numerics
  int minExponent = 1000, maxExponent = 1;
  int minExponentTiny = 0;

  // Sometimes rays have tiny values
  // These may be "necessary" for validity
  // E.g., if the ray is (r1,r2), and c is the nb obj, we need c1 * r1 + c2 * r2 >= 0 (actually -1e-7)
  // It may be that r1 is tiny, but c1 is very large, and c2 * r2 < 0 by that amount, roughly
  // Other times, if the non-tiny values satisfy the objective cut on their own, let's ignore the tiny ones
  // We will not add the tiny indices + values until the end, when we check if they are necessary
  std::vector<int> tinyIndex;
  std::vector<double> tinyElem;
  // double tinyObjOffset = 0.; // how much c . r changes when tiny values are ignored
  double nonTinyObj = 0.;

  // We loop through the *originally* non-basic variables (at v) to get their ray components
  // Some of these may now be basic...
  for (int c = 0; c < numNB; c++) {
    const int origNBVar = nonBasicVarIndex[c];
    const int origNBVarRow = rowOfOrigNBVar[c];
    const int rayMult = isNonBasicUBVar(origSolver, origNBVar) ? -1 : 1;
    const int newRayMult = tmpRayMult * rayMult;

    // If the current ray, from NBVar, is the same variable as origNBVar,
    // then we should insert at this point the ray direction and index.
    // Though we are working in the complemented non-basic space, the ray direction may
    // actually be -1.0 here, if the NBVar has switched to being tight at its other bound
    if (origNBVar == tmpNBVarOrigIndex) {
      NBRayIndex.push_back(c);
      NBRayElem.push_back(newRayMult);
      nonTinyObj += newRayMult * nonBasicReducedCost[c];
      if (0 < minExponent) {
        minExponent = 0;
      }
    } else if (origNBVarRow >= 0) {
      // If we are here, then the NB var is basic currently.
      // Thus, we need to get its component and store it.
      // For the variables that were upper-bounded, we actually go in the negative direction.
      const int flip = rayNeedsCalculation ? -1 : tmpRayMult; // Negate the use of tmpRayMult previously, if ray already calculated
//      const int rowMult =
//          (rayNeedsCalculation && (origNBVar >= origSolver->getNumCols())
//              && (origSolver->getRowSense()[origNBVarRow] == 'G')) ? -1 : 1;
      const double rayVal = flip * newRayMult * ray[origNBVarRow];
      if (!isZero(rayVal, params.get(doubleConst::RAYEPS))) {
        NBRayIndex.push_back(c);
        NBRayElem.push_back(rayVal);
        nonTinyObj += rayVal * nonBasicReducedCost[c];

        // Get exponent
        const double exp = std::log10(std::abs(rayVal));
        const int intExponent = static_cast<int>(((exp >= 0) ? 0.046 : -0.046) + exp);
        if (intExponent > maxExponent) {
          maxExponent = intExponent;
        }
        if (intExponent < minExponent) {
          minExponent = intExponent;
        }
      } // non-tiny
      else if (!isZero(rayVal * nonBasicReducedCost[c], params.get(doubleConst::RAYEPS))) {
        tinyIndex.push_back(c);
        tinyElem.push_back(rayVal);
        // tinyObjOffset += rayVal * nonBasicReducedCost[c];
      } // tiny
      /*
      else {
        if (rayVal != 0)
        printf("\nDEBUG: rejecting ray coeff on nb_var %d (variable %d) with value %.10e\n", c, origNBVar, rayVal);
      }
      */
    }
  } /* end iterating over non-basic elements from original basis */

  // Check whether tiny elements are needed
  bool useTinyElements = lessThanVal(nonTinyObj, 0., params.get(EPS))
          || !isVal(nonTinyObj, objViolation, params.get(EPS));
  const int numNonTiny = NBRayIndex.size();
  if (useTinyElements) {
    // We may not need all of the tiny elements
    // Sort the tiny elements in decreasing order and keep adding them until the criterion is reached
    const int numTinyTotal = tinyIndex.size();
    CoinPackedVector tinyVec(numTinyTotal, tinyIndex.data(), tinyElem.data());
    tinyVec.sortDecrElement();
    int numTiny = 0;
    for (int i = 0; i < numTinyTotal && useTinyElements; i++) {
      const double val = tinyVec.getElements()[i];
      numTiny++;
      nonTinyObj += val * nonBasicReducedCost[tinyVec.getIndices()[i]];
      useTinyElements = lessThanVal(nonTinyObj, 0., params.get(EPS))
          || !isVal(nonTinyObj, objViolation, params.get(EPS));

      // Get exponent
      const double exp = std::log10(std::abs(val));
      const int intExponent =
          static_cast<int>(((exp >= 0) ? 0.046 : -0.046) + exp);
      if (intExponent < minExponentTiny) {
        minExponentTiny = intExponent;
      }
    }

    int scale_exp = 0;
    const double avgExponentsTiny = 0.5 * (minExponentTiny + maxExponent);
    if (avgExponentsTiny < -0.25) {
      scale_exp = 1 - static_cast<int>(avgExponentsTiny + 0.25);
    }
    scale = std::pow(10., static_cast<double>(scale_exp));

    std::vector<int> newNBRayIndex(numNonTiny + numTiny, 0);
    std::vector<double> newNBRayElem(numNonTiny + numTiny, 0);
    int ind = 0, ind1 = 0, ind2 = 0;
    while (ind1 < numNonTiny) {
      while (ind2 < numTiny && tinyIndex[ind2] < NBRayIndex[ind1]) {
        newNBRayIndex[ind] = tinyVec.getIndices()[ind2];
        newNBRayElem[ind] = tinyVec.getElements()[ind2] * scale;
        ind++;
        ind2++;
      }
      newNBRayIndex[ind] = NBRayIndex[ind1];
      newNBRayElem[ind] = NBRayElem[ind1] * scale;
      ind++;
      ind1++;
    }
    while (ind2 < numTiny) {
      newNBRayIndex[ind] = tinyVec.getIndices()[ind2];
      newNBRayElem[ind] = tinyVec.getElements()[ind2] * scale;
      ind++;
      ind2++;
    }

    vec.setVector((int) newNBRayIndex.size(),
        newNBRayIndex.data(), newNBRayElem.data(), false);
  } else {
    // Get exponent for scaling
    // We will only scale if it is up, not down, by the heuristic that small values are bad
    int scale_exp = 0;
    const double avgExponents = 0.5 * (minExponent + maxExponent);
    if (avgExponents < -0.25) {
      scale_exp = 1 - static_cast<int>(avgExponents + 0.25);
    } else if (minExponent > 0) {
      scale_exp = -1 * minExponent;
    }
    scale = std::pow(10., static_cast<double>(scale_exp));

    if (!isVal(scale, 1.)) {
      for (int i = 0; i < numNonTiny; i++) {
        NBRayElem[i] *= scale;
      }
    }

    vec.setVector(numNonTiny, NBRayIndex.data(), NBRayElem.data(),
        false);
  }

  if (lessThanVal(nonTinyObj, 0., params.get(EPS))) {
    if (lessThanVal(nonTinyObj, 0., params.get(doubleConst::DIFFEPS))) {
      error_msg(errorstring,
          "Ray %d: dot product with obj < 0. Obj viol from solver: %.8f. Calculated: %.8f.\n",
          tmpNBVar, objViolation, nonTinyObj);
      writeErrorToLog(errorstring, params.logfile);
      exit(1);
    } else {
      warning_msg(warnstring,
          "Ray %d: dot product with obj < 0. Obj viol from solver: %.8f. Calculated: %.8f.\n",
          tmpNBVar, objViolation, nonTinyObj);
    }
  }

  if (!isVal(nonTinyObj, objViolation, params.get(EPS))) {
    if (!isVal(nonTinyObj, objViolation, params.get(doubleConst::DIFFEPS))) {
      error_msg(errorstring,
          "Ray %d: Calculated dot product with obj differs from solver's. Obj viol from solver: %.8f. Calculated: %.8f.\n",
          tmpNBVar, objViolation, nonTinyObj);
      writeErrorToLog(errorstring, params.logfile);
      exit(1);
    } else {
      warning_msg(warnstring,
          "Ray %d: Calculated dot product with obj differs from solver's. Obj viol from solver: %.8f. Calculated: %.8f.\n",
          tmpNBVar, objViolation, nonTinyObj);
    }
  }
} /* setCompNBCoor (Ray) */

/**
 * @brief Puts currColValue (necessary to provide) into complemented nonbasic space form
 *
 * Assumed to be for the same set of variables as in originalSolver
 * We also need that the right-hand side of the original solver remains unchanged (except the bounds)
 */
void setCompNBCoor(
    /// [out] Resulting vector is saved here
    CoinPackedVector& vec, 
    /// [out] Calculated objective violation
    double& objViolation,
    /// [in]  To get epsilon, diff-eps, and the logfile in case of an error
    const VPCParameters& params, 
    /// [in]  Point that is being converted (value of structural variables)
    const double* const currColValue,
    /// [in]  Slack values of point being converted
    const double* const currSlackValue,
    /// [in]  The solver that is used to determine the nonbasic space
    const OsiSolverInterface* const origSolver,
    /// [in]  The nonbasic space to which we wish to convert
    const std::vector<int>& nonBasicVarIndex,
    /// [in]  This is the objective vector in the nonbasic space
    const std::vector<double>& nonBasicReducedCost, 
    /// [in]  If nonnegative, then this is the index of a variable that was removed from origSolver
    const int deletedVar) {
  // If a variable was deleted to make tmpSolver, then we will need to adjust the index
  // It is -1 because we need the index in the *new* space
  const int deletedVarAdjustment = (deletedVar >= 0) ? -1 : 0;

  const int numCols = origSolver->getNumCols();
  const int numNB = numCols; // assuming that nb fixed vars are treated as nb

  // For keeping the packed coordinates to input as a new point
  std::vector<int> packedIndex;
  std::vector<double> packedElem;

  packedIndex.reserve(numNB);
  packedElem.reserve(numNB);

  objViolation = dotProduct(currColValue, origSolver->getObjCoefficients(), numCols) - dotProduct(origSolver->getColSolution(), origSolver->getObjCoefficients(), numCols);

  // Sometimes points have tiny values
  // These may be "necessary" for validity
  // E.g., if the point is (p1,p2), and c is the nb obj, we need c1 * p1 + c2 * p2 >= 1 (actually 1 - 1e-7)
  // It may be that p1 is tiny, but c1 is very large, and c2 * p2 < 0 by that amount, roughly
  // Other times, if the non-tiny values satisfy the objective cut on their own, let's ignore the tiny ones
  // We will not add the tiny indices + values until the end, when we check if they are necessary
  double nonTinyObj = 0.;
  std::vector<int> tinyIndex;
  std::vector<double> tinyElem;
  // double tinyObjOffset = 0.; // how much c . p changes when tiny values stop being ignored
#ifdef TRACE
  std::vector<int> superTinyIndex;
  std::vector<double> superTinyElem;
  // double superTinyObjOffset = 0.;
#endif

  // All coefficients are going to be nonnegative, since we work in the complemented NB space
  for (int i = 0; i < numNB; i++) {
    const int currVar = nonBasicVarIndex[i];
    const int currVarTmpIndex = (currVar < deletedVar) ? currVar : currVar + deletedVarAdjustment;
    double tmpVal = 0., origVal = 0.;
    if (currVar < numCols) {
      tmpVal = currColValue[currVarTmpIndex];
      origVal = origSolver->getColSolution()[currVar];
    } else {
      const int currRow = currVar - numCols;
      if (currSlackValue) {
        tmpVal = currSlackValue[currRow];
      } else {
        // Need to calculate the slack value if it is not provided to us
        const CoinShallowPackedVector vec = origSolver->getMatrixByRow()->getVector(currRow);
        tmpVal = origSolver->getRightHandSide()[currRow] - dotProduct(vec, currColValue);
      }
      origVal = 0.0;
    }
    // ***
    // For *struct* var:
    // If the variable is lower-bounded, then we need to shift it up for the comp space
    // That is, xComp = x - lb
    // If the variable is upper-bounded, we need to get the value in the complemented space
    // That is, xComp = ub - x
    // What to do if the variable is fixed NB or free NB?
    // In the former case, thinking of it as a lower-bounded variable, it will simply be zero
    // in any solution, including in the original one at v, in which case we don't need to
    // worry about it for the packed ray.
    // NB: The latter case we could have handled by preprocessing the free variables out
    // ****
    // For *slack* var:
    // In the original NB space at v, all slacks were at zero.
    // We simply look at b[row] - activity[row]
    // For >= rows, slack is non-positive
    // For <= rows, slack is non-negative
    // We take absolute value to deal with this
    const double newVal = std::abs(tmpVal - origVal);

    if (!isZero(newVal, params.get(EPS))) {
      packedIndex.push_back(i);
      packedElem.push_back(newVal);
      nonTinyObj += newVal * nonBasicReducedCost[i];
    } // non-tiny
    else if (!isZero(newVal * nonBasicReducedCost[i], params.get(EPS))) {
      tinyIndex.push_back(i);
      tinyElem.push_back(newVal);
      // tinyObjOffset += newVal * nonBasicReducedCost[i];
    } // tiny
#ifdef TRACE
    else if (!isZero(newVal * nonBasicReducedCost[i], 0.0)) {
      superTinyIndex.push_back(i);
      superTinyElem.push_back(newVal);
      // superTinyObjOffset += newVal * nonBasicReducedCost[i];
    } // super-tiny
#endif
  } // end iterating over nonbasic elements from original basis

  // Check whether tiny elements are needed
  bool useTinyElements = //lessThanVal(nonTinyObj, 0., param.getEPS())
      !isVal(nonTinyObj, objViolation, params.get(EPS));
  const int numNonTiny = packedIndex.size();
  if (useTinyElements) {
    // We may not need all of the tiny elements
    // Sort the tiny elements in decreasing order and keep adding them until the criterion is reached
    const int numTinyTotal = tinyIndex.size();
    CoinPackedVector tinyVec(numTinyTotal, tinyIndex.data(), tinyElem.data());
    tinyVec.sortDecrElement();
    int numTiny = 0;
    for (int i = 0; i < numTinyTotal && useTinyElements; i++) {
      const double val = tinyVec.getElements()[i];
      numTiny++;
      nonTinyObj += val * nonBasicReducedCost[tinyVec.getIndices()[i]];
      useTinyElements = //lessThanVal(nonTinyObj, 0., param.getEPS())
          !isVal(nonTinyObj, objViolation, params.get(EPS));
    }

    std::vector<int> newPackedIndex(numNonTiny + numTiny, 0);
    std::vector<double> newPackedElem(numNonTiny + numTiny, 0);
    int ind = 0, ind1 = 0, ind2 = 0;
    const double scale = 1.; // perhaps in a future version we will scale things
    while (ind1 < numNonTiny) {
      while (ind2 < numTiny && tinyIndex[ind2] < packedIndex[ind1]) {
        newPackedIndex[ind] = tinyVec.getIndices()[ind2];
        newPackedElem[ind] = tinyVec.getElements()[ind2] * scale;
        ind++;
        ind2++;
      }
      newPackedIndex[ind] = packedIndex[ind1];
      newPackedElem[ind] = packedElem[ind1] * scale;
      ind++;
      ind1++;
    }
    while (ind2 < numTiny) {
      newPackedIndex[ind] = tinyVec.getIndices()[ind2];
      newPackedElem[ind] = tinyVec.getElements()[ind2] * scale;
      ind++;
      ind2++;
    }

    vec.setVector((int) newPackedIndex.size(), newPackedIndex.data(),
        newPackedElem.data(), false);
  } else {
    vec.setVector((int) packedIndex.size(), packedIndex.data(),
        packedElem.data(), false);
  }

  if (!isVal(nonTinyObj, objViolation, params.get(EPS))) {
    if (!isVal(nonTinyObj, objViolation, params.get(doubleConst::DIFFEPS))) {
      error_msg(errorstring,
          "Calculated dot product with obj differs from solver's. Obj viol from solver: %.8f. Calculated: %.8f.\n",
          objViolation, nonTinyObj);
      writeErrorToLog(errorstring, params.logfile);
      exit(1);
    } else {
      warning_msg(warnstring,
          "Calculated dot product with obj differs from solver's. Obj viol from solver: %.8f. Calculated: %.8f.\n",
          objViolation, nonTinyObj);
    }
  }
} /* setCompNBCoor (generic) */

//void setStructCoor(const std::vector<int>& nonZeroColIndex, const double* const coeff,
//    const double beta, const OsiSolverInterface* const solver,
//    const std::vector<int> &nonBasicVarIndex, const bool inCompSpace,
//    const double EPS) {
//  // Useful info
//  const int num_cols = solver->getNumCols();
//  const int num_sparse_cols =
//      (nonZeroColIndex.size() > 0) ? nonZeroColIndex.size() : num_cols;
//
//  // Allocate space
//  double StructCutRHS = beta;
//  const double* rowRHS = solver->getRightHandSide();
//  const char* rowSense = solver->getRowSense();
//  const CoinPackedMatrix* mat = solver->getMatrixByRow();
//
//  std::vector<double> struct_coeff(num_cols);
//  for (int tmp_ind = 0; tmp_ind < num_sparse_cols; tmp_ind++) {
//    const int curr_ind =
//        (nonZeroColIndex.size() > 0) ? nonZeroColIndex[tmp_ind] : tmp_ind;
//    const int curr_var = nonBasicVarIndex[curr_ind];
//    const double curr_val = coeff[tmp_ind];
//
//    // If we are in the complemented space,
//    // at this point we need to do the complementing
//    // ub: y_j = u_j - x_j
//    // lb: y_j = x_j - l_j
//    if (curr_var < num_cols) {
//      struct_coeff[curr_var] = curr_val;
//      if (inCompSpace) {
//        if (isNonBasicUBCol(solver, curr_var)) {
//          StructCutRHS -= curr_val * solver->getColUpper()[curr_var];
//          struct_coeff[curr_var] *= -1;
//        } else if (isNonBasicLBCol(solver, curr_var)) {
//          StructCutRHS += curr_val * solver->getColLower()[curr_var];
//        }
//        //        if (solver->getModelPtr()->getColumnStatus(curr_var)
//        //            == ClpSimplex::atUpperBound) {
//        //          StructCutRHS -= curr_val * solver->getColUpper()[curr_var];
//        //          struct_coeff[curr_var] *= -1;
//        //        } else if (solver->getModelPtr()->getColumnStatus(curr_var)
//        //            == ClpSimplex::atLowerBound
//        //            || solver->getModelPtr()->getColumnStatus(curr_var)
//        //                == ClpSimplex::isFixed) {
//        //          StructCutRHS += curr_val * solver->getColLower()[curr_var];
//        //        }
//      }
//    } else {
//      const int curr_row = curr_var - num_cols;
//      const double mult =
//      //          (inCompSpace && isNonBasicLBSlack(solver, curr_row)) ? 1.0 : -1.0;
//          (inCompSpace && rowSense[curr_row] == 'G') ? 1.0 : -1.0; // TODO check
//
//      const CoinShallowPackedVector slackRow = mat->getVector(curr_row);
//      int slackRowNumElements = slackRow.getNumElements();
//      const int* slackRowIndices = slackRow.getIndices();
//      const double* slackRowCoeffs = slackRow.getElements();
//
//      StructCutRHS += mult * rowRHS[curr_row] * curr_val;
//      for (int j = 0; j < slackRowNumElements; j++) {
//        struct_coeff[slackRowIndices[j]] += mult * slackRowCoeffs[j] * curr_val;
//      }
//    }
//  } // iterate over columns
//
//} /* setStructCoor */

/***********************************************************************/
/**
 * @brief Set current cut by converting it from nonbasic coeffs given in coeff
 */
void setCutFromNBCoefficients(OsiRowCut* const cut,
    const std::vector<int>& nonZeroColIndex, const double* const coeff,
    const double beta, const OsiSolverInterface* const solver,
    const std::vector<int> &nonBasicVarIndex, const bool inCompSpace,
    const double EPS) {
  // Useful info
  const int num_cols = solver->getNumCols();
  const int num_sparse_cols =
      (nonZeroColIndex.size() > 0) ? nonZeroColIndex.size() : num_cols;

  // Allocate space
  double StructCutRHS = beta;
  const double* rowRHS = solver->getRightHandSide();
  const char* rowSense = solver->getRowSense();
  const CoinPackedMatrix* mat = solver->getMatrixByRow();

  std::vector<double> struct_coeff(num_cols);
  for (int tmp_ind = 0; tmp_ind < num_sparse_cols; tmp_ind++) {
    const int curr_ind = (nonZeroColIndex.size() > 0) ? nonZeroColIndex[tmp_ind] : tmp_ind;
    const int curr_var = nonBasicVarIndex[curr_ind];
    const double curr_val = coeff[tmp_ind];

    // If we are in the complemented space,
    // at this point we need to do the complementing
    // ub: y_j = u_j - x_j
    // lb: y_j = x_j - l_j
    if (curr_var < num_cols) {
      struct_coeff[curr_var] = curr_val;
      if (inCompSpace) {
        if (isNonBasicUBCol(solver, curr_var)) {
          StructCutRHS -= curr_val * solver->getColUpper()[curr_var];
          struct_coeff[curr_var] *= -1;
        } else if (isNonBasicLBCol(solver, curr_var)) {
          StructCutRHS += curr_val * solver->getColLower()[curr_var];
        }
//        if (solver->getModelPtr()->getColumnStatus(curr_var)
//            == ClpSimplex::atUpperBound) {
//          StructCutRHS -= curr_val * solver->getColUpper()[curr_var];
//          struct_coeff[curr_var] *= -1;
//        } else if (solver->getModelPtr()->getColumnStatus(curr_var)
//            == ClpSimplex::atLowerBound
//            || solver->getModelPtr()->getColumnStatus(curr_var)
//                == ClpSimplex::isFixed) {
//          StructCutRHS += curr_val * solver->getColLower()[curr_var];
//        }
      }
    } else {
      const int curr_row = curr_var - num_cols;
      const double mult =
//          (inCompSpace && isNonBasicLBSlack(solver, curr_row)) ? 1.0 : -1.0;
          (inCompSpace && rowSense[curr_row] == 'G') ? 1.0 : -1.0; // TODO check

      const CoinShallowPackedVector slackRow = mat->getVector(curr_row);
      int slackRowNumElements = slackRow.getNumElements();
      const int* slackRowIndices = slackRow.getIndices();
      const double* slackRowCoeffs = slackRow.getElements();

      StructCutRHS += mult * rowRHS[curr_row] * curr_val;
      for (int j = 0; j < slackRowNumElements; j++) {
        struct_coeff[slackRowIndices[j]] += mult * slackRowCoeffs[j] * curr_val;
      }
    }
  } // iterate over columns

  // Normalize
  /*
  const double absRHS = std::abs(StructCutRHS);
  if (!isZero(absRHS)) {
    for (int i = 0; i < numcols; i++) {
      struct_coeff[i] = struct_coeff[i] / absRHS;
    }

    if (greaterThanVal(StructCutRHS, 0.0)) {
      StructCutRHS = 1.0;
    } else if (lessThanVal(StructCutRHS, 0.0)) {
      StructCutRHS = -1.0;
    }
  }
  */
  setOsiRowCut(cut, std::vector<int>(), struct_coeff.size(), struct_coeff.data(), StructCutRHS, EPS);
} /* setCutFromJspaceCoefficients */
