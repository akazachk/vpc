#include "utility.hpp"
#include "SolverHelper.hpp"

#include <fstream>
#include <sstream>
#include <sys/stat.h> // for fexists

/** We assume it is comma separated */
double getObjValueFromFile(const VPCParameters& params) {
  const char* opt_filename = params.get(stringParam::OPTFILE).c_str();

  // Take the full filename of the instance and remove any directory information
  std::string fullfilename = params.get(stringParam::FILENAME);
  size_t slashindex = fullfilename.find_last_of("/\\");
  std::string this_inst_name = (slashindex != std::string::npos) ? fullfilename.substr(slashindex+1) : fullfilename;

  // Put string after last '.' into string in_file_ext (accounting for compressed files)
  size_t found_dot = this_inst_name.find_last_of(".");
  std::string in_file_ext(this_inst_name.substr(found_dot + 1));
  this_inst_name = (found_dot != std::string::npos) ? this_inst_name.substr(0, found_dot) : this_inst_name;
  if (in_file_ext.compare("gz") == 0 || in_file_ext.compare("bz2") == 0) {
    found_dot = this_inst_name.find_last_of(".");

    if (found_dot == std::string::npos) {
      error_msg(errorstring,
          "Other than gz or bz2, cannot find the file extension (no '.' in input file name: %s).\n",
          fullfilename.c_str());
      writeErrorToLog(errorstring, params.logfile);
      exit(1);
    }

    in_file_ext = this_inst_name.substr(found_dot + 1);
    this_inst_name = this_inst_name.substr(0, found_dot);
  }
  if (in_file_ext.compare("mps") != 0 && in_file_ext.compare("lp") != 0) {
    error_msg(errorstring, "Could not identify instance name from file %s.\n", fullfilename.c_str());
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }

  if (!opt_filename) {
    return std::numeric_limits<double>::lowest();
  }

  std::ifstream infile(opt_filename);
  if (infile.is_open()) {
    std::string line;
    while (std::getline(infile, line)) {
      std::istringstream iss(line);
      if (line.empty()) {
        continue;
      }
      std::string inst_name;
      if (!(std::getline(iss, inst_name, ','))) {
        warning_msg(warnstring,
            "Could not read instance name. String is %s.\n",
            line.c_str());
        continue;
      }
      if (inst_name == this_inst_name) {
        try {
          std::string token;
          if (!(std::getline(iss, token, ','))) {
            throw;
          }
          return std::stod(token);
        } catch (std::exception& e) {
          warning_msg(warnstring,
              "Could not read optimal value. String is %s.\n",
              line.c_str());
          continue;
        }
      }
    }
    infile.close();
  } else {
    // If we were not able to open the file, throw an error
    error_msg(errorstring, "Not able to open obj file.\n");
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }

  return std::numeric_limits<double>::lowest();
} /* getObjValueFromFile */

/**
 * Assumed to be on the same set of variables
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
  const int numNB = numCols; // TODO again assuming that nb fixed vars are treated as ub vars

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
  double tinyObjOffset = 0.; // how much c . p changes when tiny values are ignored
  double nonTinyObj = 0.;

  // All coefficients are going to be non-negative, since we work in the complemented NB space
  for (int i = 0; i < numNB; i++) {
    const int currVar = nonBasicVarIndex[i];
    const int currVarTmpIndex = (currVar < deletedVar) ? currVar : currVar + deletedVarAdjustment;
    double tmpVal = 0., origVal = 0.;
    if (currVar < numCols) {
      tmpVal = tmpSolver->getColSolution()[currVarTmpIndex];
      origVal = origSolver->getColSolution()[currVar];
    } else {
      const int currRow = currVar - numCols;
      tmpVal = tmpSolver->getRightHandSide()[currRow]
          - tmpSolver->getRowActivity()[currRow];
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
    // I think the latter case we could have handled by preprocessing the free variables out
    // ****
    // For *slack* var:
    // In the original NB space at v, all slacks were at zero.
    // We simply look at b[row] - activity[row]
    // For >= rows, slack is non-positive
    // For <= rows, slack is non-negative
    // We take absolute value to deal with this
    const double newVal = std::abs(tmpVal - origVal);

    if (!isZero(newVal)) {
      packedIndex.push_back(i);
      packedElem.push_back(newVal);
      nonTinyObj += newVal * nonBasicReducedCost[i];
    } // non-tiny
    else if (!isZero(newVal * nonBasicReducedCost[i])) {
      tinyIndex.push_back(i);
      tinyElem.push_back(newVal);
      tinyObjOffset += newVal * nonBasicReducedCost[i];
    } // tiny
  } /* end iterating over non-basic elements from original basis */

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
    const int scale = 1.; // perhaps in a future version we will scale things
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
    if (!isVal(nonTinyObj, objViolation, 1e-3)) {
      error_msg(errorstring,
          "Point: Calculated dot product with obj differs from solver's. Obj viol from solver: %.8f. Calculated: %.8f.\n",
          objViolation, nonTinyObj);
      writeErrorToLog(errorstring, params.logfile);
      exit(1);
    } else {
      warning_msg(warnstring,
          "Point: Calculated dot product with obj differs from solver's. Obj viol from solver: %.8f. Calculated: %.8f.\n",
          objViolation, nonTinyObj);
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
  double tinyObjOffset = 0.; // how much c . r changes when tiny values are ignored
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
      if (!isZero(rayVal, params.get(RAYEPS))) {
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
      else if (!isZero(rayVal * nonBasicReducedCost[c], params.get(RAYEPS))) {
        tinyIndex.push_back(c);
        tinyElem.push_back(rayVal);
        tinyObjOffset += rayVal * nonBasicReducedCost[c];
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
    if (lessThanVal(nonTinyObj, 0., params.get(DIFFEPS))) {
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
    if (!isVal(nonTinyObj, objViolation, params.get(DIFFEPS))) {
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

void setCompNBCoor(CoinPackedVector& vec, double& objViolation,
    const VPCParameters& params, const double* const currColValue,
    const double* const currSlackValue,
    const OsiSolverInterface* const origSolver,
    const std::vector<int>& nonBasicVarIndex,
    const std::vector<double>& nonBasicReducedCost, const int deletedVar) {
  // If a variable was deleted to make tmpSolver, then we will need to adjust the index
  // It is -1 because we need the index in the *new* space
  const int deletedVarAdjustment = (deletedVar >= 0) ? -1 : 0;

  const int numCols = origSolver->getNumCols();
  const int numNB = numCols; // TODO again assuming that nb fixed vars are treated as ub vars

  // For keeping the packed coordinates to input as a new point
  std::vector<int> packedIndex;
  std::vector<double> packedElem;

  packedIndex.reserve(numNB);
  packedElem.reserve(numNB);

  objViolation = 0.;

  // All coefficients are going to be non-negative, since we work in the complemented NB space
  for (int i = 0; i < (int) nonBasicVarIndex.size(); i++) {
    const int currVar = nonBasicVarIndex[i];
    const int currVarTmpIndex = (currVar < deletedVar) ? currVar : currVar + deletedVarAdjustment;
    double tmpVal = 0., origVal = 0.;
    double objCoeff = 0.;
    if (currVar < numCols) {
      tmpVal = currColValue[currVarTmpIndex];
      origVal = origSolver->getColSolution()[currVar];
      objCoeff = origSolver->getObjCoefficients()[currVar];
    } else {
      const int currRow = currVar - numCols;
      tmpVal = currSlackValue[currRow];
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
    // TODO The latter case we have handled by preprocessing the free variables out
    // ****
    // For *slack* var:
    // In the original NB space at v, all slacks were at zero.
    // We simply look at b[row] - activity[row]
    // For >= rows, slack is non-positive
    // For <= rows, slack is non-negative
    // We take absolute value to deal with this
    const double newVal = std::abs(tmpVal - origVal);
    objViolation += newVal * nonBasicReducedCost[i];

    if (!isZero(newVal) || !isZero(newVal * objCoeff)) {
      packedIndex.push_back(i);
      packedElem.push_back(newVal);
    }
  }
  vec.setVector((int) packedIndex.size(), packedIndex.data(),
      packedElem.data(), false);
} /* setCompNBCoor (generic) */

/**
 * @brief Check if file exists
 */
bool fexists(const char* filename) {
  struct stat buffer;
  return (stat (filename, &buffer) == 0);
} /* fexists */

/**
 * @brief Parses int from string using strtol
 */
bool parseInt(const char *str, int& val) {
  long tmpval;
  bool rc = parseLong(str, tmpval);
  val = static_cast<int>(tmpval);
  return rc;
} /* parseInt */

/**
 * @brief Parses long int from string using strtol
 */
bool parseLong(const char *str, long& val) {
  char *temp;
  bool rc = true;
  errno = 0;
  val = strtol(str, &temp, 10);

  if (temp == str || *temp != '\0' ||
      ((val == std::numeric_limits<long>::min() || val == std::numeric_limits<long>::max()) && errno == ERANGE))
    rc = false;

  return rc;
} /* parseLong */

/**
 * @brief Parses double from string using strtod
 */
bool parseDouble(const char *str, double& val) {
  char *temp;
  bool rc = true;
  errno = 0;
  val = strtod(str, &temp);

  if (temp == str || *temp != '\0' ||
      ((val == std::numeric_limits<double>::min() || val == std::numeric_limits<double>::lowest()
        || val == std::numeric_limits<double>::max()) && errno == ERANGE))
    rc = false;

  return rc;
} /* parseDouble */

// Last edit: 03/27/12
//
// Name:     common_definitions.cpp
// Author:   Giacomo Nannicini
//           Singapore University of Technology and Design, Singapore
//           email: nannicini@sutd.edu.sg
// Date:     04/09/11
//-----------------------------------------------------------------------------
// Copyright (C) 2011, Giacomo Nannicini.  All Rights Reserved.

double dotProduct(const CoinPackedVector& vec1, const double* vec2) {
  const int size1 = vec1.getNumElements();
  const int* ind1 = vec1.getIndices();
  const double* el1 = vec1.getElements();
  return dotProduct(size1, ind1, el1, vec2);
}

double dotProduct(const CoinPackedVector& vec1, const CoinPackedVector& vec2) {
  const int size1 = vec1.getNumElements();
  const int* ind1 = vec1.getIndices();
  const double* el1 = vec1.getElements();
  const int size2 = vec2.getNumElements();
  const int* ind2 = vec2.getIndices();
  const double* el2 = vec2.getElements();
  return dotProduct(size1, ind1, el1, size2, ind2, el2);
}

/**********************************************************/
double dotProduct(const double* a, const double* b, int dimension) {
  // Accumulate the dot product here
  double accumulator = 0.0;
  // Compensation term for low-order bits lost
  double compensation = 0.0;
  // Temporary number for storing current term
  double currterm = 0.0;
  // Temporary number for storing new accumulator
  double nextacc = 0.0;
  for (int i = 0; i < dimension; ++i) {
    // Compute current term of the dot product adding compensation
    currterm = (a[i] * b[i]) - compensation;
    // This is the value of the sum
    nextacc = accumulator + currterm;
    // Recover what we just lost adding currterm to accumulator
    compensation = (nextacc - accumulator) - currterm;
    // Now save new value of the accumulator
    accumulator = nextacc;
  }
  return accumulator;
} /* dotProduct (dense x dense) */

/**********************************************************/
double dotProduct(int sizea, const int* indexa, const double* a,
    const double* b) {
  // Accumulate the dot product here
  double accumulator = 0.0;
  // Compensation term for low-order bits lost
  double compensation = 0.0;
  // Temporary number for storing current term
  double currterm = 0.0;
  // Temporary number for storing new accumulator
  double nextacc = 0.0;
  for (int i = 0; i < sizea; ++i) {
    // Compute current term of the dot product adding compensation
    currterm = (a[i] * b[indexa[i]]) - compensation;
    // This is the value of the sum
    nextacc = accumulator + currterm;
    // Recover what we just lost adding currterm to accumulator
    compensation = (nextacc - accumulator) - currterm;
    // Now save new value of the accumulator
    accumulator = nextacc;
  }
  return accumulator;
} /* dotProduct (sparse x dense) */

/**********************************************************/
double dotProduct(int sizea, const int* indexa, const double* a, int sizeb,
    const int* indexb, const double* b) {
  // Accumulate the dot product here
  double accumulator = 0.0;
  // Compensation term for low-order bits lost
  double compensation = 0.0;
  // Temporary number for storing current term
  double currterm = 0.0;
  // Temporary number for storing new accumulator
  double nextacc = 0.0;
  // Current position in vectors a and b
  int posa = 0, posb = 0;
  while (posa < sizea && posb < sizeb) {
    // If it is the same component, compute dot product
    if (indexa[posa] == indexb[posb]) {
      // Compute current term of the dot product adding compensation
      currterm = (a[posa] * b[posb]) - compensation;
      // This is the value of the sum
      nextacc = accumulator + currterm;
      // Recover what we just lost adding currterm to accumulator
      compensation = (nextacc - accumulator) - currterm;
      // Now save new value of the accumulator
      accumulator = nextacc;
      // Increment both position indices
      ++posa;
      ++posb;
    } else if (indexa[posa] < indexb[posb]) {
      // Increment only smaller position index
      ++posa;
    } else if (indexa[posa] > indexb[posb]) {
      // Increment only smaller position index
      ++posb;
    }
  }
  return accumulator;
} /* dotProduct (sparse x sparse) */
