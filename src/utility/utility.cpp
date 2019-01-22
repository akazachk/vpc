#include "utility.hpp"

#include <fstream>
#include <sstream>
#include <sys/stat.h> // for fexists

/**
 * @brief Checks whether a solver is optimal
 *
 * Something that can go wrong (e.g., bc1 -64 sb5 tmp_ind = 14):
 *  Solver is declared optimal for the scaled problem but there are primal or dual infeasibilities in the unscaled problem
 *  In this case, secondary status will be 2 or 3
 * Sometimes cleanup helps
 * Sometimes things break after enabling factorization (e.g., secondary status = 0, but sum infeasibilities > 0)
 */
bool checkSolverOptimality(OsiSolverInterface* const solver,
    const bool exitOnDualInfeas, const double timeLimit,
    const int maxNumResolves) {
  OsiClpSolverInterface* clpsolver = NULL;
  try {
    clpsolver = dynamic_cast<OsiClpSolverInterface*>(solver);
  } catch (std::exception& e) {
    // Disregard / catch below
  }
  if (!clpsolver) {
    return solver->isProvenOptimal();
  }

  // If it is supposedly proven primal infeasible, might be good to do a resolve first
  // This is, e.g., something that arises with miplib2003/pp08aCUTS_presolved -32
  int resolve_count = 0;
  if (clpsolver->isProvenPrimalInfeasible()) {
    clpsolver->getModelPtr()->setNumberIterations(0);
    clpsolver->getModelPtr()->setMaximumSeconds(timeLimit);
    clpsolver->resolve();
    resolve_count++;
  }

  // First clean
  if (resolve_count < maxNumResolves && clpsolver) {
    int status = clpsolver->getModelPtr()->secondaryStatus();
    bool is_cleaned = false;
    if (status == 2 || !isZero(clpsolver->getModelPtr()->sumPrimalInfeasibilities())) {
      clpsolver->getModelPtr()->cleanup(1);
      is_cleaned = true;
    } else if (status == 3 || status == 9 || !isZero(clpsolver->getModelPtr()->sumDualInfeasibilities())) {
      clpsolver->getModelPtr()->cleanup(2);
      is_cleaned = true;
    } else if (status == 4) {
      clpsolver->getModelPtr()->cleanup(1);
      is_cleaned = true;
    }
    if (is_cleaned) {
      resolve_count++;
      return checkSolverOptimality(solver, exitOnDualInfeas, timeLimit, 0);
    }

    // Do some resolves, but first save whether the initial status is dual infeasible
    // The reason is that for neos15 w/str=2, we were getting a dual infeasible problem
    // turn into a primal infeasible problem after resolves
    // (probably the strengthening causes numerical issues)
    const bool oldStatusDualInfeasible = (clpsolver->isProvenDualInfeasible());
    const bool oldStatusOptimal = (clpsolver->isProvenOptimal());
    bool resolve = true;
    double infeas = std::numeric_limits<double>::max();
    while (resolve && resolve_count < maxNumResolves) {
      const double curr_infeas =
          clpsolver->getModelPtr()->sumPrimalInfeasibilities()
              + clpsolver->getModelPtr()->sumDualInfeasibilities();
      resolve = (resolve_count == 0)
          || (!isZero(curr_infeas) && (curr_infeas < infeas));
      if (resolve) {
        clpsolver->getModelPtr()->setNumberIterations(0);
        clpsolver->getModelPtr()->setMaximumSeconds(timeLimit);
        clpsolver->resolve();
        resolve_count++;
      }
      infeas = curr_infeas;
    }

    // If dual infeas -> primal infeas, do a dual resolve, which should fix the issue
    if (oldStatusDualInfeasible && clpsolver->isProvenPrimalInfeasible()) {
      clpsolver->getModelPtr()->dual(2,0); // just do values pass
      if (!clpsolver->isProvenDualInfeasible() && !clpsolver->isProvenOptimal()) {
        clpsolver->getModelPtr()->setProblemStatus(2);
        resolve_count = maxNumResolves; // stop resolving
      }
    }

    if (oldStatusOptimal && clpsolver->isProvenPrimalInfeasible()) {
      clpsolver->enableFactorization();
      clpsolver->getModelPtr()->setNumberIterations(0);
      clpsolver->getModelPtr()->setMaximumSeconds(timeLimit);
      clpsolver->resolve();
      resolve_count++;
      clpsolver->disableFactorization();
    }

    // Clean once more if needed and possible
    status = clpsolver->getModelPtr()->secondaryStatus();
    is_cleaned = false;
    if (status == 2
        || !isZero(clpsolver->getModelPtr()->sumPrimalInfeasibilities())) {
      clpsolver->getModelPtr()->cleanup(1);
      is_cleaned = true;
    } else if (status == 3 || status == 9
        || !isZero(clpsolver->getModelPtr()->sumDualInfeasibilities())) {
      clpsolver->getModelPtr()->cleanup(2);
      is_cleaned = true;
    } else if (status == 4) {
      clpsolver->getModelPtr()->cleanup(1);
      is_cleaned = true;
    }
    if (is_cleaned) {
      resolve_count++;
      return checkSolverOptimality(solver, exitOnDualInfeas, timeLimit, 0);
    }
  }

  if (solver->isProvenPrimalInfeasible()) {
    return false;
  } else if (!(solver->isProvenOptimal())) {
    // Sometimes need to resolve once more to get the correct status
    if (resolve_count < maxNumResolves) {
      solver->resolve();
    }
    if (solver->isProvenPrimalInfeasible()) {
      return false;
    } else if (solver->isProvenDualInfeasible()) {
      if (exitOnDualInfeas) {
        error_msg(errstr,
            "Solver is dual infeasible. Check why this happened!\n");
        exit(1);
      } else {
        return false;
      }
    } else if (!(solver->isProvenOptimal())) {
      return false;
    }
  }

  return true;
} /* checkSolverOptimality */

/**
 * @brief Enable factorization, and check whether some cleanup needs to be done
 */
bool enableFactorization(OsiSolverInterface* const solver, const double EPS, const int resolveFlag) {
  solver->enableFactorization();

  if (resolveFlag > 0) {
    try {
      OsiClpSolverInterface* clpsolver =
          dynamic_cast<OsiClpSolverInterface*>(solver);

      // After enabling factorization, things sometimes look different and it is worth resolving
      // This is motivated by seymour-disj-10; after adding one round of GMICs, row 5077 had negative row price
      if ((clpsolver->getModelPtr()->sumPrimalInfeasibilities() > EPS)
          || (clpsolver->getModelPtr()->sumDualInfeasibilities() > EPS)) {
        if (resolveFlag == 1) {
          clpsolver->initialSolve(); // hopefully this fixes things rather than breaks things; actually it breaks things, e.g., for rd2 SICs for coral/neos17, the variable basic in row 648 changes from 164 to 23
        } else {
          clpsolver->resolve();
        }
        clpsolver->disableFactorization();
        clpsolver->enableFactorization(); // this will hopefully solve the neos17 problem
      }
    } catch (std::exception& e) {
      // Disregard
    }
  }

  return solver->isProvenOptimal();
} /* enableFactorization */

void setLPSolverParameters(OsiSolverInterface* const solver,
    const double max_time) {
#ifndef TRACE
  solver->messageHandler()->setLogLevel(0);
  try {
    dynamic_cast<OsiClpSolverInterface*>(solver)->getModelPtr()->messageHandler()->setLogLevel(0);
  } catch (std::exception& e) {

  }
#endif
  try {
    dynamic_cast<OsiClpSolverInterface*>(solver)->getModelPtr()->setMaximumSeconds(max_time);
  } catch (std::exception& e) {

  }

  // Try turning on scaling with enableFactorization
//  solver->setSpecialOptions(solver->specialOptions() | 512);
} /* setLPSolverParameters (OsiClp) */

void setIPSolverParameters(CbcModel* const cbc_model) {
#ifdef TRACE
  cbc_model->setLogLevel(3);
  cbc_model->messagesPointer()->setDetailMessages(10, 10000, (int *) NULL);
#else
  cbc_model->setLogLevel(0);
  cbc_model->messagesPointer()->setDetailMessages(10,5,5000);
#endif
  if (cbc_model->solver()) {
    cbc_model->solver()->setHintParam(OsiDoReducePrint, true, OsiHintTry);
  }
  cbc_model->setPrintFrequency(1);
} /* setIPSolverParameters (Cbc) */

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
void setCompNBCoor(CoinPackedVector& vec, double& objViolation,
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
} /* setCompNBCoor */

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
} /* setCompNBCoor */

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
