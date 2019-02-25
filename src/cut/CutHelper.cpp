// Name:     CutHelper.cpp
// Author:   A. M. Kazachkov
// Date:     2019-Feb-20
//-----------------------------------------------------------------------------

#include <numeric> // inner_product, sqrt

#include "CutHelper.hpp"
#include "CglVPC.hpp" // For FailureType
#include "SolverHelper.hpp"
#include "VPCParameters.hpp"

/***********************************************************************/
/**
 * @brief Set current cut by converting it from non-basic coeffs given in coeff
 */
void setCutFromJSpaceCoefficients(OsiRowCut* const cut,
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
  } /* iterate over columns */

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

/**
 * The cut is stored as alpha x >= beta.
 */
void setOsiRowCut(OsiRowCut* const cut, const std::vector<int>& nonZeroColIndex,
    const int num_coeff, const double* coeff, const double rhs,
    const double EPS) {
  std::vector<int> indices;
  std::vector<double> sparse_coeff;
  indices.reserve(num_coeff);
  sparse_coeff.reserve(num_coeff);

  for (int i = 0; i < num_coeff; i++) {
    if (!isZero(coeff[i], EPS)) {
      indices.push_back((nonZeroColIndex.size() > 0) ? nonZeroColIndex[i] : i);
      sparse_coeff.push_back(coeff[i]);
    }
  }
  cut->setLb(rhs);
  cut->setRow((int) indices.size(), indices.data(),
      sparse_coeff.data());
} /* setOsiRowCut */

/**
 * Set the cut solver objective coefficients based on a packed vector
 * We do not zero out the other coefficients unless requested
 */
void addToObjectiveFromPackedVector(OsiSolverInterface* const cutSolver,
    const CoinPackedVectorBase* vec, const bool zeroOut, const double mult) {
  if (zeroOut) {
    for (int i = 0; i < cutSolver->getNumCols(); i++) {
      cutSolver->setObjCoeff(i, 0.);
    }
  }

  if (vec) {
    const double twoNorm = vec->twoNorm() / cutSolver->getNumCols();
    for (int i = 0; i < (*vec).getNumElements(); i++) {
      const int col = (*vec).getIndices()[i];
      const double elem = (*vec).getElements()[i];
      const double coeff = cutSolver->getObjCoefficients()[col];
      cutSolver->setObjCoeff(col, coeff + elem * mult / twoNorm);
    }
  }
} /* addToObjectiveFromPackedVector */

void setConstantObjectiveFromPackedVector(
    OsiSolverInterface* const cutSolver, const double val,
    const int numIndices, const int* indices) {
  if (numIndices > 0 && indices) {
    for (int i = 0; i < numIndices; i++) {
      cutSolver->setObjCoeff(indices[i], val);
    }
  } else {
    for (int i = 0; i < cutSolver->getNumCols(); i++) {
      cutSolver->setObjCoeff(i, val);
    }
  }
} /* setConstantObjectiveFromPackedVector */

/**
 * Taken from CglGMI with some minor modifications, including switching to >= cut as we use
 */
void removeSmallCoefficients(OsiRowCut* const cut, const OsiSolverInterface* const solver, const double EPS, const double EPS_COEFF) {
  CoinPackedVector vec = cut->mutableRow();
  int cutNz = vec.getNumElements();
  int* cutIndex = vec.getIndices();
  double* cutElem = vec.getElements();
  double cutRhs = cut->rhs();
  const double* colLower = solver->getColLower();
  const double* colUpper = solver->getColUpper();

  double coeff, absval;
  int currPos = 0;
  int col;
  for (int i = 0; i < cutNz; ++i) {
    col = cutIndex[i];
    coeff = cutElem[i];
    absval = std::abs(coeff);
    if (isZero(absval, EPS)
        && isZero(absval * solver->getObjCoefficients()[col], EPS)) {
      continue; // throw away this element
    }
    if (absval <= EPS_COEFF) {
      // small coefficient: remove and adjust rhs if possible
      if ((coeff < 0.0) && !isNegInfinity(colLower[col])) {
        cutRhs -= coeff * colLower[col];
      } else if ((coeff > 0.0) && !isInfinity(colUpper[col])) {
        cutRhs -= coeff * colUpper[col];
      }
    } else if (absval > EPS_COEFF) {
      if (currPos < i) {
        cutElem[currPos] = cutElem[i];
        cutIndex[currPos] = cutIndex[i];
      }
      currPos++;
    }
  }
  cutNz = currPos;
  cut->setRow(cutNz, cutIndex, cutElem);
  if (isZero(cutRhs, EPS)) {
    cutRhs = 0.;
  }
  cut->setLb(cutRhs);
} /* removeSmallCoefficients */

bool badSupport(const int cutNz, const int numCols, const double max_sup_abs, const double max_sup_rel) {
  if (cutNz > max_sup_abs + max_sup_rel * numCols) {
    return true;
  }
  return false;
} /* badSupport */

bool badViolation(const OsiRowCut* const cut, const OsiSolverInterface* const solver, const double min_viol_abs, const double min_viol_rel) {
  // Check the cut
  const double currCutNorm = cut->row().twoNorm();
  const double violation = cut->violated(solver->getColSolution());
  //    const bool cuttingSolutionFlagAbs = (inNBSpace)
  //        ? lessThanVal(cutSolver->getObjValue() - beta, 0.0, param.getMIN_VIOL_ABS())
  //        : greaterThanVal(currCut.violated(origSolver->getColSolution()), 0.0, param.getMIN_VIOL_ABS());
  //    const bool cuttingSolutionFlagRel = (inNBSpace)
  //        ? lessThanVal(cutSolver->getObjValue() - beta, 0.0, currCutNorm * param.getMIN_VIOL_REL())
  //        : greaterThanVal(currCut.violated(origSolver->getColSolution()), 0.0, currCutNorm * param.getMIN_VIOL_REL());
  const bool cuttingSolutionFlagAbs = //(inNBSpace) ||
      greaterThanVal(violation, 0.0, min_viol_abs);
  const bool cuttingSolutionFlagRel = //(inNBSpace) ||
      greaterThanVal(violation, 0.0, currCutNorm * min_viol_rel);
  if (!cuttingSolutionFlagAbs || !cuttingSolutionFlagRel) {
    return true;
  }
  return false;
} /* badViolation */

bool badDynamism(const OsiRowCut* const cut, const double minAbsCoeff, const double maxAbsCoeff, const double SENSIBLE_MAX, const double EPS) {
  const CoinPackedVector vec = cut->row();
  const int num_el = vec.getNumElements();
  const double* el = vec.getElements();

  double minAbsElem = 0.0, maxAbsElem = 0.0;
  for (int i = 0; i < num_el; i++) {
    const double absCurr = std::abs(el[i]);
    if (isZero(absCurr, EPS)) { // make sure min and max are not zero
      continue;
    }

    if ((minAbsElem > absCurr) || isZero(minAbsElem, EPS)) {
      minAbsElem = absCurr;
    }
    if (maxAbsElem < absCurr) {
      maxAbsElem = absCurr;
    }
  }
  const double absCurr = std::abs(cut->rhs());
  if (absCurr > maxAbsElem) {
    maxAbsElem = absCurr;
  }
  if (isZero(maxAbsElem, EPS)
      // || isZero(minAbsElem / maxAbsElem)
      || greaterThanVal(maxAbsElem / minAbsElem, SENSIBLE_MAX, EPS)
      || (!isZero(minAbsCoeff, EPS)
          && greaterThanVal(std::abs(maxAbsElem / minAbsCoeff),
              std::abs(SENSIBLE_MAX * maxAbsCoeff / minAbsCoeff), EPS))) {
    return true;
  }
  return false;
} /* badDynamism */

/**
 * Based on the similar method in CglGMI, we clean the cut coefficients and check if the cut is good
 * Returns 0 if no error, otherwise returns -1 * (fail index + 1)
 */
int cleanCut(OsiRowCut* const cut, const OsiSolverInterface* const solver,
    const VPCParameters& params, const double min_abs_coeff,
    const double max_abs_coeff, const double EPS, const double beta,
    const bool checkViolation) {
  const double EPS_COEFF = params.get(doubleConst::EPS_COEFF);
  const double MAX_DYN = params.get(doubleConst::MAX_DYN);
  const int MAX_SUP_ABS = params.get(intConst::MAX_SUPPORT_ABS);
  const double MAX_SUP_REL = params.get(doubleConst::MAX_SUPPORT_REL);
  const double MIN_VIOL_ABS = params.get(doubleConst::MIN_VIOL_ABS);
  const double MIN_VIOL_REL = params.get(doubleConst::MIN_VIOL_REL);

  removeSmallCoefficients(cut, solver, EPS, EPS_COEFF);
  if (badSupport(cut->row().getNumElements(), solver->getNumCols(), MAX_SUP_ABS, MAX_SUP_REL)) {
    return -1 * (static_cast<int>(CglVPC::FailureType::BAD_SUPPORT) + 1);
  }
  // checkViolation = !inNBSpace && !lessThanVal(beta, 0.) 
  if (checkViolation && badViolation(cut, solver, MIN_VIOL_ABS, MIN_VIOL_REL)) {
    return -1 * (static_cast<int>(CglVPC::FailureType::BAD_VIOLATION) + 1);
  }
  if (badDynamism(cut, min_abs_coeff, max_abs_coeff,
      CoinMax((EPS > 0) ? 1. / EPS : MAX_DYN, MAX_DYN), EPS)) {
    return -1 * (static_cast<int>(CglVPC::FailureType::BAD_DYNAMISM) + 1);
  }
  return 0;
} /* cleanCut */

/**
 * @brief Decide if two rows are the same.
 * @return
 *   0: seem same in coeff and rhs,
 *   +/-1: seem same in coeff, but diff in rhs (-1: cut1 better, +1: cut2 better),
 *   2: seem different in coeff and rhs.
 *
 * The value of eps will be used to determine whether a cut coefficient is zero or not.
 * We are roughly checking that whether cut1 = ratio * cut2 for some ratio.
 * Let ratio = rhs1 / rhs2. (If both are non-zero. Typically both are = 1.)
 * We say the cuts are "truly different" if |coeff^1_j - ratio * coeff^2_j| >= diffeps.
 * However, might be that the cut coefficients are scaled versions, but the rhs values are different.
 * So we also compute ratio_coeff = coeff^1_j / coeff^2_j for the first j where both are non-zero.
 * If |coeff^1_j - ratio * coeff^2_j| >= diffeps for some j,
 * but |coeff^1_j - ratio_coeff * coeff^2_j| < diffeps for all j,
 * then the return code will be -1 or 1 (essentially, one cut will dominate the other).
 *
 * By the way, this is all essentially checking orthogonality...
 */
int isRowDifferent(const CoinPackedVector& cut1Vec, const double cut1rhs,
    const CoinPackedVector& cut2Vec, const double cut2rhs, const double EPS) {
  const int numElem1 = cut1Vec.getNumElements();
  const int* index1 = cut1Vec.getIndices();
  const double* value1 = cut1Vec.getElements();

  const int numElem2 = cut2Vec.getNumElements();
  const int* index2 = cut2Vec.getIndices();
  const double* value2 = cut2Vec.getElements();

  double ratio = -1.0;
  double ratio_coeff = -1.0;
  double maxDiff = 0.;
  double maxDiffCoeff = 0.;
  double maxDiffZero = 0.;

  // First, we try to get the ratio from RHS
  if (isZero(cut1rhs, EPS) || isZero(cut2rhs, EPS)) {
    // If one but not both are 0, then the cuts are clearly different
    // But they may only differ in the rhs
    if (!isZero(cut1rhs, EPS)) {
      maxDiff = std::abs(cut1rhs);
    } else if (!isZero(cut2rhs, EPS)) {
      maxDiff = std::abs(cut2rhs);
    }
    // Otherwise, they're both near zero,
    // and we're not going to be able to get
    // any kind of ratio from here
  } else {
    ratio = cut1rhs / cut2rhs;
  }

  // Now check if this ratio (if we found one) holds for
  // all the coefficients. In particular, we check that
  // std::abs(cut1.coeff - ratio * cut2.coeff) < eps
  int ind2 = 0;
  for (int ind1 = 0; ind1 < numElem1; ind1++) {
    // Find matching indices
    while (ind2 < numElem2 && index1[ind1] > index2[ind2]) {
      if (!isZero(value2[ind2], EPS)) {
        return 2;
      }
      if (value2[ind2] > maxDiffZero) {
        maxDiffZero = std::abs(value2[ind2]);
      }
      ind2++;
    }

    // If we reached end of the second cut,
    // or index1[ind1] does not exist in the second cut,
    // capture maxDiffZero, or exit early if something "really different" detected
    if (ind2 >= numElem2 || index2[ind2] > index1[ind1]) {
      if (!isZero(value1[ind1], EPS)) {
        return 2;
      }
      if (value1[ind1] > maxDiffZero) {
        maxDiffZero = std::abs(value1[ind1]);
      }
      continue;
    }

    // Either one or both of the coefficients may be zero,
    // in which case both better be zero.
    if (isZero(value1[ind1], EPS) || isZero(value2[ind2], EPS)) {
      if (!isZero(value1[ind1], EPS) || !isZero(value2[ind2], EPS)) {
        return 2; // truly different
      }
      ind2++;
      continue;
    }

    // Alternatively, we make sure inequality is not scaled version of other
    // Set scale factor if it has not been set before
    if (lessThanVal(ratio, 0.0)) {
      ratio = value1[ind1] / value2[ind2];
      if (lessThanVal(ratio, 0.0)) {
        // The scale is negative
        // Potentially one cut is the negative of the other,
        // but they are thus still distinct
        // (all cuts are ax >= b)
        return 2;
      }
    } else {
      const double diff = std::abs(value1[ind1] - ratio * value2[ind2]);
      if (diff > maxDiff) {
        maxDiff = diff;
      }
    }

    // Also check the ratio_coeff in case one cut dominates the other
    if (lessThanVal(ratio_coeff, 0.0)) {
      ratio_coeff = value1[ind1] / value2[ind2];
      if (lessThanVal(ratio_coeff, 0.0)) {
        // The scale is negative
        // Potentially one cut is the negative of the other,
        // but they are thus still distinct
        // (all cuts are ax >= b)
        return 2;
      }
    } else {
      const double diff = std::abs(value1[ind1] - ratio_coeff * value2[ind2]);
      if (diff > maxDiffCoeff) {
        maxDiffCoeff = diff;
      }
    }

    if (!isZero(maxDiff, EPS) && !isZero(maxDiffCoeff, EPS)) {
      return 2;
    }

    ind2++;
  } /* iterate over cut1 coefficients */

  // Check any remaining cut2 indices that were not matched
  while (ind2 < numElem2) {
    if (!isZero(value2[ind2], EPS)) {
      return 2;
    }
    if (value2[ind2] > maxDiffZero) {
      maxDiffZero = std::abs(value2[ind2]);
    }
    ind2++;
  }

  // Are the cuts scaled versions of one another?
  if (isZero(maxDiff, EPS)) {
    return 0; // they seem the same
  }

  // Are the cut coefficients scaled versions of one another?
  if (isZero(maxDiffCoeff, EPS)) {
    if (lessThanVal(ratio_coeff, 0., EPS)) {
      return 2; // cuts face opposite directions
    }

    // Either one or the other rhs will be better
    const double rhs_diff = cut1rhs - ratio_coeff * cut2rhs;
    if (isZero(rhs_diff, EPS)) {
      return 0; // they seem the same, not sure how ratio missed it
    } else if (greaterThanVal(rhs_diff, 0., EPS)) {
      return -1;
    } else {
      return 1;
    }
  }

  return 0;
} /* isRowDifferent */

/**
 * Two vectors are parallel iff u.v/|u|*|v| = 1
 */
double getParallelism(const CoinPackedVectorBase& vec1,
    const CoinPackedVectorBase& vec2) {
  const double scalarprod = dotProduct(vec1, vec2);
  return scalarprod / (vec1.twoNorm() * vec2.twoNorm());
} /* getParallelism (packed vectors) */

/**
 * Two vectors are parallel iff u.v/|u|*|v| = 1
 */
double getParallelism(const CoinPackedVectorBase& vec1, const int numElem,
    const double* vec2) {
  const double scalarprod = vec1.dotProduct(vec2);
  const double norm1 = vec1.twoNorm();
  const double norm2 = std::sqrt(std::inner_product(vec2, vec2 + numElem, vec2, 0.0L));
  return scalarprod / (norm1 * norm2);
} /* getParallelism (packed and not packed) */

/**
 * Two vectors are parallel iff u.v/|u|*|v| = 1
 */
double getParallelism(const int numElem, const double* vec1,
    const double* vec2) {
  const double scalarprod = dotProduct(vec1, vec2, numElem);
  const double norm1 = std::sqrt(std::inner_product(vec1, vec1 + numElem, vec1, 0.0L));
  const double norm2 = std::sqrt(std::inner_product(vec2, vec2 + numElem, vec2, 0.0L));
  return scalarprod / (norm1 * norm2);
} /* getParallelism (not packed) */

/**
 * Two vectors are parallel iff u.v/|u|*|v| = 1
 */
double getOrthogonality(const CoinPackedVectorBase& vec1,
    const CoinPackedVectorBase& vec2) {
  return 1. - getParallelism(vec1, vec2);
} /* getOrthogonality (packed vectors) */

/**
 * Two vectors are parallel iff u.v/|u|*|v| = 1
 */
double getOrthogonality(const CoinPackedVectorBase& vec1, const int numElem,
    const double* vec2) {
  return 1. - getParallelism(vec1, numElem, vec2);
} /* getOrthogonality (packed and not packed) */

/**
 * Two vectors are parallel iff u.v/|u|*|v| = 1
 */
double getOrthogonality(const int numElem, const double* vec1,
    const double* vec2) {
  return 1. - getParallelism(numElem, vec1, vec2);
} /* getOrthogonality (not packed) */

/** 
 * Check whether a cut is duplicate or too orthogonal to a previous cut in the collection
 */
int howDuplicate(const OsiCuts& cuts, const OsiRowCut& tmpCut, int& duplicateCutIndex,
    int& minOrthoIndex, double& minOrtho, const double MIN_ORTHO, const double EPS) {
  duplicateCutIndex = -1;
  minOrthoIndex = -1;
  minOrtho = 2.;
  for (int i = 0; i < cuts.sizeCuts(); i++) {
    const OsiRowCut* currCut = cuts.rowCutPtr(i);
    if (MIN_ORTHO > 0) {
      const double this_ortho = getOrthogonality(tmpCut.row(), currCut->row());
      if (this_ortho < minOrtho) {
        minOrtho = this_ortho;
        minOrthoIndex = i;
      }
    }

    const int howDifferent = isRowDifferent(currCut->row(), currCut->rhs(), tmpCut.row(), tmpCut.rhs(), EPS);
    if (howDifferent != 2) {
      duplicateCutIndex = i;
      return howDifferent;
    }
  }
  return 0;
} /* howDuplicate */

