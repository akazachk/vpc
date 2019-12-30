/**
 * @file CutHelper.hpp
 * @author A. M. Kazachkov
 * @date 2019-02-20
 */
#pragma once
#include <limits>
#include <vector>

class CoinPackedVector;
class CoinPackedVectorBase;
class OsiRowCut;
class OsiCuts;
class OsiSolverInterface;

void setOsiRowCut(OsiRowCut* const cut, const std::vector<int>& nonZeroColIndex,
    const int num_coeff, const double* coeff, const double rhs,
    const double EPS);

void addToObjectiveFromPackedVector(OsiSolverInterface* const solver,
    const CoinPackedVectorBase* vec, const bool zeroOut = false,
    const double mult = 1., const std::vector<int>* const nonZeroColIndices = 0);
void setConstantObjectiveFromPackedVector(OsiSolverInterface* const solver,
    const double val = 0., const int numIndices = 0, const int* indices = 0);

int cleanCut(OsiRowCut* const cut, const OsiSolverInterface* const solver,
    const double EPS_COEFF,
    const double MAX_DYN,
    const int MAX_SUP_ABS,
    const double MAX_SUP_REL,
    const double MIN_VIOL_ABS,
    const double MIN_VIOL_REL,
    const double MIN_ABS_COEFF, 
    const double MAX_ABS_COEFF, 
    const double EPS, 
    const bool checkViolation);

/** Decide if two rows are the same */
int isRowDifferent(const CoinPackedVector& cut1Vec, const double cut1rhs,
    const CoinPackedVector& cut2Vec, const double cut2rhs, const double EPS);

/** Determine parallelism; two vectors are parallel iff u.v/|u|*|v| = 1 */
double getParallelism(const CoinPackedVectorBase& vec1, const CoinPackedVectorBase& vec2);
double getParallelism(const CoinPackedVectorBase& vec1, const int numElem, const double* vec2);
double getParallelism(const int numElem, const double* vec1, const double* vec2);
double getOrthogonality(const CoinPackedVectorBase& vec1, const CoinPackedVectorBase& vec2);
double getOrthogonality(const CoinPackedVectorBase& vec1, const int numElem, const double* vec2);
double getOrthogonality(const int numElem, const double* vec1, const double* vec2);

/** Check whether a cut is duplicate or too orthogonal to a previous cut in the collection */
int howDuplicate(const OsiCuts& cuts, const OsiRowCut& tmpCut,
    const int startIndex, int& duplicateCutIndex, int& minOrthoIndex,
    double& minOrtho, const double MIN_ORTHO, const double EPS);
