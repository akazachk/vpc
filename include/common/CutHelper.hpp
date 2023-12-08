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

/// @brief 
void setOsiRowCut(OsiRowCut* const cut, const std::vector<int>& nonZeroColIndex,
    const int num_coeff, const double* coeff, const double rhs,
    const double EPS);

/** @brief Clean cut coefficients and check if the cut is good */
int cleanCut(OsiRowCut* const cut,
    const OsiSolverInterface* const solver,
    const double EPS_COEFF,
    const double MAX_DYN,
    const int MIN_SUP_THRESH,
    const int MAX_SUP_ABS,
    const double MAX_SUP_REL,
    const double MIN_VIOL_ABS,
    const double MIN_VIOL_REL,
    const double MIN_ABS_COEFF, 
    const double MAX_ABS_COEFF, 
    const double EPS, 
    const bool checkViolation);

/** @brief Check whether the LP solution is cut away by a minimum amount **/
bool badViolation(const OsiRowCut* const cut, const OsiSolverInterface* const solver, const double min_viol_abs, const double min_viol_rel);

/** @brief Decide if two rows are the same */
int isRowDifferent(const CoinPackedVectorBase* const cut1Vec, const double cut1rhs,
    const CoinPackedVectorBase* const cut2Vec, const double cut2rhs, const double EPS);

/** @brief Determine parallelism; two vectors are parallel iff u.v/|u|*|v| = 1 */
double getParallelism(const CoinPackedVectorBase& vec1, const CoinPackedVectorBase& vec2);
/** @brief Determine parallelism; two vectors are parallel iff u.v/|u|*|v| = 1 */
double getParallelism(const CoinPackedVectorBase& vec1, const int numElem, const double* vec2);
/** @brief Determine parallelism; two vectors are parallel iff u.v/|u|*|v| = 1 */
double getParallelism(const int numElem, const double* vec1, const double* vec2);
/// @brief Determine orthogonality; 1 - #getParallelism(const CoinPackedVectorBase&, const CoinPackedVectorBase&)
double getOrthogonality(const CoinPackedVectorBase& vec1, const CoinPackedVectorBase& vec2);
/// @brief Determine orthogonality; 1 - #getParallelism(const CoinPackedVectorBase&, const int, const double*)
double getOrthogonality(const CoinPackedVectorBase& vec1, const int numElem, const double* vec2);
/// @brief Determine orthogonality; 1 - #getParallelism(const int, const double*, const double*)
double getOrthogonality(const int numElem, const double* vec1, const double* vec2);

/** @brief Check whether a cut is duplicate or too orthogonal to a previous cut in the collection */
int howDuplicate(const OsiCuts& cuts, const OsiRowCut& tmpCut,
    const int startIndex, int& duplicateCutIndex, int& minOrthoIndex,
    double& minOrtho, const double MIN_ORTHO, const double EPS);
