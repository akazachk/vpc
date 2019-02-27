// Name:     CutHelper.hpp
// Author:   A. M. Kazachkov
// Date:     2019-Feb-20
//-----------------------------------------------------------------------------
#pragma once

#include <OsiRowCut.hpp>
#include <OsiCuts.hpp>
#include <OsiSolverInterface.hpp>

struct VPCParameters; // defined in VPCParameters.hpp

/**
 * Generic way to apply a set of cuts to a solver
 * @return Solver value after adding the cuts, if they were added successfully
 */
double applyCutsCustom(OsiSolverInterface* solver, const OsiCuts& cs,
    const int startCutIndex, const int numCutsOverload, double compare_obj);
double applyCutsCustom(OsiSolverInterface* solver, const OsiCuts& cs,
    const int numCutsOverload = -1,
    double compare_obj = std::numeric_limits<double>::lowest());

void setOsiRowCut(OsiRowCut* const cut, const std::vector<int>& nonZeroColIndex,
    const int num_coeff, const double* coeff, const double rhs,
    const double EPS);

void addToObjectiveFromPackedVector(OsiSolverInterface* const solver,
    const CoinPackedVectorBase* vec, const bool zeroOut = false,
    const double mult = 1.);
void setConstantObjectiveFromPackedVector(OsiSolverInterface* const solver,
    const double val = 0., const int numIndices = 0, const int* indices = NULL);

int cleanCut(OsiRowCut* const cut, const OsiSolverInterface* const solver,
    const VPCParameters& params, const double min_abs_coeff,
    const double max_abs_coeff, const double EPS, const double beta,
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
