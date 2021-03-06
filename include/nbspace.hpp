/**
 * @file nbspace.hpp
 * @brief Methods related to working in the nonbasic space
 *
 * @author A. M. Kazachkov
 * @date 2019-02-20
 */
#pragma once

#include <vector>

// COIN-OR
#include <OsiSolverInterface.hpp>

// VPC
namespace VPCParametersNamespace {
  struct VPCParameters;
}

/**
 * @brief Convert optimal solution to \p tmpSolver into a sparse point in the complemented nonbasic space
 *
 * Complemented space = nonbasic variables at their upper bound are flipped.
 * In other words, in this space, the LP solution is the origin.
 *
 * Assumed to be on the same set of variables.
 * We also need that the right-hand side of the original solver remains unchanged (except the bounds).
 */
void setCompNBCoorPoint(CoinPackedVector& vec, double& objViolation,
    const VPCParametersNamespace::VPCParameters& params,
    const OsiSolverInterface* const tmpSolver,
    const OsiSolverInterface* const origSolver,
    const std::vector<int>& nonBasicVarIndex,
    const std::vector<double>& nonBasicReducedCost, const int deletedVar = -1);

/**
 * @details Note that we may have deleted a variable to get to tmpSolver
 * In that case, tmpNBVar is in the space of tmpSolver
 */
void setCompNBCoorRay(CoinPackedVector& vec, const double* ray, double& objViolation, double& scale,
    const VPCParametersNamespace::VPCParameters& params, const OsiSolverInterface* const tmpSolver,
    const OsiSolverInterface* const origSolver,
    const std::vector<int>& rowOfOrigNBVar,
    const std::vector<int>& nonBasicVarIndex,
    const std::vector<double>& nonBasicReducedCost,
    const int tmpNBVar, const int deletedVar = -1,
    const bool rayNeedsCalculation = false);

/**
 * @brief Puts currColValue (necessary to provide) into complemented nonbasic space form
 *
 * Complemented space = nonbasic variables at their upper bound are flipped.
 * In other words, in this space, the LP solution is the origin.
 *
 * Assumed to be for the same set of variables as in originalSolver.
 * We also need that the right-hand side of the original solver remains unchanged (except the bounds).
 */
void setCompNBCoor(CoinPackedVector& vec, double& objViolation,
    const VPCParametersNamespace::VPCParameters& params,
    const double* const currColValue, const double* const currSlackValue,
    const OsiSolverInterface* const origSolver,
    const std::vector<int>& nonBasicVarIndex,
    const std::vector<double>& nonBasicReducedCost, const int deletedVar = -1);

/**
 * @brief Set current cut by converting it from nonbasic coeffs given in coeff
 */
void setCutFromNBCoefficients(OsiRowCut* const cut,
    const std::vector<int>& nonZeroColIndex, const double* const coeff,
    const double beta, const OsiSolverInterface* const solver,
    const std::vector<int> &nonBasicVarIndex, const bool inCompSpace,
    const double EPS);
