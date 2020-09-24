/**
 * @file preprocess.hpp
 * @brief Save preprocessed instance
 *
 * @author A. M. Kazachkov
 * @date 2019-03-16
 */
#pragma once
#include <string>

class OsiSolverInterface;
namespace VPCParametersNamespace {
  struct VPCParameters;
}

/// @brief Print header when preprocessing is enabled
void printPreprocessingHeader(const VPCParametersNamespace::VPCParameters& params, const char SEP = ',');

/// @brief Perform preprocessing and get statistics
void performCleaning(const VPCParametersNamespace::VPCParameters& orig_params,
    const OsiSolverInterface* const solver, const std::string& filename,
    const double ip_obj, const int CLEANING_MODE_OPTION, const char SEP = ',');

/// @brief Makes sure no variable bounds can be tightened, including via strong branching
bool cleanProblem(const VPCParametersNamespace::VPCParameters& params, OsiSolverInterface* solver,
    int& numBoundsChanged, int& numSBFixed);
