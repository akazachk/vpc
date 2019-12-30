// Name:     preprocess.hpp
// Author:   A. M. Kazachkov
// Date:     2019-03-16
//-----------------------------------------------------------------------------
#pragma once
#include <string>

class OsiSolverInterface;
namespace VPCParametersNamespace {
  struct VPCParameters;
}

void printPreprocessingHeader(const VPCParametersNamespace::VPCParameters& params, const char SEP = ',');
void performCleaning(const VPCParametersNamespace::VPCParameters& orig_params,
    OsiSolverInterface* const solver, const std::string& filename,
    const double ip_obj, const int CLEANING_MODE_OPTION, const char SEP = ',');
bool cleanProblem(const VPCParametersNamespace::VPCParameters& params, OsiSolverInterface* solver,
    int& numBoundsChanged, int& numSBFixed);
