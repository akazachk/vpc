/**
 * @file DisjunctionHelper.hpp
 * @author A. M. Kazachkov
 * @date 2018-02-25
 */
#pragma once

#include <vector>

class Disjunction;
class OsiSolverInterface;
namespace VPCParametersNamespace {
  struct VPCParameters;
}

#include "CglVPC.hpp" // ExitReason, VPCMode

/// Set disjunctions; if integer-optimal solution is found, delete all but one disjunction, which will have that solution
CglVPC::ExitReason setDisjunctions(std::vector<Disjunction*>& disjVec,
    const OsiSolverInterface* const si, const VPCParametersNamespace::VPCParameters& params,
    CglVPC::VPCMode mode);

/// @brief Generate SplitDisjunction for all variables that are fractional, using hot starts
int generateSplitDisjunctions(std::vector<Disjunction*>& disjVec,
    const OsiSolverInterface* const si, const VPCParametersNamespace::VPCParameters& params);
