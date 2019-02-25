// Name:     DisjunctionHelper.hpp
// Author:   A. M. Kazachkov
// Date:     2018-02-25
//-----------------------------------------------------------------------------
#pragma once

#include <vector>

class Disjunction;
class OsiSolverInterface;
struct VPCParameters;
enum class ExitReason;

#include "CglVPC.hpp" // ExitReason, VPCMode

ExitReason setDisjunctions(std::vector<Disjunction*>& disjVec,
    const OsiSolverInterface* const si, const VPCParameters& params,
    CglVPC::VPCMode mode);
int generateSplitDisjunctions(std::vector<Disjunction*>& disjVec,
    const OsiSolverInterface* const si, const VPCParameters& params);
