/**
 * @file CbcHelper.hpp
 * @author A. M. Kazachkov
 * @date 2021-Oct-3
 */
#pragma once

#include <limits>
#include <string>
#include <vector>

class OsiSolverInterface;
class OsiCuts;

struct BBInfo;
namespace VPCParametersNamespace {
  struct VPCParameters;
}

#ifdef USE_CBC
void presolveModelWithCbc(const VPCParametersNamespace::VPCParameters& params, int strategy,
    const char* f_name, double& presolved_lp_opt, std::string& presolved_name,
    const double best_bound);
void presolveModelWithCbc(const VPCParametersNamespace::VPCParameters& params, int strategy,
    const OsiSolverInterface* const solver, double& presolved_lp_opt, std::string& presolved_name);
void doBranchAndBoundWithCbc(const VPCParametersNamespace::VPCParameters& params, int strategy, 
    const char* f_name, BBInfo& info,
    const double best_bound = std::numeric_limits<double>::max(),
    std::vector<double>* const solution = NULL);
void doBranchAndBoundWithCbc(const VPCParametersNamespace::VPCParameters& params, int strategy,
    const OsiSolverInterface* const solver, BBInfo& info,
    const double best_bound = std::numeric_limits<double>::max(),
    std::vector<double>* const solution = NULL);
void doBranchAndBoundWithUserCutsCbc(const VPCParametersNamespace::VPCParameters& params, int strategy,
    const char* f_name, const OsiCuts* cuts, BBInfo& info,
    const double best_bound, const bool addAsLazy = false);
void doBranchAndBoundWithUserCutsCbc(const VPCParametersNamespace::VPCParameters& params, int strategy,
    const OsiSolverInterface* const solver, const OsiCuts* cuts, BBInfo& info,
    const double best_bound, const bool addAsLazy = false);

/*
void doBranchAndBoundNoCuts(const VPCParametersNamespace::VPCParameters& params, const OsiSolverInterface* const solver, BBInfo& info);
void doBranchAndBoundYesCuts(const VPCParametersNamespace::VPCParameters& params, const OsiSolverInterface* const solver,
    BBInfo& info, const OsiCuts& structCuts, const bool doCutSelection,
    const int numCutsToAddPerRound, const int maxRounds,
    const std::string logstring);
*/
#endif /* USE_CBC */
