/**
 * @file GurobiHelper.hpp
 * @author A. M. Kazachkov
 * @date 2019-Feb-28
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

#ifdef USE_GUROBI
void presolveModelWithGurobi(const VPCParametersNamespace::VPCParameters& params, int strategy,
    const char* f_name, double& presolved_lp_opt, std::string& presolved_name,
    const double best_bound);
void presolveModelWithGurobi(const VPCParametersNamespace::VPCParameters& params, int strategy,
    const OsiSolverInterface* const solver, double& presolved_lp_opt, std::string& presolved_name);
void doBranchAndBoundWithGurobi(const VPCParametersNamespace::VPCParameters& params, int strategy, 
    const char* f_name, BBInfo& info,
    const double best_bound = std::numeric_limits<double>::max(),
    std::vector<double>* const solution = NULL);
void doBranchAndBoundWithGurobi(const VPCParametersNamespace::VPCParameters& params, int strategy,
    const OsiSolverInterface* const solver, BBInfo& info,
    const double best_bound = std::numeric_limits<double>::max(),
    std::vector<double>* const solution = NULL);
void doBranchAndBoundWithUserCutsGurobi(const VPCParametersNamespace::VPCParameters& params, int strategy,
    const char* f_name, const OsiCuts* cuts, BBInfo& info,
    const double best_bound, const bool addAsLazy = false);
void doBranchAndBoundWithUserCutsGurobi(const VPCParametersNamespace::VPCParameters& params, int strategy,
    const OsiSolverInterface* const solver, const OsiCuts* cuts, BBInfo& info,
    const double best_bound, const bool addAsLazy = false);
#endif /* USE_GUROBI */
