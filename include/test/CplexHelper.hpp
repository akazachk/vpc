/**
 * @file CplexHelper.hpp
 * @author A. M. Kazachkov
 * @date 2020-May-13
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

#ifdef USE_CPLEX
// Callable
void presolveModelWithCplexCallable(const VPCParametersNamespace::VPCParameters& params, int strategy,
    const char* f_name, double& presolved_opt, std::string& presolved_name, 
    const double best_bound);
void presolveModelWithCplexCallable(const VPCParametersNamespace::VPCParameters& params, int strategy,
    const OsiSolverInterface* const solver, double& presolved_opt, std::string& presolved_name);
void doBranchAndBoundWithCplexCallable(const VPCParametersNamespace::VPCParameters& params, int strategy,
    const char* f_name, BBInfo& info,
    const double best_bound = std::numeric_limits<double>::max(),
    std::vector<double>* const solution = NULL);
void doBranchAndBoundWithCplexCallable(const VPCParametersNamespace::VPCParameters& params, int strategy,
    const OsiSolverInterface* const solver, BBInfo& info,
    const double best_bound = std::numeric_limits<double>::max(),
    std::vector<double>* const solution = NULL);
void doBranchAndBoundWithUserCutsCplexCallable(const VPCParametersNamespace::VPCParameters& params, int strategy,
    const char* f_name, const OsiCuts* cuts, BBInfo& info,
    const double best_bound, const bool addAsLazy = false);
void doBranchAndBoundWithUserCutsCplexCallable(const VPCParametersNamespace::VPCParameters& params, int strategy,
    const OsiSolverInterface* const solver, const OsiCuts* cuts, BBInfo& info,
    const double best_bound, const bool addAsLazy = false);

#ifdef USE_CPLEX_CONCERT
// Concert
void doBranchAndBoundWithCplexConcert(const VPCParametersNamespace::VPCParameters& params, int strategy,
    const char* f_name, BBInfo& info,
    const double best_bound = std::numeric_limits<double>::max(),
    std::vector<double>* const solution = NULL);
void doBranchAndBoundWithCplexConcert(const VPCParametersNamespace::VPCParameters& params, int strategy,
    const OsiSolverInterface* const solver, BBInfo& info,
    const double best_bound = std::numeric_limits<double>::max(),
    std::vector<double>* const solution = NULL);
void doBranchAndBoundWithUserCutsCplexConcert(const VPCParametersNamespace::VPCParameters& params, int strategy,
    const char* f_name, const OsiCuts* cuts, BBInfo& info,
    const double best_bound, const bool addAsLazy = false);
void doBranchAndBoundWithUserCutsCplexConcert(const VPCParametersNamespace::VPCParameters& params, int strategy,
    const OsiSolverInterface* const solver, const OsiCuts* cuts, BBInfo& info,
    const double best_bound, const bool addAsLazy = false);
#endif /* USE_CPLEX_CONCERT */
#endif /* USE_CPLEX */
