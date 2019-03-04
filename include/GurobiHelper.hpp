// Name:     GurobiHelper.hpp
// Author:   A. M. Kazachkov
// Date:     2019-Feb-28
//-----------------------------------------------------------------------------
#pragma once

#include "CbcModel.hpp"

struct BBInfo;
struct VPCParameters;

#ifdef USE_GUROBI
// Gurobi stuff
void presolveModelWithGurobi(const VPCParameters& params, int strategy,
    const char* f_name, double& presolved_opt, std::string& presolved_name,
    const double best_bound);
void presolveModelWithGurobi(const VPCParameters& params, int strategy,
    const OsiSolverInterface* const solver, double& presolved_opt,
    std::string& presolved_name);
void doBranchAndBoundWithGurobi(const VPCParameters& params, int strategy, const char* f_name,
    BBInfo& info, const double best_bound, std::vector<double>* const solution = NULL);
void doBranchAndBoundWithGurobi(const VPCParameters& params, int strategy,
    const OsiSolverInterface* const solver, BBInfo& info,
    const double best_bound, std::vector<double>* const solution = NULL);
void doBranchAndBoundWithUserCutsGurobi(const VPCParameters& params, int strategy,
    const char* f_name, const OsiCuts* cuts, BBInfo& info,
    const double best_bound, const bool addAsLazy = false);
void doBranchAndBoundWithUserCutsGurobi(const VPCParameters& params, int strategy,
    const OsiSolverInterface* const solver, const OsiCuts* cuts, BBInfo& info,
    const double best_bound, const bool addAsLazy = false);
#endif /* USE_GUROBI */
