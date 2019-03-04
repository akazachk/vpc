// Name:     CplexHelper.hpp
// Author:   A. M. Kazachkov
// Date:     2019-Mar-01
//-----------------------------------------------------------------------------
#pragma once

#ifdef VPC_USE_CPLEX
// CPLEX stuff
// Callable
void presolveModelWithCplexCallable(const char* f_name, double& presolved_opt, std::string& presolved_name);
void presolveModelWithCplexCallable(const OsiSolverInterface* const solver, double& presolved_opt, std::string& presolved_name);
void doBranchAndBoundWithCplexCallable(const char* f_name, BBInfo& info);
void doBranchAndBoundWithCplexCallable(const OsiSolverInterface* const solver, BBInfo& info);
void doBranchAndBoundWithUserCutsCplexCallable(const char* f_name,
    const OsiCuts* cuts, BBInfo& info, const bool addAsLazy = false);
void doBranchAndBoundWithUserCutsCplexCallable(const OsiSolverInterface* const solver,
    const OsiCuts* cuts, BBInfo& info, const bool addAsLazy = false);

#ifdef VPC_USE_CPLEX_CONCERT
// Concert
void doBranchAndBoundWithCplexConcert(const char* f_name, BBInfo& info);
void doBranchAndBoundWithCplexConcert(const OsiSolverInterface* const solver, BBInfo& info);
void doBranchAndBoundWithUserCutsCplexConcert(const char* f_name,
    const OsiCuts* cuts, BBInfo& info, const bool addAsLazy = false);
void doBranchAndBoundWithUserCutsCplexConcert(const OsiSolverInterface* const solver,
    const OsiCuts* cuts, BBInfo& info, const bool addAsLazy = false);
#endif /* VPC_USE_CPLEX_CONCERT */
#endif /* VPC_USE_CPLEX */
