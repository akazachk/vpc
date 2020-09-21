/**
 * @file SolverInterface.hpp
 * @author A. M. Kazachkov
 * @brief SolverInterface typedef
 */
#pragma once

/// Define SolverInterface that we can later change in a way that will be used in all VPC files
#ifdef USE_CLP_SOLVER
  #include <OsiClpSolverInterface.hpp>
  using SolverInterface = OsiClpSolverInterface; ///< SolverInterface is set to OsiClpSolverInterface
#elif USE_CPLEX_SOLVER
  #include <OsiCpxSolverInterface.hpp>
  #include <ilcplex/cplex.h>
  using SolverInterface = OsiCpxSolverInterface; ///< SolverInterface is set to OsiCpxSolverInterface (this will not work at the moment)
#else
  #include <OsiSolverInterface.hpp>
  using SolverInterface = OsiSolverInterface; ///< SolverInterface is set to OsiSolverInterface (this will not work)
#endif
