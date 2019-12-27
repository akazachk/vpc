#pragma once

// Define SolverInterface that we can later change in a way that will be used in all VPC files
#ifdef USE_CLP_SOLVER
  #include <OsiClpSolverInterface.hpp>
  using SolverInterface = OsiClpSolverInterface;
#elif USE_CPLEX_SOLVER
  #include <OsiCpxSolverInterface.hpp>
  #include <ilcplex/cplex.h>
  using SolverInterface = OsiCpxSolverInterface;
#else
  #include <OsiSolverInterface.hpp>
  using SolverInterface = OsiSolverInterface;
#endif
