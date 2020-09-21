/**
 * @file OsiProblemData.hpp
 * @author A. M. Kazachkov
 * @brief Useful enums and methods for converting a CPLEX instance into Osi-friendly format
 */
#pragma once

/// Data structure to pass to/from Disjunctive cut generator the problem over which to generate the cuts.
/// This is to conform to OsiSolverInterface way of passing a problem.
enum OsiVarType {
  CONTINUOUS = 0,
  BINARY = 1,
  INTEGER = 2,
  SEMICONT = 3,
  SEMIINT = 4,
  NUM_VAR_TYPES = 5
};
/// Possible variable types are continuous ('C'), binary ('B'), integer ('I'), semicontinuous ('S'), and 'N'
const char OsiVarTypeChar[OsiVarType::NUM_VAR_TYPES] {
  'C', 'B', 'I', 'S', 'N'
};
/// Make data fit for entry into Osi
struct OsiProblemData
{
   int numcols;
   int numrows;
   int objsense;
   double objoffset;
   int *start;
   int *index;
   double *value;
   double *collb;
   double *colub;
   double *obj;
   double *rowlb;
   double *rowub;
   char* vartype;
}; /* OsiProblemData */

#ifdef USE_CPLEX
#include <ilcplex/cplex.h>
/// Read CPLEX instance into OsiProblemData
int ConvertCPX2Data(CPXENVptr env, CPXLPptr lp, OsiProblemData* pdata);
#endif
/// Set OsiProblemData memory
void MemSetProbData( OsiProblemData* pdata, int ncols, int nrows, int nz );
/// Free OsiProblemData memory
void FreeProbData( OsiProblemData* pdata );
