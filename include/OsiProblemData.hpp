#pragma once

/* Data structure to pass to/from Disjunctive cut generator the problem over which to generate the cuts.
   This is to conform to OsiSolverInterface way of passing a problem */
enum OsiVarType {
  CONTINUOUS = 0,
  BINARY = 1,
  INTEGER = 2,
  SEMICONT = 3,
  SEMIINT = 4,
  NUM_VAR_TYPES = 5
};
const char OsiVarTypeChar[OsiVarType::NUM_VAR_TYPES] {
  'C', 'B', 'I', 'S', 'N'
};
struct OsiProblemData
{
   int numcols;
   int numrows;
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
#include <ilclpex/cplex.h>
int ConvertCPX2Data(CPXENVptr env, CPXLPptr lp, OsiProblemData* pdata);
#endif
void MemSetProbData( OsiProblemData* pdata, int ncols, int nrows, int nz );
void MemSetProbData( OsiProblemData* pdata, int ncols, int nrows, int nz );
