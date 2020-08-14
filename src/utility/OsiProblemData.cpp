#include "OsiProblemData.hpp"
#include <cstdlib> // malloc

// Initializes and allocates pdata. pdata is assumed to not have been initialized yet

#ifdef USE_CPLEX
int ConvertCPX2Data(CPXENVptr env, CPXLPptr lp, OsiProblemData* pdata)
{
   int rval = 0;

   int nzcnt;
   int* cmatbeg = 0;
   int* cmatind = 0;
   double* cmatval = 0;
   double* obj = 0;
   double* rhs = 0;
   double* lb = 0;
   double* ub = 0;
   char* sense = 0;
   char* ctype = 0;
   int cmatspace = 0;
   int surplus = 0;
   int ncols, nrows;
   int objsense;
   double objoffset = 0.;

   ncols = CPXgetnumcols( env, lp );
   nrows = CPXgetnumrows( env, lp );
   objsense = CPXgetobjsen( env, lp );
   CPXgetobjoffset( env, lp, &objoffset );

   cmatbeg = (int*) malloc( sizeof(int) * ncols );
   obj = (double*) malloc( sizeof(double) * ncols );
   lb = (double*) malloc( sizeof(double) * ncols );
   ub = (double*) malloc( sizeof(double) * ncols );
   rhs = (double*) malloc( sizeof(double) * nrows );
   sense = (char*) malloc( sizeof(char) * nrows );
   ctype = (char*) malloc( sizeof(char) * ncols );


   rval = CPXgetrhs( env, lp, rhs, 0, nrows-1);

   rval = CPXgetobj( env, lp, obj, 0, ncols-1);

   rval = CPXgetctype( env, lp, ctype, 0, ncols-1);
   rval = CPXgetlb( env, lp, lb, 0, ncols-1);
   rval = CPXgetub( env, lp, ub, 0, ncols-1);

   rval = CPXgetsense( env, lp, sense, 0, nrows-1);


   rval = CPXgetcols( env, lp, &nzcnt, cmatbeg, cmatind, cmatval, cmatspace, &surplus, 0, ncols-1);
   // first call should return an error
   if( rval == CPXERR_NEGATIVE_SURPLUS )
   {
       cmatspace = -surplus;
       cmatind = (int*) malloc( sizeof(int) * cmatspace );
       cmatval = (double*) malloc( sizeof(double) * cmatspace );
       rval = CPXgetcols( env, lp, &nzcnt, cmatbeg, cmatind, cmatval, cmatspace, &surplus, 0, ncols-1);
       if( rval )
       {
          fprintf(stderr," EXPECTED RVAL = 0 in line %d of file %s\n", __LINE__, __FILE__ );
       }
       else 
       {
         MemSetProbData( pdata, ncols, nrows, nzcnt, objsense, objoffset );

         // Now put all cols in the corresponding data structures
         int currOSIidx = 0;
         for( int j = 0; j < ncols; j++ )
         {
           int begcpxidx, endcpxidx;
           begcpxidx = cmatbeg[j];
           if( j < ncols - 1 )
             endcpxidx = cmatbeg[ j+1 ];
           else
             endcpxidx = nzcnt;

           // Watch for infinity conversion
           pdata->collb[j] = lb[j];
           pdata->colub[j] = ub[j];

           pdata->obj[j] = obj[j];
           pdata->start[j] = currOSIidx;
           pdata->vartype[j] = ctype[j];
           for( int spj = begcpxidx; spj < endcpxidx; spj++ )
           {
             pdata->index[ currOSIidx ] = cmatind[spj];
             pdata->value[ currOSIidx ] = cmatval[spj];

             currOSIidx++;
           }
         }
         pdata->start[ncols] = currOSIidx;
         if( currOSIidx != nzcnt )
         {
           rval = 1;
           fprintf(stderr," Problem in line %d of file %s\n", __LINE__, __FILE__ );
         }
         else
         {
           for( int i = 0; i < nrows; i++)
           {
             if( sense[i] == 'G')
             {
               pdata->rowub[i] = CPX_INFBOUND;
               pdata->rowlb[i] = rhs[i];
             }
             else if( sense[i] == 'L')
             {
               pdata->rowub[i] = rhs[i];
               pdata->rowlb[i] = -CPX_INFBOUND;
             }
             else if( sense[i] == 'E')
             {
               pdata->rowub[i] = rhs[i];
               pdata->rowlb[i] = rhs[i];
             }
             else
             {
               rval = 1;
               fprintf(stderr," Problem in line %d of file %s\n", __LINE__, __FILE__ );
             }
           }
         }
       }
   } /* if( rval == CPXERR_NEGATIVE_SURPLUS ) */
   else
   {
       fprintf(stderr," EXPECTED RVAL = %d (actual = %d) in line %d of file %s\n", CPXERR_NEGATIVE_SURPLUS, rval, __LINE__, __FILE__ );
   }

   if( cmatbeg )
     free( cmatbeg );
   if( cmatind )
     free( cmatind );
   if( cmatval )
     free( cmatval );
   if( obj )
     free (obj);
   if( rhs )
     free (rhs);
   if( sense )
     free (sense);
   if( ctype )
     free (ctype);
   if( lb )
     free (lb);
   if( ub )
     free (ub);
   return rval;
} /* ConvertCPX2Data */
#endif /* USE_CPLEX */

void MemSetProbData( OsiProblemData* pdata, int ncols, int nrows, int nz, int objsense, double objoffset )
{
  pdata->numcols = ncols;
  pdata->numrows = nrows;
  pdata->objsense = objsense;
  pdata->objoffset = objoffset;
  pdata->start = (int*) malloc( sizeof(int) * (ncols+1) );
  pdata->index = (int*) malloc( sizeof(int) * nz );
  pdata->value = (double*) malloc( sizeof(double) * nz );
  pdata->collb = (double*) malloc( sizeof(double) * ncols );
  pdata->colub = (double*) malloc( sizeof(double) * ncols );
  pdata->obj = (double*) malloc( sizeof(double) * ncols );
  pdata->rowlb = (double*) malloc( sizeof(double) * nrows );
  pdata->rowub = (double*) malloc( sizeof(double) * nrows );
  pdata->vartype = (char*) malloc( sizeof(char)* ncols );
} /* MemSetProbData */

void FreeProbData( OsiProblemData* pdata )
{
  free( pdata->start );
  free( pdata->index );
  free( pdata->value );
  free( pdata->collb );
  free( pdata->colub );
  free( pdata->obj );
  free( pdata->rowlb );
  free( pdata->rowub );
  free( pdata->vartype );
} /* FreeProbData */
