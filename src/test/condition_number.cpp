/**
 * @file condition_number.cpp
 * @author A. M. Kazachkov, based on code provided by S. Nadarajah
 * @date 2023-07-02
 */
#include "condition_number.hpp"

// COIN-OR files
#include <OsiSolverInterface.hpp>

// Singular Value Decomposition (from Lapack);
// needed to compute condition number in norm 2
extern "C" int dgesvd_(char *jobu, char *jobvt, int *m, int *n, double *a,
    int *lda, double *s, double *u, int *ldu, double *vt, int *ldvt,
    double *work, int *lwork, int *info);

double compute_condition_number_norm2(const OsiSolverInterface* const solver) {
  double condition = 0.0;
  if (!solver->basisIsAvailable()) {
    return condition;
  }
  int m = solver->getNumRows();
  int n = solver->getNumCols();
  int* rstat = new int[m];
  int* cstat = new int[n];
  solver->getBasisStatus(cstat, rstat);
  double* basis = new double[m * m];
  memset(basis, 0, m * m * sizeof(double));
  const CoinPackedMatrix* byCol = solver->getMatrixByCol();
  const int * rowPos = byCol->getIndices();
  const CoinBigIndex * columnStart = byCol->getVectorStarts();
  const int * columnLength = byCol->getVectorLengths();
  const double * columnElements = byCol->getElements();
  int bindex = 0;
  // fill the structural columns of the basis
  for (int i = 0; i < n; ++i) {
    if (cstat[i] == 1) {
      for (int j = columnStart[i]; j < columnStart[i] + columnLength[i];
          ++j) {
        basis[bindex * m + rowPos[j]] = columnElements[j];
      }
      bindex++;
    }
  }
  // fill the artificial columns of the basis
  for (int i = 0; i < m; ++i) {
    if (rstat[i] == 1) {
      basis[bindex * m + i] = 1.0;
      bindex++;
    }
  }
  // prepare data for SVD; get the optimal size of the work array
  int info = 0;
  char job = 'N';
  double* sv = new double[m];
  double sizeWork = 0.0;
  // we first call the SVD routine with lwork = -1; this way, the function
  // returns the optimal size of the work array that is needed for the
  // computation (the value is written in sizeWork)
  int lwork = -1;
  dgesvd_(&job, &job, &m, &m, basis, &m, sv, NULL, &m, NULL, &m, &sizeWork,
      &lwork, &info);
  if (info != 0) {
    printf(
        "### WARNING: can not obtain optimal size of work array for SVD\n");
    return condition;
  }
  // now we know the optimal size, so we allocate work array and do SVD
  lwork = (int) sizeWork;
  double* work = new double[lwork];
  dgesvd_(&job, &job, &m, &m, basis, &m, sv, NULL, &m, NULL, &m, work, &lwork,
      &info);
  if (info != 0) {
    printf("### WARNING: can not compute condition number\n");
    return condition;
  }
  condition = sv[0] / sv[m - 1];
  delete[] rstat;
  delete[] cstat;
  delete[] basis;
  delete[] sv;
  delete[] work;
  return condition;
} /* compute_condition_number_norm2 */

double compute_condition_number_norm1(const OsiSolverInterface* const ref_solver) {
  OsiSolverInterface* solver = ref_solver->clone();
  double condition = 0.0;
  if (!solver->basisIsAvailable()) {
    return condition;
  }
  int m = solver->getNumRows();
  int n = solver->getNumCols();
  int* rstat = new int[m];
  int* cstat = new int[n];
  solver->getBasisStatus(cstat, rstat);
  const CoinPackedMatrix* byCol = solver->getMatrixByCol();
  const CoinBigIndex * columnStart = byCol->getVectorStarts();
  const int * columnLength = byCol->getVectorLengths();
  const double * columnElements = byCol->getElements();
  double norm1B = 0.0;
  double currnorm;
  // norm of structural columns
  for (int i = 0; i < n; ++i) {
    if (cstat[i] == 1) {
      currnorm = 0.0;
      for (int j = columnStart[i]; j < columnStart[i] + columnLength[i];
          ++j) {
        currnorm += std::abs(columnElements[j]);
      }
      if (currnorm > norm1B) {
        norm1B = currnorm;
      }
    }
  }
  // norm of artificial columns
  if (norm1B < 1.0) {
    norm1B = 1.0;
  }
  solver->enableSimplexInterface(true);
  solver->enableFactorization();
  double norm1Binv = 0.0;
  double* binvcol = new double[m];
  for (int i = 0; i < m; ++i) {
    solver->getBInvCol(i, binvcol);
    currnorm = 0.0;
    for (int j = 0; j < m; ++j) {
      currnorm += std::abs(binvcol[j]);
    }
    if (currnorm > norm1Binv) {
      norm1Binv = currnorm;
    }
  }
  solver->disableFactorization();
  solver->disableSimplexInterface();
  condition = norm1B * norm1Binv;
  delete[] rstat;
  delete[] cstat;
  delete[] binvcol;
  return condition;
} /* compute_condition_number_norm1 */

// compute the condition number in the given norm (1 or 2)
double compute_condition_number(const OsiSolverInterface* const solver, const int norm) {
  if (norm < 1 || norm > 2) {
    printf("### WARNING: can not compute condition number in norm %d\n",
        norm);
    return 0.0;
  }
  if (norm == 1) {
    return compute_condition_number_norm1(solver);
  }
  return compute_condition_number_norm2(solver);
} /* compute_condition_number */
