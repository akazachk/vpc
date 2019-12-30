/**
 * @file debug.cpp
 *
 * @author A. M. Kazachkov
 * @date 2018-12-25
 */
#include "debug.hpp"

#include "utility.hpp"

#include <CoinPackedMatrix.hpp>
#include <CoinPackedVectorBase.hpp>
#include <CoinPackedVector.hpp>
#include <CoinShallowPackedVector.hpp>

/**
 * Print sparse vector
 *
 * @param vec           Vector to be printed 
 * @param use_newline   Whether to print the entire vector on one line or not
 */
void printVector(const CoinPackedVectorBase& vec, const bool use_newline) {
  int numElems = vec.getNumElements();
  const int* index = vec.getIndices();
  const double* element = vec.getElements();
  if (use_newline)
    fprintf(stdout, "Num elements is %d.", numElems);
  for (int i = 0; i < numElems; i++) {
    if (use_newline)
      fprintf(stdout, "\n\t");
    else if (i > 0)
      fprintf(stdout, ", ");
    fprintf(stdout, "(%d, %g)", index[i], element[i]);
  }
  fprintf(stdout, "\n");
} /* printVector (CoinPackedVectorBase) */

/** Print several vectors in CoinPackedVectorBase form */
void printVectors(const std::vector<CoinPackedVector>& vecs, const bool use_newline) {
  for (auto& vec : vecs)
    printVector(vec, use_newline);
} /* printVectors */

/**
 * Print dense array
 *
 * @param n     Number of elements in the vector
 * @param vec   Vector to be printed
 */
template <typename T>
void printVector(const int n, const T* vec, const bool use_newline) {
  for (int i = 0; i < n; ++i) {
    if (use_newline)
      fprintf(stdout, "\n");
    else if (i > 0)
      fprintf(stdout, ", ");

    if (std::abs(vec[i]-std::floor(vec[i])) != 0.0)
      fprintf(stdout, "\t(%d, %g)", i, static_cast<double>(vec[i]));
    else
      fprintf(stdout, "\t(%d, %d)", i, static_cast<int>(vec[i]));
  }
  fprintf(stdout, "\n");
} /* printVector (int, T*) */

/**
 * Print vector
 */
template <typename T>
void printVector(const std::vector<T>& vec, const bool use_newline) {
  printVector(vec.size(), vec.data(), use_newline);
} /* printVector (std::vector<T>) */

/**
 * Print matrix (row-wise)
 *
 * @param mx_in   Matrix to be printed
 */
void printMatrix(const CoinPackedMatrix& mx_in) {
  const CoinPackedMatrix* mx; 
  if (mx_in.isColOrdered()) {
    CoinPackedMatrix* mx_rows = new CoinPackedMatrix;
    mx_rows->reverseOrderedCopyOf(mx_in);
    mx = mx_rows;
  } else {
    mx = &mx_in;
  }
  for (int i = 0; i < mx->getNumRows(); i++) {
    const CoinShallowPackedVector& vec = mx->getVector(i);
    fprintf(stdout, "Row %d: ", i);
    printVector(vec, false);
  }
  if (mx_in.isColOrdered()) {
    delete mx;
  }
} /* printMatrix */
