/**
 * @file utility.cpp
 * @author A. M. Kazachkov
 * @date 2018-Dec-24
 */
#include "utility.hpp"
#include "SolverHelper.hpp"

#include <algorithm> // tolower
#include <fstream>
#include <sstream>
#include <sys/stat.h> // for fexists
//#include <cstdio> // for tmpnam

#include <CoinPackedVectorBase.hpp>
#include <CoinPackedVector.hpp>
#include <CoinPackedMatrix.hpp>
#include <OsiSolverInterface.hpp>

/**
 * @details Creates temporary file (in /tmp) so that it can be accessed later.
 * It does not delete the file.
 */
void createTmpFilename(std::string& f_name,
    const std::string add_ext, const std::string tmp_folder) {
//  if (f_name.empty()) {
    // Generate temporary file name
    std::string str_name = tmp_folder.empty() ? "/tmp/tmpvpcXXXXXX" : tmp_folder + "/tmpvpcXXXXXX";
    char template_name[str_name.size() + 1];
    strcpy(template_name, str_name.c_str());

    mkstemp(template_name);
    f_name = template_name;
    if (f_name.empty()) {
      error_msg(errorstring, "Could not generate temp file.\n");
      throw errorstring;
    }
    f_name += add_ext;
//  }
} /* createTmpFilename */

/// \return 0 if successful, 1 if error
int parseFilename(
    /// [out] parent directory
    std::string& dir,
    /// [out] stub of filename without extension
    std::string& instname,
    /// [out] only the extension (after stripping gz/bz2)
    std::string& in_file_ext,
    /// [in] filename to be parsed
    const std::string& fullfilename,
    /// [in] where to write errors
    FILE* logfile) {
  // Get file name stub
  size_t found_dot = fullfilename.find_last_of(".");
  std::string filename = fullfilename.substr(0, found_dot);

  // Put string after last '.' into string in_file_ext
  if (found_dot >= fullfilename.length()) {
    warning_msg(warnstring, "Cannot find the file extension (no '.' in input file name: \"%s\").\n", fullfilename.c_str());
//    writeErrorToLog(errorstring, logfile);
//    exit(1);
    return 1;
  }

  // Check if archived file
  in_file_ext = fullfilename.substr(found_dot + 1);
  if (in_file_ext.compare("gz") == 0 || in_file_ext.compare("bz2") == 0) {
    unsigned found_dot_tmp = filename.find_last_of('.');

    // Put string after last '.' into string in_file_ext
    if (found_dot_tmp < filename.length()) {
      in_file_ext = filename.substr(found_dot_tmp + 1);
      filename = filename.substr(0, found_dot_tmp);
    } else {
      warning_msg(warnstring,
          "Other than gz or bz2, cannot find the file extension (no '.' in input file name: \"%s\").\n",
          fullfilename.c_str());
//      writeErrorToLog(errorstring, logfile);
//      exit(1);
      return 1;
    }
  }

  //  const std::string filestub = (slashindex != std::string::npos) ? fullfilename.substr(slashindex+1) : fullfilename;
  size_t slashindex = filename.find_last_of("/\\");
  dir = (slashindex != std::string::npos) ? filename.substr(0,slashindex) : ".";
  instname = (slashindex != std::string::npos) ? filename.substr(slashindex+1) : filename;
  return 0;
} /* parseFilename (name and logfile given) */

/** @details We assume file is comma separated */
double getObjValueFromFile(std::string opt_filename, std::string fullfilename, FILE* logfile) {
  std::string dir, instname, in_file_ext;
  parseFilename(dir, instname, in_file_ext, fullfilename, logfile);

  if (opt_filename.empty()) {
    return std::numeric_limits<double>::lowest();
  }

  std::ifstream infile(opt_filename.c_str());
  if (infile.is_open()) {
    std::string line;
    while (std::getline(infile, line)) {
      std::istringstream iss(line);
      if (line.empty()) {
        continue;
      }
      std::string curr_instname;
      if (!(std::getline(iss, curr_instname, ','))) {
        warning_msg(warnstring,
            "Could not read instance name. String is %s.\n",
            line.c_str());
        continue;
      }
      if (curr_instname == instname) {
        try {
          std::string token;
          if (!(std::getline(iss, token, ','))) {
            throw;
          }
          return std::stod(token);
        } catch (std::exception& e) {
          warning_msg(warnstring,
              "Could not read optimal value. String is %s.\n",
              line.c_str());
          continue;
        }
      }
    }
    infile.close();
  } else {
    // If we were not able to open the file, throw an error
    error_msg(errorstring, "Not able to open obj file %s.\n", opt_filename.c_str());
    writeErrorToLog(errorstring, logfile);
    exit(1);
  }

  return std::numeric_limits<double>::lowest();
} /* getObjValueFromFile */

void getSolFromFile(
    ///> [in] File with lines "varname value" (space-separated) and comments starting with # or *.
    const char* filename,
    ///> [out] Solution is stored here; space is allocated based on how many variables are listed in \p filename. Needs to be N to match the variables in the linearized model.
    std::vector<double>& sol) {
  if (!filename) {
    return;
  }

  std::ifstream infile(filename);
  if (infile.is_open()) {
    // Reset sol
    sol.clear();

    std::string line;
    while (std::getline(infile, line)) {
      std::istringstream iss(line);
      if (line.empty() || line[0] == '#' || line[0] == '*') {
        continue;
      }
      std::string var_name;
      if (!(std::getline(iss, var_name, ' '))) {
        warning_msg(warnstring,
            "Could not read variable name. String is %s.\n",
            line.c_str());
        continue;
      }
      try {
        std::string token;
        if (!(std::getline(iss, token, ' '))) {
          continue; // unable to find value on this line
        }
        if (token.empty() || token == " ") {
          sol.push_back(0);
        }
        const double val = std::stod(token);
        sol.push_back(val);
      } catch (std::exception& e) {
        warning_msg(warnstring,
            "Could not read value. String is %s.\n",
            line.c_str());
        continue;
      }
    }
    infile.close();
  } else {
    // If we were not able to open the file, throw an error
    error_msg(errorstring, "Not able to open solution file %s.\n", filename);
    exit(1);
  }

  return;
} /* getSolFromFile */

/**
 * @brief Check if file exists
 */
bool fexists(const char* filename) {
  struct stat buffer;
  return (stat (filename, &buffer) == 0);
} /* fexists */

/**
 * @brief Parses int from string using strtol
 */
bool parseInt(const char *str, int& val) {
  long tmpval;
  bool rc = parseLong(str, tmpval);
  val = static_cast<int>(tmpval);
  return rc;
} /* parseInt */

/**
 * @brief Parses long int from string using strtol
 */
bool parseLong(const char *str, long& val) {
  char *temp;
  bool rc = true;
  errno = 0;
  val = strtol(str, &temp, 10);

  if (temp == str || *temp != '\0' ||
      ((val == std::numeric_limits<long>::min() || val == std::numeric_limits<long>::max()) && errno == ERANGE))
    rc = false;

  return rc;
} /* parseLong */

/**
 * @brief Parses double from string using strtod
 */
bool parseDouble(const char *str, double& val) {
  char *temp;
  bool rc = true;
  errno = 0;
  val = strtod(str, &temp);

  if (temp == str || *temp != '\0' ||
      ((val == std::numeric_limits<double>::min() || val == std::numeric_limits<double>::lowest()
        || val == std::numeric_limits<double>::max()) && errno == ERANGE))
    rc = false;

  return rc;
} /* parseDouble */

std::string lowerCaseString(const std::string& tmpData) {
  std::string data = tmpData;
  std::transform(data.begin(), data.end(), data.begin(), ::tolower);
  return data;
}

std::vector<std::string> lowerCaseStringVector(const std::vector<std::string>& strVec) {
  std::vector<std::string> tmp;
  tmp.resize(strVec.size());
  for (int i = 0; i < (int) strVec.size(); i++) {
    tmp[i] = lowerCaseString(strVec[i]);
//    std::transform(tmp[i].begin(), tmp[i].end(), tmp[i].begin(), ::tolower);
  }
  return tmp;
}

std::string upperCaseString(const std::string& tmpData) {
  std::string data = tmpData;
  std::transform(data.begin(), data.end(), data.begin(), ::toupper);
  return data;
}

bool remove_underscore(char c) {
  return (c == '_');
}

std::string upperCaseStringNoUnderscore(const std::string& tmpData) {
  std::string data = upperCaseString(tmpData);
  data.erase(std::remove_if(data.begin(), data.end(), remove_underscore), data.end());
  return data;
}

double dotProductNoCompensation(const CoinPackedVector& vec1, const double* vec2) {
  const int size1 = vec1.getNumElements();
  const int* ind1 = vec1.getIndices();
  const double* el1 = vec1.getElements();
  return dotProductNoCompensation(size1, ind1, el1, vec2);
}

double dotProductNoCompensation(const CoinPackedVector& vec1, const CoinPackedVector& vec2) {
  const int size1 = vec1.getNumElements();
  const int* ind1 = vec1.getIndices();
  const double* el1 = vec1.getElements();
  const int size2 = vec2.getNumElements();
  const int* ind2 = vec2.getIndices();
  const double* el2 = vec2.getElements();
  return dotProductNoCompensation(size1, ind1, el1, size2, ind2, el2);
}

double dotProductNoCompensation(const double* a, const double* b, int dimension) {
  // Accumulate the dot product here
  double accumulator = 0.0;
  // Compensation term for low-order bits lost
  double compensation = 0.0;
  // Temporary number for storing current term
  double currterm = 0.0;
  // Temporary number for storing new accumulator
  double nextacc = 0.0;
  for (int i = 0; i < dimension; ++i) {
    // Compute current term of the dot product adding compensation
    currterm = (a[i] * b[i]) - compensation;
    // This is the value of the sum
    nextacc = accumulator + currterm;
    // Recover what we just lost adding currterm to accumulator
//    compensation = (nextacc - accumulator) - currterm;
    // Now save new value of the accumulator
    accumulator = nextacc;
  }
  return accumulator;
} /* dotProductNoCompensation (dense x dense) */

double dotProductNoCompensation(int sizea, const int* indexa, const double* a,
    const double* b) {
  // Accumulate the dot product here
  double accumulator = 0.0;
  // Compensation term for low-order bits lost
  double compensation = 0.0;
  // Temporary number for storing current term
  double currterm = 0.0;
  // Temporary number for storing new accumulator
  double nextacc = 0.0;
  for (int i = 0; i < sizea; ++i) {
    // Compute current term of the dot product adding compensation
    currterm = (a[i] * b[indexa[i]]) - compensation;
    // This is the value of the sum
    nextacc = accumulator + currterm;
    // Recover what we just lost adding currterm to accumulator
//    compensation = (nextacc - accumulator) - currterm;
    // Now save new value of the accumulator
    accumulator = nextacc;
  }
  return accumulator;
} /* dotProductNoCompensation (sparse x dense) */

double dotProductNoCompensation(int sizea, const int* indexa, const double* a, int sizeb,
    const int* indexb, const double* b) {
  // Accumulate the dot product here
  double accumulator = 0.0;
  // Compensation term for low-order bits lost
  double compensation = 0.0;
  // Temporary number for storing current term
  double currterm = 0.0;
  // Temporary number for storing new accumulator
  double nextacc = 0.0;
  // Current position in vectors a and b
  int posa = 0, posb = 0;
  while (posa < sizea && posb < sizeb) {
    // If it is the same component, compute dot product
    if (indexa[posa] == indexb[posb]) {
      // Compute current term of the dot product adding compensation
      currterm = (a[posa] * b[posb]) - compensation;
      // This is the value of the sum
      nextacc = accumulator + currterm;
      // Recover what we just lost adding currterm to accumulator
//      compensation = (nextacc - accumulator) - currterm;
      // Now save new value of the accumulator
      accumulator = nextacc;
      // Increment both position indices
      ++posa;
      ++posb;
    } else if (indexa[posa] < indexb[posb]) {
      // Increment only smaller position index
      ++posa;
    } else if (indexa[posa] > indexb[posb]) {
      // Increment only smaller position index
      ++posb;
    }
  }
  return accumulator;
} /* dotProductNoCompensation (sparse x sparse) */


// Last edit: 03/27/12
//
// Name:     common_definitions.cpp
// Author:   Giacomo Nannicini
//           Singapore University of Technology and Design, Singapore
//           email: nannicini@sutd.edu.sg
// Date:     04/09/11
//-----------------------------------------------------------------------------
// Copyright (C) 2011, Giacomo Nannicini.  All Rights Reserved.

double dotProduct(const CoinPackedVector& vec1, const double* vec2) {
  const int size1 = vec1.getNumElements();
  const int* ind1 = vec1.getIndices();
  const double* el1 = vec1.getElements();
  return dotProduct(size1, ind1, el1, vec2);
}

double dotProduct(const CoinPackedVector& vec1, const CoinPackedVector& vec2) {
  const int size1 = vec1.getNumElements();
  const int* ind1 = vec1.getIndices();
  const double* el1 = vec1.getElements();
  const int size2 = vec2.getNumElements();
  const int* ind2 = vec2.getIndices();
  const double* el2 = vec2.getElements();
  return dotProduct(size1, ind1, el1, size2, ind2, el2);
}

/**********************************************************/
double dotProduct(const double* a, const double* b, int dimension) {
  // Accumulate the dot product here
  double accumulator = 0.0;
  // Compensation term for low-order bits lost
  double compensation = 0.0;
  // Temporary number for storing current term
  double currterm = 0.0;
  // Temporary number for storing new accumulator
  double nextacc = 0.0;
  for (int i = 0; i < dimension; ++i) {
    // Compute current term of the dot product adding compensation
    currterm = (a[i] * b[i]) - compensation;
    // This is the value of the sum
    nextacc = accumulator + currterm;
    // Recover what we just lost adding currterm to accumulator
    compensation = (nextacc - accumulator) - currterm;
    // Now save new value of the accumulator
    accumulator = nextacc;
  }
  return accumulator;
} /* dotProduct (dense x dense) */

/**********************************************************/
double dotProduct(int sizea, const int* indexa, const double* a,
    const double* b) {
  // Accumulate the dot product here
  double accumulator = 0.0;
  // Compensation term for low-order bits lost
  double compensation = 0.0;
  // Temporary number for storing current term
  double currterm = 0.0;
  // Temporary number for storing new accumulator
  double nextacc = 0.0;
  for (int i = 0; i < sizea; ++i) {
    // Compute current term of the dot product adding compensation
    currterm = (a[i] * b[indexa[i]]) - compensation;
    // This is the value of the sum
    nextacc = accumulator + currterm;
    // Recover what we just lost adding currterm to accumulator
    compensation = (nextacc - accumulator) - currterm;
    // Now save new value of the accumulator
    accumulator = nextacc;
  }
  return accumulator;
} /* dotProduct (sparse x dense) */

/**********************************************************/
double dotProduct(int sizea, const int* indexa, const double* a, int sizeb,
    const int* indexb, const double* b) {
  // Accumulate the dot product here
  double accumulator = 0.0;
  // Compensation term for low-order bits lost
  double compensation = 0.0;
  // Temporary number for storing current term
  double currterm = 0.0;
  // Temporary number for storing new accumulator
  double nextacc = 0.0;
  // Current position in vectors a and b
  int posa = 0, posb = 0;
  while (posa < sizea && posb < sizeb) {
    // If it is the same component, compute dot product
    if (indexa[posa] == indexb[posb]) {
      // Compute current term of the dot product adding compensation
      currterm = (a[posa] * b[posb]) - compensation;
      // This is the value of the sum
      nextacc = accumulator + currterm;
      // Recover what we just lost adding currterm to accumulator
      compensation = (nextacc - accumulator) - currterm;
      // Now save new value of the accumulator
      accumulator = nextacc;
      // Increment both position indices
      ++posa;
      ++posb;
    } else if (indexa[posa] < indexb[posb]) {
      // Increment only smaller position index
      ++posa;
    } else if (indexa[posa] > indexb[posb]) {
      // Increment only smaller position index
      ++posb;
    }
  }
  return accumulator;
} /* dotProduct (sparse x sparse) */

/**
 * Get the two norm of a row, using its packed form
 */
double getRowTwoNorm(const int row, const CoinPackedMatrix* const mat) {
  const int start = mat->getVectorFirst(row);
  const int size = mat->getVectorSize(row);
  double prod = 0.;
  for (int i = start; i < start + size; i++) {
    const double elem = mat->getElements()[i];
    prod += elem * elem;
  }
  prod = std::sqrt(prod);
  return prod;
} /* getRowTwoNorm */

/** Assume sorted vector **/
void packedSortedVectorSum(CoinPackedVector& sum, const double mult1,
    const CoinPackedVectorBase& vec1, const double mult2,
    const CoinPackedVectorBase& vec2, const double eps) {
  if (isZero(mult1, eps) || isZero(mult2, eps)) {
    if (isZero(mult1, eps) && isZero(mult2, eps)) {
      return;
    }
    double mult = isZero(mult1, eps) ? mult2 : mult1;
    if (isZero(mult1, eps)) {
      sum = vec2;
    } else {
      sum = vec1;
    }
    for (int i = 0; i < sum.getNumElements(); i++) {
      sum.setElement(i, mult * sum.getElements()[i]);
    }
    return;
  }

  const int numVec1Elems = vec1.getNumElements();
  const int numVec2Elems = vec2.getNumElements();
  const int* vec1Index = vec1.getIndices();
  const int* vec2Index = vec2.getIndices();
  const double* vec1Val = vec1.getElements();
  const double* vec2Val = vec2.getElements();

  std::vector<int> sumIndex;
  std::vector<double> sumVal;
  sumIndex.reserve(numVec1Elems + numVec2Elems);
  sumVal.reserve(numVec1Elems + numVec2Elems);

  int vec2_el = 0;
  for (int el = 0; el < numVec1Elems; el++) {
    const int vec1_ind = vec1Index[el];

    double new_val = mult1 * vec1Val[el]; 
    while (vec2_el < numVec2Elems) {
      const int vec2_ind = vec2Index[vec2_el];
      if (vec2_ind < vec1_ind) {
        if (!isZero(vec2Val[vec2_el], eps)) {
          sumIndex.push_back(vec2_ind);
          sumVal.push_back(mult2 * vec2Val[vec2_el]);
        }
        vec2_el++;
      } else if (vec2_ind == vec1_ind) {
        new_val += mult2 * vec2Val[vec2_el];
        vec2_el++;
        break;
      } else {
        break;
      }
    }

    if (!isZero(new_val, eps)) {
      sumIndex.push_back(vec1_ind);
      sumVal.push_back(new_val);
    }
  }

  while (vec2_el < numVec2Elems) {
    if (!isZero(vec2Val[vec2_el], eps)) {
      sumIndex.push_back(vec2Index[vec2_el]);
      sumVal.push_back(mult2 * vec2Val[vec2_el]);
    }
    vec2_el++;
  }

  sum.setVector(sumIndex.size(), sumIndex.data(), sumVal.data(), false);
} /* packedSortedVectorSum */

/** @details determine if the variable bounds of solver1 contain the variable
 * bounds of solver2 */
bool variableBoundsContained(const OsiSolverInterface* const solver1,
                             const OsiSolverInterface* const solver2) {

  for (int i = 0; i < solver1->getNumCols(); i++) {
    if (solver1->getColLower()[i] > solver2->getColLower()[i] + 1e-7) {
      return false;
    }
    if (solver1->getColUpper()[i] < solver2->getColUpper()[i] - 1e-7) {
      return false;
    }
  }
  return true;
} /* variableBoundsContain */

/// @brief create a mutable solver interface
std::shared_ptr<SolverInterface> getSolver(const OsiSolverInterface* const si, FILE* logfile) {
  // setup the solver for getting term bases
  std::shared_ptr<SolverInterface> solver;
  try {
    solver = std::make_shared<SolverInterface>(*dynamic_cast<SolverInterface*>(si->clone()));
  } catch (std::exception& e) {
    error_msg(errorstring, "Unable to clone solver into desired SolverInterface.\n");
    writeErrorToLog(errorstring, logfile);
    exit(1);
  }
#ifdef USE_CLP
  try {
    setupClpForStrongBranching(dynamic_cast<OsiClpSolverInterface*>(solver.get()));
  } catch (std::exception& e) {
    // It's okay, we can continue
  }
#endif
  // make sure solver is optimal
  if (!solver->isProvenOptimal()) {
    solver->initialSolve();
    if (!solver->isProvenOptimal()) {
      error_msg(errorstring, "Solver must be optimal to create disjunction.\n");
      writeErrorToLog(errorstring, logfile);
      exit(1);
    }
  }
  solver->enableFactorization();
  solver->markHotStart();
  return solver;
}

/// @brief Find the branching decisions leading to the child that are not present in the parent (i.e. fixed by strong branching)
std::vector<int> findIndicesOfDifference(std::vector<int> child_var, std::vector<int> parent_var,
                                         std::vector<int> child_bound, std::vector<int> parent_bound,
                                         std::vector<double> child_value, std::vector<double> parent_value) {

  // vector to hold which indices of child branching decisions are not present in parent
  std::vector<int> indices;

  // Iterate through branching decisions leading to child node
  for (int i = 0; i < child_var.size(); ++i) {
    bool found = false;

    // iterate over branching decisions leading to parent node
    for (int j = 0; j < parent_var.size(); ++j) {

      // if child branching decision is found, move on to check the next decision
      if (child_var[i] == parent_var[j] && child_bound[i] == parent_bound[j] && child_value[i] == parent_value[j]) {
        found = true;
        break;
      }
    }
    if (!found) {
      // Element from child_var not found in parent_var
      indices.push_back(i);
    }
  }

  return indices;
}

/// @brief Return message if condition is not true
void verify(bool condition, const std::string& msg) {
  if (!condition) {
    std::cerr << msg << std::endl;
    assert(false);  // cause failure with a stack trace
  }
}

/**
 * @details Set the <bound> of <col> to <val> in <solver>.
 *
 * @param solver The solver to modify.
 * @param col The column to bound.
 * @param bound The bound to set. 0 is lower bound, 1 is upper bound.
 * @param val The value to set the bound to.
 *
 * @return void
 */
void addVarBound(OsiSolverInterface* solver, const int col, const int bound, const double val){
  checkColumnAndBound(solver, col, bound);
  if (bound <= 0) {
    solver->setColLower(col, val);
  } else {
    solver->setColUpper(col, val);
  }
}

/**
 * @details Set the <bound> of <col> to <val> in <solver>. Additionally, append
 * the variable, bound, and value to the end of the respective vectors.
 *
 * @param solver The solver to modify.
 * @param col The column to bound.
 * @param bound The bound to set. 0 is lower bound, 1 is upper bound.
 * @param val The value to set the bound to.
 *
 * @return void
 */
void addVarBound(OsiSolverInterface* solver, const int col, const int bound,
                 const double val, std::vector<int>& fixed_var,
                 std::vector<int>& fixed_bound, std::vector<double>& fixed_value){
  addVarBound(solver, col, bound, val);
  fixed_var.push_back(col);
  fixed_bound.push_back(bound);
  fixed_value.push_back(val);
}

/**
 * @details Check that the bounds encoded in fixed_var, fixed_bound, and fixed_value
 * exist in solver if isTrue is true. Otherwise, they should not exist
 *
 * @param solver The solver to check.
 * @param fixed_var The variables to check.
 * @param fixed_bound The bounds to check.
 * @param fixed_value The values to check.
 * @param isTrue Whether the bounds should be exist in the solver already or not
 *
 */
void checkBounds(const OsiSolverInterface* const solver, const std::vector<int>& fixed_var,
                 const std::vector<int>& fixed_bound, const std::vector<double>& fixed_value,
                 bool isTrue){

  verify((fixed_var.size() == fixed_bound.size()) && (fixed_bound.size() == fixed_value.size()),
         "fixed_var, fixed_bound, and fixed_value must be the same size.");

  for (int idx = 0; idx < fixed_var.size(); ++idx) {
    const int col = fixed_var[idx];
    const int bound = fixed_bound[idx];
    const double val = fixed_value[idx];

    checkColumnAndBound(solver, col, bound);
    // checking for satisfaction and not equality because general integers may be branched on multiple times
    if (bound <= 0) {
      verify((solver->getColLower()[col] >= val) == isTrue,
             "The branching decisions in fixed_var, fixed_bound, and fixed_value "
             "are inconsistent with the expectations on the solver.");
    } else {
      verify((solver->getColUpper()[col] <= val) == isTrue,
             "The branching decisions in fixed_var, fixed_bound, and fixed_value "
             "are inconsistent with the expectations on the solver.");
    }
  }
}

/**
 * @details Check that the bound is either lower (0) or upper (1) and that the
 * column is required to be integer
 *
 * @param solver The solver to check
 * @param col The column to check
 * @param bound The bound to check (0 for lower (up branch) and 1 for upper (down branch))
 */
void checkColumnAndBound(const OsiSolverInterface* const solver, const int col,
                         const int bound){
  verify(bound == 0 || bound == 1, "bound must be 0 or 1.");
  verify(0 <= col && col < solver->getNumCols(), "col must be a valid column index.");
  verify(solver->isInteger(col), "col must be an integer variable.");
}

/**
 * @details sort the encoding of branching decisions by the variable index
 *
 * @param vars the indices of variables branched on
 * @param bounds the directions of branching that occurred
 * @param vals the values that variables were branched on
 * @return void
 */
void sortBranchingDecisions(std::vector<int>& vars, std::vector<int>& bounds,
                            std::vector<double>& vals){

  // Create a vector of tuples to store the original values and their indices
  std::vector<std::tuple<int, int, double>> combined;

  // Fill the combined vector with tuples of vars and their corresponding bounds and vals
  for (size_t i = 0; i < vars.size(); i++) {
    combined.emplace_back(vars[i], bounds[i], vals[i]);
  }

  // Sort the combined vector by the first element (the vars)
  std::sort(combined.begin(), combined.end(), [](const auto& a, const auto& b) {
    return std::get<0>(a) < std::get<0>(b);
  });

  // Extract the sorted values into the original vectors
  for (size_t i = 0; i < vars.size(); ++i) {
    vars[i] = std::get<0>(combined[i]);
    bounds[i] = std::get<1>(combined[i]);
    vals[i] = std::get<2>(combined[i]);
  }
}

/**
 * @details Check if the objective function is minimization. If not, negate it.
 *
 * @param solver the solver to check
 */
void ensureMinimizationObjective(SolverInterface* solver){
  if (solver->getObjSense() < 1e-3) {
    printf(
        "\n## Detected maximization problem. Negating objective function to make it minimization. ##\n");
    solver->setObjSense(1.0);
    const double* obj = solver->getObjCoefficients();
    for (int col = 0; col < solver->getNumCols(); col++) {
      solver->setObjCoeff(col, -1. * obj[col]);
    }
    double objOffset = 0.;
    solver->getDblParam(OsiDblParam::OsiObjOffset, objOffset);
    if (objOffset != 0.) {
      solver->setDblParam(OsiDblParam::OsiObjOffset, -1. * objOffset);
    }
  }
}

/** check if sol is feasible for solver */
bool isFeasible(
    /// [in] problem
    const OsiSolverInterface& solver,
    /// [in] solution
    const std::vector<double>& sol) {

  // get bounds and constraint coefficients
  const double* rowLower = solver.getRowLower();
  const double* rowUpper = solver.getRowUpper();
  const double* colLower = solver.getColLower();
  const double* colUpper = solver.getColUpper();
  const CoinPackedMatrix* mat = solver.getMatrixByRow();
  const int* cols = mat->getIndices();
  const double* elem = mat->getElements();

  // check dimensions match
  verify(sol.size() == solver.getNumCols(), "solution has wrong dimension");

  // make sure constraints are valid
  for (int row = 0; row < solver.getNumRows(); row++) {
    double sum = 0.0;
    const int start = mat->getVectorFirst(row);
    const int end = mat->getVectorLast(row);
    for (int ind = start; ind < end; ind++) {
      const double j = cols[ind];
      const double el = elem[ind];
      sum += el * sol[j];
    }
    if (lessThanVal(sum, rowLower[row]) || greaterThanVal(sum, rowUpper[row])) {
      return false;
    }
  }

  // make sure variables are valid
  for (int col = 0; col < solver.getNumCols(); col++) {
    if (lessThanVal(sol[col], colLower[col]) || greaterThanVal(sol[col], colUpper[col])) {
      return false;
    }
    if (solver.isInteger(col) && !isInteger(sol[col])) {
      return false;
    }
  }
  return true;
} /* isFeasible */

/** check if a value is an integer */
bool isInteger(double val){
  return isZero(min(val - std::floor(val), std::ceil(val) - val));
}

/** take a min of two values */
double min(double a, double b){
  if (a < b){
    return a;
  } else {
    return b;
  }
}

/** take a max of two values */
double max(double a, double b){
  if (a > b){
    return a;
  } else {
    return b;
  }
}