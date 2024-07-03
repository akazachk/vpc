/**
 * @file utility.cpp
 * @author A. M. Kazachkov
 * @date 2018-Dec-24
 */
#include "utility.hpp"
#include "SolverHelper.hpp"

#include <algorithm> // tolower, min, max
#include <fstream>
#include <sstream>
#include <sys/stat.h> // for fexists
//#include <cstdio> // for tmpnam

#include <CoinPackedVectorBase.hpp>
#include <CoinPackedVector.hpp>
#include <CoinPackedMatrix.hpp>

std::vector<std::string> getStatNames() {
  const int num_stats = static_cast<int>(Stat::num_stats);
  std::vector<std::string> stat_names(num_stats);
  stat_names[static_cast<int>(Stat::total)] = "total";
  stat_names[static_cast<int>(Stat::avg)] = "avg";
  stat_names[static_cast<int>(Stat::min)] = "min";
  stat_names[static_cast<int>(Stat::max)] = "max";
  stat_names[static_cast<int>(Stat::stddev)] = "stddev";
  stat_names[static_cast<int>(Stat::num)] = "num";
  stat_names[static_cast<int>(Stat::num_nonzero)] = "num_nonzero";
  return stat_names;
} /* getStatNames */

void printStats(const StatVector& stats, const bool printNames, const char SEP, FILE* logfile) {
  const int num_stats = static_cast<int>(Stat::num_stats);
  const std::vector<std::string> stat_names = getStatNames();
  if (logfile == NULL) { logfile = stdout; }
  for (int i = 0; i < num_stats; i++) {
    if (printNames) {
      fprintf(logfile, "\t%s: ", stat_names[i].c_str());
    }
    fprintf(logfile, "%s%c", stringValue(stats[i], "%.10g").c_str(), SEP);
  }
  //fprintf(logfile, "\n");
} /* printStats */

template StatVector computeStats<int>(const std::vector<int>& vals);
template StatVector computeStats<double>(const std::vector<double>& vals);

template <typename T>
StatVector computeStats(const std::vector<T>& vals) {
  StatVector stats(static_cast<int>(Stat::num_stats), 0.);
  if (vals.empty()) {
    return stats;
  }
  initializeStats(stats);
  updateAndFinalizeStats(stats, vals, 0);
  return stats;
} /* computeStats */

/// @details Initialize most values to 0, except min to MAX_POSSIBLE_VAL, max to MIN_POSSIBLE_VAL.
void initializeStats(StatVector& stats, const double MIN_POSSIBLE_VAL, const double MAX_POSSIBLE_VAL) {
  stats.assign(static_cast<int>(Stat::num_stats), 0);
  stats[static_cast<int>(Stat::min)] = MAX_POSSIBLE_VAL;
  stats[static_cast<int>(Stat::max)] = MIN_POSSIBLE_VAL;
} /* initializeStats */

template void updateAndFinalizeStats<int>(StatVector& stats, const std::vector<int>& vals, const int prev_size);
template void updateAndFinalizeStats<double>(StatVector& stats, const std::vector<double>& vals, const int prev_size);

/// @details If \p prev_size is positive, then need to adjust previous average/stddev. New total size is \p prev_size + \p vals size.
/// If \p prev_size is negative, use #Stat::num as divisor.
template <typename T>
void updateAndFinalizeStats(StatVector& stats, const std::vector<T>& vals, const int prev_size) {
  if (vals.empty()) { return; }
  const int new_size = prev_size + static_cast<int>(vals.size());

  if (prev_size != 0) {
    // Return previous average to total and stddev to sum of squares
    const int OLD_SIZE = (prev_size > 0) ? prev_size : stats[static_cast<int>(Stat::num)];
    unfinalizeStats(stats, OLD_SIZE);
  }

  for (const auto& val : vals) {
    updateStatsBeforeFinalize(stats, val);
  }
  finalizeStats(stats, new_size);
} /* updateAndFinalizeStats */

template void updateStatsBeforeFinalize<int>(StatVector& stats, const int& val, const int size);
template void updateStatsBeforeFinalize<double>(StatVector& stats, const double& val, const int size);

/// @details Parameter \p size is optional and only applied when > 0; when == 0, then we return immediately, and when < 0, we assume that size will be set later.
template <typename T>
void updateStatsBeforeFinalize(StatVector& stats, const T& val, const int size) {
  const double divisor = (size >= 0) ? static_cast<double>(size) : 1.;
  stats[static_cast<int>(Stat::total)] += val;
  stats[static_cast<int>(Stat::avg)] = stats[static_cast<int>(Stat::total)] / divisor;
  stats[static_cast<int>(Stat::min)] = std::min(stats[static_cast<int>(Stat::min)], static_cast<double>(val));
  stats[static_cast<int>(Stat::max)] = std::max(stats[static_cast<int>(Stat::max)], static_cast<double>(val));

  // For standard deviation, this is not the final value
  // We are computing Var[X] = E[X^2] - E[X]^2
  // Below, we are only computing E[X^2]; we will subtract E[X]^2 later in finalizeStats
  // And afterwards we will take the square root to get the standard deviation
  stats[static_cast<int>(Stat::stddev)] += (double) val * val / divisor;

  stats[static_cast<int>(Stat::num)] += 1;
  if (!isZero(val)) {
    stats[static_cast<int>(Stat::num_nonzero)] += 1;
  }
} /* updateStatsBeforeFinalize */

/// @details If \p size is positive, then assume that Stat::avg = Stat::total so we need to divide by size.
/// Similarly, for stddev if \p size > 0, we need to divide by size (assume Stat::stddev holds sum of squares, not E[X^2] as in case of \p size <= 0).
/// To finish stddev calculation = sqrt(E[X^2] - E[X]^2), need to subtract average squared then take square root.
/// If \p size == 0, then we assume that num has already been applied for divisor everywhere.
/// If \p size < 0, then we use size based on #Stat::num.
void finalizeStats(StatVector& stats, const int size) {
  const int avg_ind = static_cast<int>(Stat::avg);
  const int stddev_ind = static_cast<int>(Stat::stddev);
  const double divisor = (size > 0) ? static_cast<double>(size) : (size == 0 ? 1. : stats[static_cast<int>(Stat::num)]);

  stats[avg_ind] = stats[static_cast<int>(Stat::total)] / divisor;
  stats[stddev_ind] = (stats[stddev_ind]) / divisor - stats[avg_ind] * stats[avg_ind];
  stats[stddev_ind] = (stats[stddev_ind] > 0) ? std::sqrt(stats[stddev_ind]) : 0.;
} /* finalizeStats */

/// @details Need to adjust stddev for previous values.
/// We have currently stats[stddev] = sqrt( E[X_{old}^2] - E[X_{old}]^2 ).
/// Can square to get OLD_VARIANCE = E[X_{old}^2] - E[X_{old}]^2.
/// OLD_AVG = E[X_{old}] = old_total / old_size (should equal previous value in stats[avg]).
/// E[X_{old}^2] = OLD_VARIANCE + E[X_{old}]^2 = OLD_VARIANCE + OLD_AVG^2.
/// OLD_TOTAL_SOS = OLD_VARIANCE + OLD_AVG^2.
/// Then we add new square values to OLD_TOTAL_SOS and subtract new average squared to get new variance
void unfinalizeStats(StatVector& stats, const int prev_size) {
  const int OLD_SIZE = (prev_size > 0) ? prev_size : stats[static_cast<int>(Stat::num)];
  if (OLD_SIZE == 0) {
    return;
  }

  const int avg_ind = static_cast<int>(Stat::avg);
  const int stddev_ind = static_cast<int>(Stat::stddev);

  // Save previous total, avg, stddev
  const double OLD_TOTAL = stats[static_cast<int>(Stat::total)];
  const double OLD_AVG = stats[avg_ind];
  const double OLD_STDDEV = stats[stddev_ind]; // sqrt(E[X_{old}^2] - E[X_{old}]^2)

  assert ( isVal( OLD_TOTAL / OLD_SIZE, OLD_AVG, 1e-3) );

  const double OLD_VARIANCE = OLD_STDDEV * OLD_STDDEV; // E[X_{old}^2] - E[X_{old}]^2
  const double OLD_TOTAL_SOS = (OLD_VARIANCE + OLD_AVG * OLD_AVG) * OLD_SIZE; // total sum X_{old}^2

  stats[avg_ind] *= OLD_SIZE;
  stats[stddev_ind] = OLD_TOTAL_SOS;
} /* finalizeStats */

/**
 * @details Creates temporary file (in /tmp) so that it can be accessed later.
 * It does not delete the file.
 */
void createTmpFilename(std::string& f_name,
    const std::string add_ext) {
//  if (f_name.empty()) {
    // Generate temporary file name
    char template_name[] = "/tmp/tmpvpcXXXXXX";

    int fd;
    if ((fd = mkstemp(template_name)) == -1) {
      error_msg(errorstring,
          "Could not generate temp file with error %s.\n",
          strerror(errno));
      throw std::logic_error(errorstring);
    }

    f_name = template_name;
    if (f_name.empty()) {
      error_msg(errorstring, "Could not generate temp file.\n");
      throw std::logic_error(errorstring);
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
