/**
 * @file VPCSolverInterface.cpp
 * @author A. M. Kazachkov
 * @brief Wrapper for VPC-related solver calls
 */
#include "VPCSolverInterface.hpp"
#include "utility.hpp"
#include "CglVPC.hpp"
#include "Disjunction.hpp"
#include "OsiProblemData.hpp"
#include "SolverInterface.hpp"
#include "VPCParameters.hpp"

#ifdef USE_COIN
#include <OsiSolverInterface.hpp>

std::vector<SparseCut> convertCutsToSparseCuts(const OsiCuts* const cuts) {
  if (!cuts) {
    return std::vector<SparseCut>();
  }
  std::vector<SparseCut> cuts_out(cuts->sizeCuts());
  for (int cut_ind = 0; cut_ind < cuts->sizeCuts(); cut_ind++) {
    const OsiRowCut* cut = cuts->rowCutPtr(cut_ind);
    int num_elem = cut->row().getNumElements();
    const int* ind = cut->row().getIndices();
    const double* coeff = cut->row().getElements();
    cuts_out[cut_ind].num_elem = num_elem;
    cuts_out[cut_ind].index.assign(ind, ind + num_elem);
    cuts_out[cut_ind].coeff.assign(coeff, coeff + num_elem);
    cuts_out[cut_ind].rhs = cut->rhs();
  }
  return cuts_out;
} /* convertCutsToSparseCuts */

std::vector<DenseCut> convertCutsToDenseCuts(const OsiCuts* const cuts, const int n) {
  if (!cuts) {
    return std::vector<DenseCut>();
  }
  std::vector<DenseCut> cuts_out(cuts->sizeCuts());
  for (int cut_ind = 0; cut_ind < cuts->sizeCuts(); cut_ind++) {
    const OsiRowCut* cut = cuts->rowCutPtr(cut_ind);
    int num_elem = cut->row().getNumElements();
    const int* ind = cut->row().getIndices();
    const double* coeff = cut->row().getElements();
    cuts_out[cut_ind].coeff.resize(n, 0.0);
    for (int i = 0; i < num_elem; i++) {
      const int curr_ind = ind[i];
      const int curr_val = coeff[i];
      cuts_out[cut_ind].coeff[curr_ind] = curr_val;
    }
    cuts_out[cut_ind].rhs = cut->rhs();
  }
  return cuts_out;
} /* convertCutsToDenseCuts */
#endif // USE_COIN

VPCSolverInterface::VPCSolverInterface() {
  initialize();
} /* default constructor */

VPCSolverInterface::VPCSolverInterface(const VPCSolverInterface &source) {
  initialize(&source);
} /* copy constructor */

VPCSolverInterface::~VPCSolverInterface() {
	free();
} /* destructor */

VPCSolverInterface &VPCSolverInterface::operator=(const VPCSolverInterface &source) {
  if (this != &source) {
    initialize(&source);
  }
	return *this;
} /* assignment operator */

VPCSolverInterface* VPCSolverInterface::clone() const {
	return new VPCSolverInterface(*this);
} /* clone */

void VPCSolverInterface::setParams(const VPCParametersNamespace::VPCParameters* const param) {
  if (params)
    delete params;
  this->params = new VPCParametersNamespace::VPCParameters(*param);
} /* setParams */

void VPCSolverInterface::flipObj(const int sense) {
  // Make sure we are doing a minimization problem; this is just to make later
  // comparisons simpler (i.e., a higher LP obj after adding the cut is better).
  if (!isVal(solver->getObjSense(), sense)) {
    if (sense == 1) {
      printf(
          "\n## Detected maximization problem. Negating objective function to make it minimization. ##\n");
    } else {
      printf(
          "\n## Detected minimization problem. Negating objective function to make it maximization. ##\n");
    }

    solver->setObjSense((double) sense);
    const double *obj = solver->getObjCoefficients();
    for (int col = 0; col < solver->getNumCols(); col++) {
      solver->setObjCoeff(col, -1. * obj[col]);
    }
    double objOffset = 0.;
    solver->getDblParam(OsiDblParam::OsiObjOffset, objOffset);
    if (objOffset != 0.) {
      solver->setDblParam(OsiDblParam::OsiObjOffset, -1. * objOffset);
    }
  }
} /* setObjSense */

void VPCSolverInterface::setDisjunction(const Disjunction *const disj) {
  if (disj)
    delete disj;
  this->disj = disj->clone();
} /* setDisjunction */

void VPCSolverInterface::load(std::string filename) {
  std::string dir, instname, in_file_ext;
  parseFilename(dir, instname, in_file_ext, filename, params->logfile);

  int status = 0;
  if (in_file_ext.compare("lp") == 0) {
#ifdef TRACE
    printf("\n## Reading LP file. ##\n");
#endif
    status = solver->readLp(filename.c_str());
  } else {
    if (in_file_ext.compare("mps") == 0) {
#ifdef TRACE
      printf("\n## Reading MPS file. ##\n");
#endif
      status = solver->readMps(filename.c_str());
    } else {
      try {
#ifdef TRACE
        printf("\n## Reading MPS file. ##\n");
#endif
        status = solver->readMps(filename.c_str());
      } catch (std::exception &e) {
        error_msg(errorstring, "Unrecognized extension: %s.\n",
            in_file_ext.c_str());
        writeErrorToLog(errorstring, params->logfile);
        exit(1);
      }
    }
  } // read file
  if (status < 0) {
    error_msg(errorstring, "Unable to read in file %s.\n", filename.c_str());
    writeErrorToLog(errorstring, params->logfile);
    exit(1);
  }

  // Make sure we are doing a minimization problem; this is just to make later
  // comparisons simpler (i.e., a higher LP obj after adding the cut is better).
  this->flipObj();
} /* load (filename) */

#ifdef USE_COIN
void VPCSolverInterface::load(const OsiSolverInterface* const si) {
  if (solver)
    delete solver;
  solver = si->clone();

  // Make sure we are doing a minimization problem; this is just to make later
  // comparisons simpler (i.e., a higher LP obj after adding the cut is better).
  this->flipObj();
} /* load (OsiSolverInterface) */
#endif // USE_COIN

/**
 * @details If we wish to replace Osi solver with something else,
 * we need to implement loadProblem, setInteger, setColLower, setColUpper, setInteger
 */
void VPCSolverInterface::load(OsiProblemData *data) {
  if (!solver)
    solver = new SolverInterface;
  solver->loadProblem(data->numcols, data->numrows, data->start, data->index,
      data->value, data->collb, data->colub, data->obj, data->rowlb,
      data->rowub);
  solver->setObjSense(data->objsense);
  for (int i = 0; i < data->numcols; i++) {
    // setColumnType is not implemented in earlier versions of OsiSolverInterface, so we revert to older methods
//      solver->setColumnType(i, data->vartype[i]);
    const char vartype = data->vartype[i];

    if (vartype == OsiVarTypeChar[OsiVarType::BINARY]) {
      solver->setInteger(i);
      solver->setColLower(i, 0);
      solver->setColUpper(i, 1);
    }
    if (vartype == OsiVarTypeChar[OsiVarType::INTEGER]) {
      solver->setInteger(i);
    }
  }

  // Make sure we are doing a minimization problem; this is just to make later
  // comparisons simpler (i.e., a higher LP obj after adding the cut is better).
  this->flipObj();
} /* load (OsiProblemData) */

void VPCSolverInterface::solve() {
#ifdef USE_COIN
  solver->initialSolve();
#else
  fprintf(stderr, "No solver defined.\n");
#endif
} /* solve */

void VPCSolverInterface::resolve() {
#ifdef USE_COIN
  solver->resolve();
#else
  fprintf(stderr, "No solver defined.\n");
#endif
} /* resolve */

void VPCSolverInterface::generateCuts(
    /// [in] optional: set of user objectives for PRLP (in structural space), to be tried before default objectives
    const std::vector<std::vector<double> >* const user_objectives) {
  if (!solver->isProvenOptimal()) {
    error_msg(errorstring, "VPCSolverInterface::generateCuts: Solver not proven optimal.\n");
    throw(errorstring);
  }
  CglVPC gen(*params);
  if (user_objectives) { gen.setUserObjectives(*user_objectives); }
  gen.generateCuts(*solver, *cuts); // solution may change slightly due to enable factorization called in getProblemData...
  if (gen.disj()) {
    if (this->disj) { delete this->disj; }
    this->disj = gen.disj()->clone();
  }
} /* generateCuts(user_objectives) */

int VPCSolverInterface::applyCuts(
  bool* cutsToAdd,
  bool clear_cuts
) {
  int num_applied = 0;
  int num_cuts = 0;
  if (cutsToAdd != NULL) {
    int num_before = solver->getNumRows();
    for (int i = 0; i < cuts->sizeRowCuts(); ++i) {
      if (cutsToAdd[i]) {
        num_cuts++;
        solver->applyRowCuts(1, cuts->rowCutPtr(i));
      }
    }
    num_applied = solver->getNumRows() - num_before;
    if (clear_cuts) {
      delete cuts;
      cuts = new OsiCuts;
    }
  } else {
    num_cuts = cuts->sizeCuts();
    OsiSolverInterface::ApplyCutsReturnCode code;
    code = solver->applyCuts(*(cuts));
    num_applied = code.getNumApplied();
    if (clear_cuts) {
      delete cuts;
      cuts = new OsiCuts;
    }
  }
#ifdef TRACE
    printf("\n## Applied %d/%d cuts. ##\n", num_applied, num_cuts);
#endif
  return num_applied;
} /* applyCuts */

void VPCSolverInterface::save(std::string filename) {
  solver->writeMps(filename.c_str());
} /* save (to filename) */

int VPCSolverInterface::numCuts() {
#ifdef USE_COIN
  if (cuts)
    return cuts->sizeCuts();
  else
#endif
    return 0;
} /* numCuts */

#ifdef USE_COIN
const OsiCuts* const VPCSolverInterface::getCuts() {
  return cuts;
} /* getCuts */
#endif // USE_COIN

std::vector<SparseCut> VPCSolverInterface::getCutsSparse() {
  return convertCutsToSparseCuts(cuts);
} /* getCutsSparse */

std::vector<DenseCut> VPCSolverInterface::getCutsDense() {
  return convertCutsToDenseCuts(cuts, solver->getNumCols());
} /* getCutsDense */

double VPCSolverInterface::getObjValue() {
  return solver->getObjValue();
}

void VPCSolverInterface::initialize(const VPCSolverInterface *const source,
    const VPCParametersNamespace::VPCParameters *const param) {
  if (param != NULL) {
    setParams(param);
  }
  if (source != NULL) {
    free();
#ifdef USE_COIN
    solver = source->solver->clone();
    cuts = new OsiCuts(*source->cuts);
#else
    fprintf(stderr, "No solver defined.\n");
#endif
    disj = source->disj->clone();
    setParams(source->params);
  } else {
#ifdef USE_COIN
    solver = new SolverInterface;
    cuts = new OsiCuts;
#else
    fprintf(stderr, "No solver defined.\n");
#endif
#ifdef USE_CPLEX_SOLVER
    CPXENVptr env = dynamic_cast<OsiCpxSolverInterface*>(solver)->getEnvironmentPtr();
    CPXsetintparam(env, CPXPARAM_Threads, 1);
#ifndef TRACE
    CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_OFF);
#endif
#endif
    disj = NULL;
    params = new VPCParametersNamespace::VPCParameters;
  }
} /* initialize */

void VPCSolverInterface::free() {
#ifdef USE_COIN
  if (solver) {
    delete solver;
    solver = NULL;
  }
  if (cuts) {
    delete cuts;
    cuts = NULL;
  }
#endif
  if (disj) {
    delete disj;
    disj = NULL;
  }
  if (params) {
    delete params;
    params = NULL;
  }
} /* free */
