#pragma once

#include <string>
#include <vector>

/** For passing out cuts, we have dense and sparse versions */
struct SparseCut
{
  int num_elem;
  std::vector<int> index;
  std::vector<double> coeff;
  double rhs;
};

struct DenseCut
{
  std::vector<double> coeff;
  double rhs;
};

// Forward declarations
class Disjunction;
namespace VPCParametersNamespace {
  struct VPCParameters;
}
struct OsiProblemData;

#ifdef USE_COIN
class OsiCuts;
class OsiSolverInterface;

std::vector<SparseCut> convertCutsToSparseCuts(const OsiCuts* const cuts);
std::vector<DenseCut> convertCutsToDenseCuts(const OsiCuts* const cuts, const int n);
#endif

/**
 * Class VPCSolverInterface
 * Provides way of generating rounds of cuts 
 * that can easily be interacted with by external code
 */
class VPCSolverInterface {
public:
#ifdef USE_COIN
	OsiSolverInterface *solver;
	OsiCuts *cuts;
#endif
  Disjunction *disj;
  VPCParametersNamespace::VPCParameters *params;

	/** Default constructor */
	VPCSolverInterface();

	/** Copy constructor */
	VPCSolverInterface(const VPCSolverInterface& source);

	/** Destructor */
	virtual ~VPCSolverInterface();

	/** Assignment operator */
	VPCSolverInterface& operator=(const VPCSolverInterface& source);

	/** Clone */
	virtual VPCSolverInterface *clone() const;

  /** Set params based on VPCParameters */
  void setParams(const VPCParametersNamespace::VPCParameters* const param);

  /** Flip objective */
  void flipObj(const int sense = 1);

  /** Set disjunction */
  virtual void setDisjunction(const Disjunction* const disj);

  /** Load problem methods */
	virtual void load(std::string fullfilename);
#ifdef USE_COIN
  virtual void load(const OsiSolverInterface* const si);
#endif
  virtual void load(OsiProblemData* data);

  /** Solve LP relaxation */
  virtual void solve();
  /** Resolve LP relaxation */
  virtual void resolve();

  /** Generate cuts and put them in `OsiCuts* cuts` */
  virtual void generateCuts();

  /** Add the (previously generated) cuts to the solver */
  virtual int applyCuts(bool* cutsToAdd = NULL, bool clear_cuts = false);

  /** Write LP with cuts */
  virtual void save(std::string filename);

  /** Return cuts in non-OSI format */
  virtual int numCuts();
#ifdef USE_COIN
  virtual const OsiCuts* const getCuts();
#endif
  virtual std::vector<SparseCut> getCutsSparse();
  virtual std::vector<DenseCut> getCutsDense();

  /** Query solver */
  virtual double getObjValue();

protected:
  void initialize(const VPCSolverInterface *const source = NULL,
      const VPCParametersNamespace::VPCParameters *const param = NULL);
  void free();
}; /* VPCSolverInterface */
