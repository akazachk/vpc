#pragma once

#include <string>
#include <vector>

/* Data structure to pass to/from Disjunctive cut generator the problem over which to generate the cuts.
   This is to conform to OsiSolverInterface way of passing a problem */
enum OsiVarType {
  CONTINUOUS = 0,
  BINARY = 1,
  INTEGER = 2
};
const char OsiVarTypeChar[3] {
  'c', 'b', 'i'
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
   OsiVarType* vartype;
};

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
class OsiSolverInterface;
class OsiCuts;
class Disjunction;
struct VPCParameters;

std::vector<SparseCut> convertCutsToSparseCuts(const OsiCuts* const cuts);
std::vector<DenseCut> convertCutsToDenseCuts(const OsiCuts* const cuts, const int n);

/**
 * Class VPCSolverInterface
 * Provides way of generating rounds of cuts 
 * that can easily be interacted with by external code
 */
class VPCSolverInterface {
public:
	OsiSolverInterface *solver;
	OsiCuts *cuts;
  Disjunction *disj;
  VPCParameters *params;

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
  void setParams(const VPCParameters* const param);

  /** Load problem methods */
	virtual void load(std::string fullfilename);
  virtual void load(OsiSolverInterface* si);
  virtual void load(OsiProblemData* data);

  /** Write LP with cuts */
  virtual void save(std::string filename);

  /** Return cuts in non-OSI format */
  virtual int numCuts();
  virtual const OsiCuts* const getCuts();
  virtual std::vector<SparseCut> getCutsSparse();
  virtual std::vector<DenseCut> getCutsDense();

  /** Set disjunction */
  virtual void setDisjunction(const Disjunction* const disj);

	/** Generate cuts and put them in `OsiCuts* cuts` */
	virtual void generateCuts();

	/** 
	 * Add the (previously generated) cuts to the solver
	 */
	virtual int applyCuts(bool* cutsToAdd, bool clear_cuts = false);

protected:
  void initialize(const VPCSolverInterface *const source = NULL,
      const VPCParameters *const param = NULL);
  void free();
}; /* VPCSolverInterface */
