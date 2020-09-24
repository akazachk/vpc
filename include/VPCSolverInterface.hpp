/**
 * @file VPCSolverInterface.hpp
 * @author A. M. Kazachkov
 * @brief Wrapper for VPC-related solver calls
 */
#pragma once

#include <string>
#include <vector>

/// @brief For passing sparse cuts out from VPC code
struct SparseCut
{
  int num_elem;
  std::vector<int> index;
  std::vector<double> coeff;
  double rhs;
};

/// @brief For passing dense cuts out from VPC code
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
 * @brief Provides way of generating rounds of cuts 
 * that can easily be interacted with by external code
 */
class VPCSolverInterface {
public:
#ifdef USE_COIN
	OsiSolverInterface *solver; ///< pointer to underlying OsiSolver
	OsiCuts *cuts; ///< pointer to any cuts we generate
#endif
  Disjunction *disj; ///< disjunction used for cut generation
  VPCParametersNamespace::VPCParameters *params; ///< instance parameters

  /// @brief Default constructor
	VPCSolverInterface();

	/// @brief Copy constructor
	VPCSolverInterface(const VPCSolverInterface& source);

	/// @brief Destructor
	virtual ~VPCSolverInterface();

	/// @brief Assignment operator
	VPCSolverInterface& operator=(const VPCSolverInterface& source);

	/// @brief Clone
	virtual VPCSolverInterface *clone() const;

  /// @brief Set params based on VPCParameters
  void setParams(const VPCParametersNamespace::VPCParameters* const param);

  /// @brief Flip objective
  void flipObj(const int sense = 1);

  /// @brief Set disjunction
  virtual void setDisjunction(const Disjunction* const disj);

  ///@{
  /// @name Load problem methods

  /// @brief Load from string
	virtual void load(std::string fullfilename);
#ifdef USE_COIN
  /// @brief Load from OsiSolverInterface
  virtual void load(const OsiSolverInterface* const si);
#endif
  /// @brief Load from OsiProblemData
  virtual void load(OsiProblemData* data);
  ///@}

  /// @brief Solve LP relaxation
  virtual void solve();
  /// @brief Resolve LP relaxation
  virtual void resolve();

  /// @brief Generate cuts and put them in #cuts
  virtual void generateCuts();

  /// @brief Add the (previously generated) cuts to the solver
  virtual int applyCuts(bool* cutsToAdd = NULL, bool clear_cuts = false);

  /// @brief Write LP with cuts
  virtual void save(std::string filename);

  /// @brief Return cuts in non-Osi format
  virtual int numCuts();
#ifdef USE_COIN
  /// @brief Return cuts in Osi format
  virtual const OsiCuts* const getCuts();
#endif
  /// @brief Get vector of sparse cuts
  virtual std::vector<SparseCut> getCutsSparse();
  /// @brief Get vector of dense cuts
  virtual std::vector<DenseCut> getCutsDense();

  /// @brief Query solver
  virtual double getObjValue();

protected:
  /// @brief Initialize class members (copy from source if provided)
  void initialize(const VPCSolverInterface *const source = NULL,
      const VPCParametersNamespace::VPCParameters *const param = NULL);
  /// @brief Free memory
  void free();
}; /* VPCSolverInterface */
