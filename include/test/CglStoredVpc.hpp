/**
 * @file CglStoredVpc.hpp
 * @author Sean Kelley
 * @date 2023-05-15
 */


// coin-or
#include "OsiSolverInterface.hpp"
#include "OsiCuts.hpp"
#include "CglTreeInfo.hpp"
#include "CglStored.hpp"


/** Stored Cut Generator Class */
class CglStoredVpc : public CglStored {

public:
  /**@name Generate Cuts */
  /** Generate Mixed Integer Stored cuts for the model of the
      solver interface, si.

      Insert the generated cuts into OsiCut, cs.

      This generator just looks at previously stored cuts
      and inserts any that are violated by enough
  */
  virtual void generateCuts(const OsiSolverInterface &si, OsiCuts &cs,
                            const CglTreeInfo info = CglTreeInfo()) const;

  /**@name Cut stuff */
  OsiRowCut *mutableRowCutPointer(int index);

  /// Default constructor
  CglStoredVpc();

  /// Copy constructor
  CglStoredVpc(const CglStoredVpc &rhs);

  /// Clone
  virtual CglCutGenerator *clone() const;

  /// Assignment operator
  CglStoredVpc &
  operator=(const CglStoredVpc &rhs);

  /// Destructor
  virtual ~CglStoredVpc();

  /// Access Cuts
  OsiCuts getCuts();

};