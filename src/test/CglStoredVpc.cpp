/**
 * @file CglStoredVpc.cpp
 * @author Sean Kelley
 * @date 2023-05-15
 */


#include "CglStoredVpc.hpp"

//-------------------------------------------------------------------
// Apply Stored cuts
//-------------------------------------------------------------------
void CglStoredVpc::generateCuts(const OsiSolverInterface &si, OsiCuts &cs,
                                 const CglTreeInfo info) const
{
  // Get basic problem information
  const double *solution = si.getColSolution();
  int numberRowCuts = cuts_.sizeRowCuts();
  for (int i = 0; i < numberRowCuts; i++) {
    const OsiRowCut *rowCutPointer = cuts_.rowCutPtr(i);
    double violation = rowCutPointer->violated(solution);
    if (violation >= requiredViolation_)
      cs.insert(*rowCutPointer);
  }
}

/**@name Cut stuff */
OsiRowCut * CglStoredVpc::mutableRowCutPointer(int index)
{
  return cuts_.rowCutPtr(index);
}

//-------------------------------------------------------------------
// Default Constructor
//-------------------------------------------------------------------
CglStoredVpc::CglStoredVpc()
    : CglStored()
{
}

//-------------------------------------------------------------------
// Copy constructor
//-------------------------------------------------------------------
CglStoredVpc::CglStoredVpc(const CglStoredVpc &source)
    : CglStored(source)
{
}

//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
CglCutGenerator * CglStoredVpc::clone() const
{
  return new CglStoredVpc(*this);
}

//-------------------------------------------------------------------
// Destructor
//-------------------------------------------------------------------
CglStoredVpc::~CglStoredVpc()
{
}

//----------------------------------------------------------------
// Assignment operator
//-------------------------------------------------------------------
CglStoredVpc & CglStoredVpc::operator=(const CglStoredVpc &rhs)
{
  if (this != &rhs) {
    CglStored::operator=(rhs);
  }
  return *this;
}

OsiCuts CglStoredVpc::getCuts() {
  return cuts_;
}