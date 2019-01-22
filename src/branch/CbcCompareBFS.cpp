// $Id: CbcCompareBFS.cpp 1899 2013-04-09 18:12:08Z stefan $
// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

//Edwin 11/24/09 carved out of CbcCompareActual
//AMK 04/14/17

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <cfloat>
//#define CBC_DEBUG

#include "CbcMessage.hpp"
#include "CbcModel.hpp"
#include "CbcTree.hpp"
#include "CbcCompareActual.hpp"
#include "CoinError.hpp"
#include "CbcCompareBFS.hpp"
/** Default Constructor

*/
CbcCompareBFS::CbcCompareBFS ()
        : CbcCompareBase()
{
    test_ = this;
}

// Copy constructor
CbcCompareBFS::CbcCompareBFS ( const CbcCompareBFS & rhs)
        : CbcCompareBase(rhs)

{
}

// Clone
CbcCompareBase *
CbcCompareBFS::clone() const
{
    return new CbcCompareBFS(*this);
}

// Assignment operator
CbcCompareBFS &
CbcCompareBFS::operator=( const CbcCompareBFS & rhs)
{
    if (this != &rhs) {
        CbcCompareBase::operator=(rhs);
    }
    return *this;
}

// Destructor
CbcCompareBFS::~CbcCompareBFS ()
{
}

// Returns true if y better than x
bool
CbcCompareBFS::test (CbcNode * x, CbcNode * y)
{
    int testX = x->depth();
    int testY = y->depth();
    if (testX != testY)
        return testX > testY;
    else
        return equalityTest(x, y); // so ties will be broken in consistent manner
}
// Create C++ lines to get to current state
void
CbcCompareBFS::generateCpp( FILE * fp)
{
    fprintf(fp, "0#include \"CbcCompareActual.hpp\"\n");
    fprintf(fp, "3  CbcCompareBFS compare;\n");
    fprintf(fp, "3  cbcModel->setNodeComparison(compare);\n");
}

