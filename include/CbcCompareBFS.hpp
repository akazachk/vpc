/**
 * @file CbcCompareBFS.hpp
 *
 * @brief Breadth-first search through the tree
 *
 * @author Edwin (edited by A. M. Kazachkov)
 * @date 2009-11-24 (amk edited 2017-04-14)
 */

// The original license for this code is as follows:
// $Id: CbcCompareBFS.hpp 1899 2013-04-09 18:12:08Z stefan $
// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

//Edwin 11/24/09 carved out of CbcCompareActual

#ifndef CbcCompareBFS_H
#define CbcCompareBFS_H


//#############################################################################
/*  These are alternative strategies for node traversal.
    They can take data etc for fine tuning

    At present the node list is stored as a heap and the "test"
    comparison function returns true if node y is better than node x.

*/
#include "CbcNode.hpp"
#include "CbcCompareBase.hpp"
#include "CbcCompare.hpp"
class CbcModel;
// This is default before first solution
class CbcCompareBFS : public CbcCompareBase {
public:
    // Default Constructor
    CbcCompareBFS () ;

    ~CbcCompareBFS();
    // Copy constructor
    CbcCompareBFS ( const CbcCompareBFS &rhs);

    // Assignment operator
    CbcCompareBFS & operator=( const CbcCompareBFS& rhs);

    /// Clone
    virtual CbcCompareBase * clone() const;
    /// Create C++ lines to get to current state
    virtual void generateCpp( FILE * fp);

    // This returns true if the depth of node y is greater than depth of node x
    virtual bool test (CbcNode * x, CbcNode * y);
};

#endif

