/**
 * @file TestCompleteDisjunction.cpp
 * @author Sean Kelley
 * @date 2023-05-04
 */

#define CATCH_CONFIG_MAIN

// standard library
#include <cstdlib> // abs

// unit test library
#include "catch.hpp"

// coin-or modules
#include "OsiClpSolverInterface.hpp" // OsiClpSolverInterface

// project modules
#include "CompleteDisjunction.hpp" // CompleteDisjunction
#include "Disjunction.hpp" // DisjExitReason


// --------------------------- constructor tests -------------------------------
TEST_CASE("Test default constructor", "[VwsSolverInterface::VwsSolverInterface]") {

  // Ensure we set up as expected
  CompleteDisjunction disj;
  REQUIRE(!disj.partialDisj);
}

