/**
 * @file TestVPCEventHandler.cpp
 * @author Sean Kelley
 * @date 2023-09-20
 */

#define CATCH_CONFIG_MAIN

// standard library
#include <cstdlib> // abs
#include <vector> // vector

// unit test library
#include "catch.hpp"

// coin-or modules
#include "OsiClpSolverInterface.hpp" // OsiClpSolverInterface
#include "OsiCuts.hpp" // OsiCuts

// project modules
#include "Disjunction.hpp" // DisjExitReason
#include "PartialBBDisjunction.hpp" // PartialBBDisjunction
#include "SolverInterface.hpp" // SolverInterface


// --------------------- test current behavior remains -------------------------
TEST_CASE("Test saveInformation", "[VPCEventHandler::saveInformation]") {

  // parameters
  VPCParametersNamespace::VPCParameters vpc_params;
  vpc_params.set(VPCParametersNamespace::DISJ_TERMS, 4);
  vpc_params.set(VPCParametersNamespace::MODE, 0);  // partial BB tree

  // solver
  OsiClpSolverInterface si;
  SolverInterface* solver;
  si.readMps("../test/p0033.mps");
  si.initialSolve();
  solver = const_cast<SolverInterface*>(dynamic_cast<const SolverInterface*>(&si));

  SECTION( "Test without strong branching terms" ) {

    // make disjunction
    PartialBBDisjunction disj = PartialBBDisjunction(vpc_params);
    disj.prepareDisjunction(solver);

    // check that we have the right number of terms
    REQUIRE(disj.num_terms == 4);
    REQUIRE(disj.terms.size() == 4);
    REQUIRE(disj.name == "(-x15 >= 0;  V (+x10 >= 1)) "
                         "V (+x15 >= 1;  V (+x10 >= 1)) "
                         "V (-x9 >= 0;  V (+x0 >= 1; -x10 >= 0; -x30 >= 0)) "
                         "V (+x9 >= 1;  V (+x0 >= 1; -x10 >= 0; -x30 >= 0))");
    REQUIRE(disj.common_changed_var.size() == 1);
    REQUIRE(disj.common_changed_bound.size() == 1);
    REQUIRE(disj.common_changed_value.size() == 1);
    REQUIRE(disj.common_changed_var[0] == 14);
    REQUIRE(disj.common_changed_bound[0] == 1);
    REQUIRE(disj.common_changed_value[0] == 0);
    REQUIRE(disj.common_ineqs.size() == 0);
    REQUIRE(disj.integer_sol.size() == 0);
    REQUIRE(disj.integer_obj > 1e300);
  }

  SECTION( "Test with strong branching terms" ) {

    // make sure to save the full tree
    vpc_params.set(VPCParametersNamespace::SAVE_FULL_TREE, 1);

    // make disjunction
    PartialBBDisjunction disj = PartialBBDisjunction(vpc_params);
    disj.prepareDisjunction(solver);

    // check that we have the right number of terms
    REQUIRE(disj.num_terms == 7);
    REQUIRE(disj.terms.size() == 7);

    // check that we have no common terms anymore
    REQUIRE(disj.common_changed_var.size() == 0);
    REQUIRE(disj.common_changed_bound.size() == 0);
    REQUIRE(disj.common_changed_value.size() == 0);

    // check terms
    std::vector<std::vector<int> > var = {
        {14}, {14, 10, 15}, {14, 10, 15}, {14, 10, 0}, {14, 10, 0, 30},
        {14, 0, 10, 30, 9}, {14, 0, 10, 30, 9}
    };
    std::vector<std::vector<int> > bound = {
        {0}, {1, 0, 1}, {1, 0, 0}, {1, 1, 1}, {1, 1, 0, 0}, {1, 0, 1, 1, 1}, {1, 0, 1, 1, 0}
    };
    std::vector<std::vector<double> > value = {
        {1}, {0, 1, 0}, {0, 1, 1}, {0, 0, 0}, {0, 0, 1, 1}, {0, 1, 0, 0, 0}, {0, 1, 0, 0, 1}
    };
    for (int i = 0; i < disj.num_terms; i++) {
      REQUIRE(disj.terms[i].changed_var == var[i]);
      REQUIRE(disj.terms[i].changed_bound == bound[i]);
      REQUIRE(disj.terms[i].changed_value == value[i]);
    }

    // check misc
    REQUIRE(disj.common_ineqs.size() == 0);
    REQUIRE(disj.integer_sol.size() == 0);
    REQUIRE(disj.integer_obj > 1e300);
  }
  // need to check best and worst objective values
  // need to test that we still get the same cuts afterwards
}


