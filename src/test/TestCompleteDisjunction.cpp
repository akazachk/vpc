/**
 * @file TestCompleteDisjunction.cpp
 * @author Sean Kelley
 * @date 2023-05-04
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
#include "CompleteDisjunction.hpp" // CompleteDisjunction
#include "CglVpc.hpp" // CglVpc
#include "Disjunction.hpp" // DisjExitReason
#include "PartialBBDisjunction.hpp" // PartialBBDisjunction


void setup(OsiClpSolverInterface & si, PartialBBDisjunction & partialDisj){
  // parameters
  VPCParametersNamespace::VPCParameters vpc_params;
  vpc_params.set(VPCParametersNamespace::DISJ_TERMS, 4);
  vpc_params.set(VPCParametersNamespace::MODE, 0);

  // solver
  si.readMps("../test/bm23.mps");
  si.initialSolve();

  // cuts
  OsiCuts disjCuts;

  // generator
  CglVPC gen = CglVPC(vpc_params);
  gen.generateCuts(si, disjCuts);

  // partial disjunction - 0, 13, 14, 16, 22, 25 ()
  partialDisj = *dynamic_cast<PartialBBDisjunction*>(gen.disj());
  std::vector<int> vars({14, 22});
  partialDisj.common_changed_var = vars;
  std::vector<int> bounds({0, 1});
  partialDisj.common_changed_bound = std::vector<int>({0, 1});
  std::vector<double> vals({1, 0});
  partialDisj.common_changed_value = vals;
}

// --------------------------- constructor tests -------------------------------
TEST_CASE("Test default constructor", "[CompleteDisjunction::CompleteDisjunction]") {
  // Ensure we set up as expected
  CompleteDisjunction disj;
  REQUIRE(!disj.partialDisj);
}

// ----------------------- prepare disjunction tests ---------------------------
TEST_CASE("Test prepareDisjunction", "[CompleteDisjunction::prepareDisjunction]") {

  OsiClpSolverInterface si;
  PartialBBDisjunction partialDisj;
  setup(si, partialDisj);

  CompleteDisjunction disj;
  disj.prepareDisjunction(&si, &partialDisj);

  // check that we have the right number of terms - 2 leaves for fixes and 4 for partial BB tree
  REQUIRE(disj.num_terms == 6);
  REQUIRE(disj.terms.size() == 6);
  REQUIRE(disj.name == "(x14 >= 1; x22 >= 1) V (x14 <= 0) V "
                       "(x14 >= 1; x22 <= 0; x16 >= 1; x13 <= 0) "
                       "V (x14 >= 1; x22 <= 0; x16 >= 1; x13 >= 1) "
                       "V (x14 >= 1; x22 <= 0; x16 <= 0; x2 >= 1) "
                       "V (x14 >= 1; x22 <= 0; x16 <= 0; x2 <= 0)");
  REQUIRE(disj.common_changed_var.size() == 0);
  REQUIRE(disj.common_changed_bound.size() == 0);
  REQUIRE(disj.common_changed_value.size() == 0);
  REQUIRE(disj.common_ineqs.size() == 0);
  REQUIRE(disj.integer_sol.size() == 0);
  REQUIRE(disj.integer_obj > 1e300);
}

TEST_CASE("Test Adding Terms", "[CompleteDisjunction::addCommonFixedTerms]"
                               "[CompleteDisjunction::addDisjunctionFixedTerms]") {

  OsiClpSolverInterface si;
  PartialBBDisjunction partialDisj;
  setup(si, partialDisj);

  CompleteDisjunction disj;
  disj.partialDisj = partialDisj.clone();
  std::shared_ptr<SolverInterface> solver = getSolver(&si);
  disj.addCommonFixedTerms(solver.get());

  // two fixes means three leaves, but one is left out because it will be further expanded
  REQUIRE(disj.num_terms == 2);
  REQUIRE(disj.terms.size() == 2);

  disj.addDisjunctionFixedTerms(solver.get());

  // two fixes each => three leaves each => 12 leaves across disjunctions => 14 leaves total
  REQUIRE(disj.num_terms == 14);
  REQUIRE(disj.terms.size() == 14);

  solver->unmarkHotStart();
  solver->disableFactorization();
}
