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
#include "CglVPC.hpp" // CglVPC
#include "Disjunction.hpp" // DisjExitReason
#include "PartialBBDisjunction.hpp" // PartialBBDisjunction
#include "SolverInterface.hpp" // SolverInterface


// --------------------- test current behavior remains -------------------------
TEST_CASE("Test saveInformation", "[VPCEventHandler::saveInformation]") {

  // parameters
  VPCParametersNamespace::VPCParameters vpc_params;
  vpc_params.set(VPCParametersNamespace::DISJ_TERMS, 4);
  vpc_params.set(VPCParametersNamespace::MODE, 0);  // partial BB tree
  vpc_params.set(VPCParametersNamespace::PARTIAL_BB_KEEP_PRUNED_NODES, 1);

  // solver
  OsiClpSolverInterface si;
  SolverInterface* solver;
  si.readMps("../test/bm23.mps");
  si.initialSolve();
  std::vector<double> sol(si.getColSolution(), si.getColSolution() + si.getNumCols());
  solver = const_cast<SolverInterface*>(dynamic_cast<const SolverInterface*>(&si));

  SECTION( "Test without strong branching or pruned terms" ) {

    // make disjunction
    PartialBBDisjunction disj = PartialBBDisjunction(vpc_params);
    disj.prepareDisjunction(solver);

    // perturb the problem
    for (int col_idx=0; col_idx < si.getNumCols(); col_idx++){
      bool perturb = col_idx == 25 || col_idx == 0 || col_idx == 13 || col_idx == 22;
      si.setObjCoeff(
          col_idx, perturb ? si.getObjCoefficients()[col_idx] * 5 :
          si.getObjCoefficients()[col_idx]);
    }

    // update the constraint bounds
    for (int row_idx=0; row_idx < si.getNumRows(); row_idx++){
      double b = si.getRowUpper()[row_idx];
      double looser_b = b > 0 ? b * 1.5 : b / 1.5;
      bool perturb = row_idx == 2 || row_idx == 7 || row_idx == 12 || row_idx == 17;
      si.setRowUpper(row_idx, perturb ? looser_b: b);
    }
    si.resolve();

    PartialBBDisjunction param_disj = disj.parameterize(&si);

    // check that we have the right metadata
    REQUIRE(param_disj.num_terms == 4);
    REQUIRE(param_disj.terms.size() == 4);
    REQUIRE(param_disj.common_changed_var.size() == 0);
    REQUIRE(param_disj.common_changed_bound.size() == 0);
    REQUIRE(param_disj.common_changed_value.size() == 0);
    REQUIRE(param_disj.common_ineqs.size() == 0);
    REQUIRE(param_disj.integer_sol.size() == 0);
    REQUIRE(param_disj.integer_obj > 1e300);
    REQUIRE(param_disj.root_obj == si.getObjValue());
    REQUIRE(isVal(param_disj.best_obj, 19.33, .1));
    REQUIRE(isVal(param_disj.worst_obj, 49.17, .1));

    for (int term_idx = 0; term_idx < 4; term_idx++){

      // these attributes should be the same
      REQUIRE(param_disj.terms[term_idx].type == disj.terms[term_idx].type);
      REQUIRE(param_disj.terms[term_idx].changed_var == disj.terms[term_idx].changed_var);
      REQUIRE(param_disj.terms[term_idx].changed_bound == disj.terms[term_idx].changed_bound);
      REQUIRE(param_disj.terms[term_idx].changed_value == disj.terms[term_idx].changed_value);
      REQUIRE(param_disj.terms[term_idx].is_feasible);

      // these attributes should be different
      REQUIRE(param_disj.terms[term_idx].obj != disj.terms[term_idx].obj);
      REQUIRE(param_disj.terms[term_idx].basis != disj.terms[term_idx].basis);
    }
  }

}


