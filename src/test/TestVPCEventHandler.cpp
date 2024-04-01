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

  // solver
  OsiClpSolverInterface si;
  SolverInterface* solver;
  si.readMps("../test/p0033.mps");
  si.initialSolve();
  solver = const_cast<SolverInterface*>(dynamic_cast<const SolverInterface*>(&si));

  SECTION( "Test without strong branching or pruned terms" ) {

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

    // make sure to save the full tree and update to keep same active leaves at end
    vpc_params.set(VPCParametersNamespace::PARTIAL_BB_KEEP_PRUNED_NODES, 1);
    vpc_params.set(VPCParametersNamespace::DISJ_TERMS, 7);

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
        {14, 10, 0, 30, 9}, {14, 10, 0, 30, 9}
    };
    std::vector<std::vector<int> > bound = {
        {0}, {1, 0, 1}, {1, 0, 0}, {1, 1, 1}, {1, 1, 0, 0}, {1, 1, 0, 1, 1}, {1, 1, 0, 1, 0}
    };
    std::vector<std::vector<double> > value = {
        {1}, {0, 1, 0}, {0, 1, 1}, {0, 0, 0}, {0, 0, 1, 1}, {0, 0, 1, 0, 0}, {0, 0, 1, 0, 1}
    };
    std::vector<std::string> term_type = {
        "complement", "leaf", "leaf", "complement", "complement", "leaf", "leaf"
    };
    for (int i = 0; i < disj.num_terms; i++) {
      REQUIRE(disj.terms[i].changed_var == var[i]);
      REQUIRE(disj.terms[i].changed_bound == bound[i]);
      REQUIRE(disj.terms[i].changed_value == value[i]);
      REQUIRE(disj.terms[i].type == term_type[i]);
    }

    // check misc
    REQUIRE(disj.common_ineqs.size() == 0);
    REQUIRE(disj.integer_sol.size() == 0);
    REQUIRE(disj.integer_obj > 1e300);
  }

  // load a different problem to test pruned terms
  si.readMps("../test/bm23.mps");
  si.initialSolve();
  solver = const_cast<SolverInterface*>(dynamic_cast<const SolverInterface*>(&si));

  // create 64 terms so we get a couple of pruned ones
  vpc_params.set(VPCParametersNamespace::DISJ_TERMS, 69);

  SECTION( "Test with pruned terms and strong branching terms" ) {

    // make sure to save the full tree
    vpc_params.set(VPCParametersNamespace::PARTIAL_BB_KEEP_PRUNED_NODES, 1);

    // make disjunction
    PartialBBDisjunction disj = PartialBBDisjunction(vpc_params);
    disj.prepareDisjunction(solver);

    std::vector<int> pruned, complement, leaf;
    for (int node_idx = 0; node_idx < disj.num_terms; node_idx++) {
      if (disj.terms[node_idx].type == "pruned") {
        pruned.push_back(node_idx);
      } else if (disj.terms[node_idx].type == "leaf") {
        leaf.push_back(node_idx);
      } else if (disj.terms[node_idx].type == "complement") {
        complement.push_back(node_idx);
      } else {
        // it should be one of the three above
        REQUIRE(false);
      }
    }
    // should have 2 pruned terms and 64 leaves
    REQUIRE(pruned.size() == 2);
    REQUIRE(complement.size() == 3);
    REQUIRE(leaf.size() == 64);

    // check that we have no common terms anymore
    REQUIRE(disj.common_changed_var.size() == 0);
    REQUIRE(disj.common_changed_bound.size() == 0);
    REQUIRE(disj.common_changed_value.size() == 0);

    // check the two pruned terms come through
    std::vector<std::vector<int> > var = {
        {16, 13, 10, 7, 12, 21, 17, 22, 11, 15},
        {16, 13, 10, 7, 12, 21, 17, 22, 11, 8}
    };
    std::vector<std::vector<int> > bound = {
        {0, 1, 1, 1, 1, 1, 1, 1, 1, 1},
        {0, 1, 1, 1, 1, 1, 1, 1, 0, 1}
    };
    std::vector<std::vector<double> > value = {
        {1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {1, 0, 0, 0, 0, 0, 0, 0, 1, 0}
    };
    for (int i = 0; i < pruned.size(); i++) {
      REQUIRE(disj.terms[pruned[i]].changed_var == var[i]);
      REQUIRE(disj.terms[pruned[i]].changed_bound == bound[i]);
      REQUIRE(disj.terms[pruned[i]].changed_value == value[i]);
    }

    // make sure we have no same terms
    for (int i = 0; i < disj.terms.size(); i++) {
      for (int j = i+1; j < disj.terms.size(); j++) {
        // if the same variables are changed, make sure they had different bounds or values
        if (disj.terms[i].changed_var == disj.terms[j].changed_var){
          REQUIRE((disj.terms[i].changed_bound != disj.terms[j].changed_bound ||
                   disj.terms[i].changed_value != disj.terms[j].changed_value));
        }
      }
    }

    // check misc
    REQUIRE(disj.common_ineqs.size() == 0);
    REQUIRE(disj.integer_sol.size() == 0);
    REQUIRE(disj.integer_obj > 1e300);
  }

  // make sure that we see the same number of cuts and bound improvement
  // regardless of saving the full tree for the same number of total active leaves
  SECTION( "Test integrations" ) {

    // create parameters
    VPCParametersNamespace::VPCParameters partial_disj_params = vpc_params;
    partial_disj_params.set(VPCParametersNamespace::DISJ_TERMS, 64);
    VPCParametersNamespace::VPCParameters full_disj_params = vpc_params;
    full_disj_params.set(VPCParametersNamespace::PARTIAL_BB_KEEP_PRUNED_NODES, 1);

    // create containers for cuts
    OsiCuts partial_tree_vpcs;
    OsiCuts full_tree_vpcs;

    // clone solvers
    SolverInterface* partial_tree_solver =
      dynamic_cast<SolverInterface*>(solver->clone());
    SolverInterface* full_tree_solver =
      dynamic_cast<SolverInterface*>(solver->clone());

    // create cut generators
    CglVPC partial_gen = CglVPC(partial_disj_params);
    CglVPC full_gen = CglVPC(full_disj_params);

    // make cuts
    partial_gen.generateCuts(*partial_tree_solver, partial_tree_vpcs);
    full_gen.generateCuts(*full_tree_solver, full_tree_vpcs);

    // apply cuts
    partial_tree_solver->applyCuts(partial_tree_vpcs);
    full_tree_solver->applyCuts(full_tree_vpcs);

    // resolve
    partial_tree_solver->resolve();
    full_tree_solver->resolve();

    // now make our checks
    // check that the full tree disjunction has more terms, but the same number of active leaves
    REQUIRE(full_gen.disj()->num_terms > partial_gen.disj()->num_terms);
    REQUIRE(full_gen.disj()->num_terms - full_gen.disj()->num_pruned_terms == partial_gen.disj()->num_terms);

    // check that we have the same number of cuts
    REQUIRE(partial_tree_vpcs.sizeCuts() == full_tree_vpcs.sizeCuts());

    // check that we have the same bound improvement
    REQUIRE(isVal(partial_tree_solver->getObjValue(),
                  full_tree_solver->getObjValue(), .1));
  }
}

// --------------------- test current behavior remains -------------------------
TEST_CASE("Test removeContinuousVariableTightenedBounds",
          "[VPCEventHandler::removeContinuousVariableTightenedBounds]") {

  // parameters
  VPCParametersNamespace::VPCParameters vpc_params;
  vpc_params.set(VPCParametersNamespace::DISJ_TERMS, 2);
  vpc_params.set(VPCParametersNamespace::MODE, 0);  // partial BB tree
  vpc_params.set(VPCParametersNamespace::PARTIAL_BB_KEEP_PRUNED_NODES, 1);  // save full tree

  // solver
  OsiClpSolverInterface si;
  SolverInterface* solver;
  si.readMps("../test/dcmulti_7.mps");
  si.initialSolve();
  solver = const_cast<SolverInterface*>(dynamic_cast<const SolverInterface*>(&si));

  // make disjunction
  PartialBBDisjunction disj = PartialBBDisjunction(vpc_params);
  disj.prepareDisjunction(solver);

  // ensure all bounds tightened by the disjunctive term belong to integer variables
  // this would fail without calling removeContinuousVariableTightenedBounds
  for (const DisjunctiveTerm& term : disj.terms) {
    for (const int& v : term.changed_var) {
      REQUIRE(si.isInteger(v));
    }
  }
}