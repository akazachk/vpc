/**
 * debug.cpp
 * A. M. Kazachkov
 * 2018-Dec-25
 */
#include "debug.hpp"
#include <vector>

#include "utility.hpp"
#include "VPCParameters.hpp"
#include "VPCEventHandler.hpp"
#include "SolverHelper.hpp"
#include "PartialBBDisjunction.hpp" // for generationPartialTree

void printVector(const int n, const double* vec) {
  for (int i = 0; i < n; ++i) {
    if (vec[i] != 0.0)
      printf("\t[%d] = %e\n", i, vec[i]);
    else
      printf("\t[%d] = 0\n", i);
  }
}

#ifdef USE_CBC
#include <CbcModel.hpp>
void printTree(PartialBBDisjunction* const orig_owner,
    OsiSolverInterface* solver, OsiCuts* vpcs, OsiCuts* gmics) {
  const int TEMP_VAL= orig_owner->params.get(intParam::TEMP);
  if (std::abs(TEMP_VAL) >= static_cast<int>(TempOptions::GEN_TIKZ_STRING_WITH_VPCS)
      && std::abs(TEMP_VAL) <= static_cast<int>(TempOptions::GEN_TIKZ_STRING_NO_CUTS)) {
    const bool shouldApplyVPCs =
        (vpcs != NULL)
            && (TEMP_VAL == static_cast<int>(TempOptions::GEN_TIKZ_STRING_WITH_VPCS)
                || TEMP_VAL == static_cast<int>(TempOptions::GEN_TIKZ_STRING_WITH_VPCS_AND_GMICS));
    const bool shouldApplySICs =
        (gmics != NULL)
            && (TEMP_VAL == static_cast<int>(TempOptions::GEN_TIKZ_STRING_WITH_GMICS)
                || TEMP_VAL == static_cast<int>(TempOptions::GEN_TIKZ_STRING_WITH_VPCS_AND_GMICS));

    VPCParameters tmp_params = orig_owner->params;
    const int old_param_val = tmp_params.get(PARTIAL_BB_STRATEGY);
    int num_strong = tmp_params.get(intParam::PARTIAL_BB_NUM_STRONG);
    if (num_strong == -1) {
      num_strong = solver->getNumCols();
    }
    if (num_strong == -2) {
      num_strong = static_cast<int>(std::ceil(std::sqrt(solver->getNumCols())));
    }

    //const std::vector<int> testset = {000, 001, 004, 010, 011, 014, 020, 021, 024, 030, 031, 034};
//    const std::vector<int> test_num_trust = {-1, std::numeric_limits<int>::max()};
//    const std::vector<int> testset = {000, 010, 030, 004, 014, 034};
//      const std::vector<int> test_paramset = { 104, 204, 304, 404, 504, -104, -204, -304, -404, -504, -604 };
    const std::vector<int> test_num_trust = { std::numeric_limits<int>::max() };
    const std::vector<int> test_paramset = { 000 };
    for (int num_before_trusted : test_num_trust) {
      for (int test_param : test_paramset) {
        tmp_params.set(PARTIAL_BB_STRATEGY, test_param);
        SolverInterface* newBBSolver =
            dynamic_cast<SolverInterface*>(solver->clone());
        newBBSolver->disableFactorization();
        if (shouldApplyVPCs) {
          newBBSolver->applyCuts(*vpcs);
        }
        if (shouldApplySICs) {
          newBBSolver->applyCuts(*gmics);
        }
        newBBSolver->resolve();
        printf(
            "\n## Generating partial branch-and-bound tree after adding cuts (VPCs: %d, SICs: %d, total new num rows: %d, new obj value: %.3f). ##\n",
            shouldApplyVPCs, shouldApplySICs,
            newBBSolver->getNumRows() - solver->getNumRows(),
            newBBSolver->getObjValue());
        setupClpForCbc(newBBSolver);
        CbcModel* new_cbc_model = new CbcModel;
        new_cbc_model->swapSolver(newBBSolver);
        new_cbc_model->setModelOwnsSolver(true); // solver will be deleted with cbc object
        setIPSolverParameters(new_cbc_model);

        PartialBBDisjunction* owner = orig_owner->clone();
        owner->params = tmp_params;
        generatePartialBBTree(owner, new_cbc_model, newBBSolver,
            std::numeric_limits<int>::max(), num_strong, num_before_trusted);

        VPCEventHandler* newEventHandler;
        try {
          newEventHandler =
              dynamic_cast<VPCEventHandler*>(new_cbc_model->getEventHandler());
        } catch (std::exception& e) {
          error_msg(errstr, "Could not get event handler.\n");
          writeErrorToLog(errstr, owner->params.logfile);
          exit(1);
        }
        printNodeStatistics(newEventHandler->getStatsVector(), false);
        if (newEventHandler->getPrunedStatsVector().size() > 0) {
          printf("\n");
          printNodeStatistics(newEventHandler->getPrunedStatsVector(), false);
        }
        generateTikzTreeString(newEventHandler, tmp_params, old_param_val,
            owner->best_obj, true);
        if (new_cbc_model->getNodeCount() < 70) {
          printf(
              "Number nodes: %d.\nCut strategy: %d.\n# leafs: %d.\nSolve strategy: %d.\nDo you wish to exit? (y/n) ",
              new_cbc_model->getNodeCount(), old_param_val, owner->num_terms,
              test_param);
          std::string response;
          std::getline(std::cin, response);
          if (response.compare("y") == 0 || response.compare("yes") == 0) {
            exit(1);
          }
        }
        if (new_cbc_model) {
          delete new_cbc_model;
          new_cbc_model = NULL;
        }
        if (owner) { delete owner; }

        tmp_params.set(PARTIAL_BB_STRATEGY, old_param_val);
      } // loop over param test set
    } // loop over num trusted to be tested
  } // check temp for whether to apply cuts and print post-cuts tree
} /* printTree */
#endif // USE_CBC

/************************************************************/
/**
 * Create a string that can be fed into Mathematica to plot the tree
 */
std::string generateTreePlotString(const VPCEventHandler* eventHandler, const VPCParameters& params,
    const bool saveToFile) {
  const std::vector<NodeStatistics>& stats = eventHandler->getStatsVector();
  const std::vector<NodeStatistics>& pruned_stats = eventHandler->getPrunedStatsVector();
  int max_node_ind = eventHandler->getNumNodes() - 1;

  std::string prepend = "TreePlot[{";
  std::string append = "}, Automatic, 0, VertexLabeling -> True]";
  std::string str;

  // First handle the nodes already explored
  for (int i = 0; i < (int) stats.size(); i++) {
    if (stats[i].orig_id != stats[i].id || stats[i].parent_id < 0) {
      continue;
    }
    if (!str.empty()) {
      str += ", ";
    }
    str += "{";
    str += std::to_string(stats[stats[i].parent_id].number);
    str += "->";
    str += std::to_string(stats[i].number);
    str += ", \"" + std::to_string(stats[stats[i].parent_id].variable) + "\"}";
  }

  // Add the infeasible nodes
  for (int i = 0; i < (int) pruned_stats.size(); i++) {
    if (!str.empty()) {
      str += ", ";
    }
    str += "{";
    str += std::to_string(stats[pruned_stats[i].parent_id].number);
    str += "->";
    str += std::to_string(pruned_stats[i].number) + "X";
    str += ", \"" + std::to_string(stats[pruned_stats[i].parent_id].variable)
        + "\"}";
  }

  // Now add the leaf nodes
  for (int tmp_ind = 0; tmp_ind < eventHandler->getNumNodesOnTree(); tmp_ind++) {
    const int i = eventHandler->getNodeIndex(tmp_ind);
    const int branch = stats[i].branch_index;
    for (int b = branch; b < 2; b++) {
      max_node_ind++;
      if (!str.empty()) {
        str += ", ";
      }
      str += "{";
      str += std::to_string(stats[i].number);
      str += "->";
      str += std::to_string(max_node_ind);
      str += ", \"" + std::to_string(stats[i].variable) + "\"}";
    }
  }

#ifdef TRACE
  printf("%s%s%s\n", prepend.c_str(), str.c_str(), append.c_str());
#endif

  if (saveToFile) {
    std::string filename = params.get(stringParam::FILENAME) + "-Tree.alex";
    FILE* myfile = fopen(filename.c_str(), "a");
    fprintf(myfile, "%s%s%s\n", prepend.c_str(), str.c_str(), append.c_str());
    printNodeStatistics(stats, false, myfile);
    printNodeStatistics(pruned_stats, false, myfile);
    fprintf(myfile, "\n");
    fclose(myfile);
  }

  return prepend + str + append;
} /* generateTreePlotString */

/************************************************************/
/**
 * Helper for generateTikzTreeString
 *
 * loc = 0: normal node
 * loc = 1: node pruned by integrality
 * loc = 2: node pruned by infeasibility or bound?
 * loc = 3: unexplored node
 * loc = 4: integer solution found during strong branching
 */
std::string generateTikzTreeStringHelper(const int node_ind, const int depth,
    const int numNodes, const std::vector<std::vector<int> >& locationOfChild,
    const std::vector<std::vector<int> >& children,
    const std::vector<int>& variable,
    const std::vector<double>& obj, const bool print_obj,
    const int child_ind = 0, const int loc = 0) {
  if (node_ind < 0) {
   return "";
  }
  const int curr_depth = 2 * (depth + 1);
  const int num_children = (node_ind < numNodes && loc != 4) ? (int) locationOfChild[node_ind].size() : 0;
  std::string str = std::string(curr_depth, ' ');
  str += std::to_string(node_ind);
  if (print_obj && node_ind < static_cast<int>(obj.size())) {
    str += "/\"" + stringValue(obj[node_ind], "%.2f") + "\"";
  }

  // Add information to print at node
  if (loc > 0 || node_ind > 0) {
    str += " [";
  }
  if (loc == 1 || loc == 4) {
    str += "green!75!black";
  } else if (loc == 2) {
    str += "red";
  } else if (loc == 3) {
    str += "blue";
  }
  if (node_ind > 0) {
    if (loc > 0) {
      str += ", ";
    }
    if (loc != 4) {
      str += ">\"{";
      str += std::to_string(variable[node_ind]);
      str += "}\"";
      if (child_ind == 0) {
        str += ", >swap";
      }
    } else {
      str += ">dashed";
    }
  }
  if (loc > 0 || node_ind > 0) {
    str += "]";
  }

  // Continue recursive process
  if (node_ind < numNodes) {
    bool addedString = false;
    for (int j = 0; j < num_children; j++) {
      if (children[node_ind][j] < 0)
        continue;
      const int loc = locationOfChild[node_ind][j];
      if (loc < 0) {
        continue;
      }

      if (!addedString) {
        addedString = true;
        str += " -- {\n";
      } else {
        str += ",\n";
      }
      str += generateTikzTreeStringHelper(children[node_ind][j], depth + 1,
          numNodes, locationOfChild, children, variable, obj, print_obj, j, loc);
    }
    if (addedString) {
      str += "\n";
      str += std::string(curr_depth, ' ') + "}";
    }
  }
  if (node_ind == 0) {
    str += "\n";
  }
  return str;
} /* generateTikzTreeStringHelper */

void setChildForTikzTreeString(std::vector<std::vector<int> >& children,
    std::vector<std::vector<int> >& locationOfChild, const int loc,
    const int parent_node_ind, const int child_ind, const int child_node_ind) {
  if (locationOfChild[parent_node_ind][child_ind] >= 0) {
    // Resize because that child already exists
    const int size = locationOfChild[parent_node_ind].size();
    locationOfChild[parent_node_ind].resize(size + 1);
    children[parent_node_ind].resize(size + 1);
    children[parent_node_ind][size] = child_node_ind;
    locationOfChild[parent_node_ind][size] = loc;
  } else {
    children[parent_node_ind][child_ind] = child_node_ind;
    locationOfChild[parent_node_ind][child_ind] = loc;
  }
} /* setChildForTikzTreeString */

/**
 * Create a string that can be fed into LuaLaTeX to plot the tree
 */
std::string generateTikzTreeString(const VPCEventHandler* eventHandler,
    const VPCParameters& params,
    const int orig_strategy, const double branching_lb, const bool saveToFile) {
  const std::vector<NodeStatistics>& stats = eventHandler->getStatsVector();
  const std::vector<NodeStatistics>& pruned_stats = eventHandler->getPrunedStatsVector();
  const int numNodes = eventHandler->getNumNodes();
  int numNodesTotal = numNodes;
  int numSBFeasNodes = 0;

  std::string prepend = "\\begin{figure}[htp!]\n";
  prepend += "\\centering\n";
  prepend += "\\makebox[\\textwidth][c]{\n";
  prepend += "\\begin{tikzpicture} [tree layout, scale = 0.5, font=\\tiny\\bfseries]\n";
  prepend += "\\graph [\n";
  prepend += "\%  nodes={inner sep=0pt, minimum size=2.5pt, circle, fill, draw}, empty nodes,\n";
  prepend += "  edges={style={pos=0.5, inner sep=0.5pt, font=\\tiny}}\n";
  prepend += "] {\n";
  std::string append = "};\n\\end{tikzpicture}\n}\n";
  std::string endfigure = "\\end{figure}\n";
//  std::vector<int> depth(numNodes);
  std::vector<std::vector<int> > children(numNodes); // the children of this node (node indices)
  std::vector<std::vector<int> > locationOfChild(numNodes); // 0: children, 1: feas_children, 2: pruned_children, 3: unexplored_children, 4: feas children found during strong branching
  std::vector<int> variable(numNodes); // variable that was branched on to get to this node
  std::vector<double> obj(numNodes); // objective value at this node
  variable[0] = -1;
  obj[0] = eventHandler->getOriginalSolver()->getObjValue();

  // Assume there are at most 3 children
  // (middle "child" reserved for feas sol found during strong branching at that node)
  // There may be more / less in general
  for (int i = 0; i < numNodes; i++) {
    children[i].resize(3, -1);
    locationOfChild[i].resize(3, -1);
  }

  // First handle the nodes already explored
  for (int i = 0; i < (int) stats.size(); i++) {
    if (stats[i].orig_id != stats[i].id || stats[i].parent_id < 0) {
      continue; // skip if the node has already been handled
    }
    const int node_ind = stats[i].number;
    const int parent_id = stats[i].parent_id;
    const int parent_node_ind = stats[parent_id].number;
    const int child_ind = (stats[parent_id].way <= 0) ? 0 : 2;
    int loc = 0;
    if (stats[i].found_integer_solution && isVal(stats[i].obj, stats[i].integer_obj)) {
      loc = 1;
    }
    setChildForTikzTreeString(children, locationOfChild, loc,
        parent_node_ind, child_ind, node_ind);
    variable[node_ind] = stats[parent_id].variable;
    obj[node_ind] = stats[i].obj;
    if (stats[i].found_integer_solution && !isVal(stats[i].obj, stats[i].integer_obj)) {
      // An integer solution was found during strong branching
      // Add this as a child with loc set to 4
      children.resize(numNodes+numSBFeasNodes+1);
      locationOfChild.resize(numNodes+numSBFeasNodes+1);
      children[numNodes+numSBFeasNodes].resize(3,-1);
      locationOfChild[numNodes+numSBFeasNodes].resize(3,-1);
      loc = 4;
      setChildForTikzTreeString(children, locationOfChild, loc,
          node_ind, 1, numNodes+numSBFeasNodes); // place it in the middle
      variable.push_back(-1);
      obj.push_back(stats[i].integer_obj);
      numSBFeasNodes++;
    }
  }

  // Add the leaf nodes
  int numActiveNodes = eventHandler->getNumNodesOnTree();
  for (int tmp_ind = 0; tmp_ind < numActiveNodes; tmp_ind++) {
    const int i = eventHandler->getNodeIndex(tmp_ind);
    const int branch = stats[i].branch_index;
    const int first_child_ind = (stats[i].way <= 0);
    const int parent_node_ind = stats[i].number;
    for (int b = branch; b < 2; b++) {
      const int tmp_child_ind = (b == 0) ? first_child_ind : !first_child_ind;
      const int child_ind = (tmp_child_ind == 0) ? 0 : 2; // reserve index 1 for strong branching "child"
      setChildForTikzTreeString(children, locationOfChild, 3, parent_node_ind, child_ind, numNodesTotal);
      variable.push_back(stats[i].variable);
      obj.push_back(-1);
      numNodesTotal++;
    }
  }

  // Add the pruned nodes
  int numInfeasNodes = 0, numFeasNodes = 0;
  for (int i = 0; i < (int) pruned_stats.size(); i++) {
    const int node_ind = pruned_stats[i].number;
    assert(children[node_ind][0] == -1 && children[node_ind][1] == -1);
    const int parent_id = pruned_stats[i].parent_id;
    const int parent_node_ind = stats[parent_id].number;
    const int child_ind = (stats[parent_id].way <= 0) ? 0 : 2;
    variable[node_ind] = stats[parent_id].variable;
    obj[node_ind] = pruned_stats[i].obj;
    if (pruned_stats[i].found_integer_solution) {
      int loc = (pruned_stats[i].obj == pruned_stats[i].integer_obj);
      setChildForTikzTreeString(children, locationOfChild, loc, parent_node_ind, child_ind, node_ind);
      numFeasNodes++;
      if (loc == 0) {
        // An integer solution was found during strong branching
        // Add this as a child with loc set to 4
        children.resize(numNodes+numSBFeasNodes+1);
        locationOfChild.resize(numNodes+numSBFeasNodes+1);
        children[numNodes+numSBFeasNodes].resize(3,-1);
        locationOfChild[numNodes+numSBFeasNodes].resize(3,-1);
        loc = 4;
        setChildForTikzTreeString(children, locationOfChild, loc,
            node_ind, 1, numNodes+numSBFeasNodes); // place it in the middle
        variable.push_back(-1);
        obj.push_back(pruned_stats[i].integer_obj);
        numSBFeasNodes++;
      }
    } else { // although it says infeasible, it may have been pruned by bound...
      setChildForTikzTreeString(children, locationOfChild, 2, parent_node_ind, child_ind, node_ind);
      numInfeasNodes++;
    }
  }

  // In the case that there are no leaf nodes, move the dangling ones to pruned
  if (numActiveNodes == 0) {
    std::vector<bool> isLeafNode(numNodes, false);
    for (int i = numNodes - 1; i >= 0; i--) {
      int num_children = locationOfChild[i].size();
      int num_empty = 0;
      for (int j = 0; j < num_children; j++) {
        const int loc = locationOfChild[i][j];
        if (loc < 0) {
          num_empty++;
        } else if (loc == 0) {
          const int child_node_ind = children[i][j];
          if (isLeafNode[child_node_ind]) {
            locationOfChild[i][j] = 2;
          }
        }
      }
      if (num_children - num_empty == 0) {
        isLeafNode[i] = true;
      }
    }
  }

  // Finally make the string
  const bool print_obj = params.get(intParam::TEMP) >= 0;
  std::string graph = generateTikzTreeStringHelper(0, 0, numNodes,
      locationOfChild, children, variable, obj, print_obj, 0,
      stats[0].found_integer_solution);

  // Add a caption
  int numLeafNodes = eventHandler->getNumLeafNodes();
  if (numLeafNodes == 0) {
    numLeafNodes = numNodesTotal - numNodes + numFeasNodes; // this may differ from the real num leaf nodes due to pruning we do at the end based on the best solution
  }
  std::string caption = "\\caption{\n";
  caption += "  \\detokenize{" + params.get(stringParam::FILENAME) + "}:";
  if (numActiveNodes > 0) {
    caption += " Partial tree";
    caption +=
        " (used for cut generation; generated with strategy "
            + stringValue(orig_strategy, "%03d") + " and "
            + std::to_string(-1 * params.get(intParam::DISJ_TERMS))
            + " leafs)";
  } else {
    caption += " Optimal tree";
    caption += " from strategy "
        + stringValue(params.get(intParam::PARTIAL_BB_STRATEGY), "%03d");
    caption +=
         " (partial tree: strategy "
            + stringValue(orig_strategy, "%03d") + " and "
            + std::to_string(-1 * params.get(intParam::DISJ_TERMS))
            + " leafs).";
  }
  caption += + "\\\\\n";
  caption += "  Initial objective value: "
      + stringValue(eventHandler->getOriginalSolver()->getObjValue(), "%.3f") + ".";
  if (numActiveNodes == 0) {
    caption += " Branching lower bound: " + stringValue(branching_lb, "%.3f") + ".";
  }
  caption += + "\\\\\n";
  caption += "  Nodes: " + std::to_string(numNodesTotal) + ".";
  caption += "   Active nodes: " + std::to_string(numActiveNodes) + ".";
  caption += "   Feasible leafs: " + std::to_string(numLeafNodes) + ".";
  caption += "   Integer-feasible leafs: " + std::to_string(numFeasNodes) + ".";
  caption += "   Pruned nodes: " + std::to_string(numInfeasNodes) + ".\n";
  caption += "}\n";

  // Collect all the text in one string
  std::string str = prepend + graph + append + caption + endfigure;

#ifdef TRACE
  printf("%s\n", str.c_str());
#endif

  if (saveToFile) {
    std::string dir, instname, ext;
    parseFilename(dir, instname, ext, params);
    std::string filename = dir + "/" + instname + "-Tree.alex";
    FILE* myfile = fopen(filename.c_str(), "a");
    if (!myfile) {
      error_msg(errorstring, "Failed to open %s.\n", filename.c_str());
      writeErrorToLog(errorstring, params.logfile);
      exit(1);
    }
    fprintf(myfile, "## Tree for %s with vpc_depth = %d ##\n",
        params.get(stringParam::FILENAME).c_str(),
        params.get(intParam::DISJ_TERMS));
    fprintf(myfile, "%s\n", str.c_str());
    printNodeStatistics(stats, false, myfile);
    printNodeStatistics(pruned_stats, false, myfile);
    fprintf(myfile, "\n");
    fclose(myfile);

    filename = dir + "/" + instname + ".tex";
    myfile = fopen(filename.c_str(), "w");
    if (!myfile) {
      error_msg(errorstring, "Failed to open %s.\n", filename.c_str());
      writeErrorToLog(errorstring, params.logfile);
      exit(1);
    }
    fprintf(myfile, "%% !TEX TS-program = lualatex\n");
    fprintf(myfile, "\\documentclass[11pt,landscape]{article}\n");
    fprintf(myfile, "\\usepackage[margin=1in]{geometry}\n");
    fprintf(myfile, "\\usepackage[T1]{fontenc}\n");
    fprintf(myfile, "\\usepackage{caption}\n");
    fprintf(myfile, "\\usepackage{tikz}\n");
    fprintf(myfile, "\\usetikzlibrary{graphs, positioning, fit, quotes}\n");
    fprintf(myfile, "\\usetikzlibrary{graphdrawing}\n");
    fprintf(myfile, "\\usegdlibrary{trees}\n");
    fprintf(myfile, "\\begin{document}\n");
    fprintf(myfile, "%s\n", str.c_str());
    fprintf(myfile, "\\end{document}\n");
    fclose(myfile);
  }

  return str;
} /* generateTikzTreeString */

///**
// * @brief Checks if a point computed in the complemented nonbasic space is equivalent to the original point in the structural space
// */
//void checkPoint(CoinPackedVector point, OsiSolverInterface* origSolver, const double* struct_point) {
//  const int num_cols = origSolver->getNumCols();
//  const int num_rows = origSolver->getNumRows();
//
//  // Get nonbasic variables
//  std::vector<int> NBVarIndex;
//  for (int var = 0; var < num_cols+num_rows; var++) {
//    if (!isBasicVar(origSolver, var)) {
//
//    }
//  }
//
//  // Check objective value matches up
//  double struct_obj = 0.;
//  for (int i = 0; i < num_cols; i++) {
//    struct_obj += struct_point[i] * origSolver->getObjCoefficients()[i];
//  }
//
////  double nb_obj = 0.;
//
//} /* checkPoint */
