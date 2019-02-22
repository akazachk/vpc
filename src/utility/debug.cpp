/**
 * debug.cpp
 * A. M. Kazachkov
 * 2018-Dec-25
 */
#include "debug.hpp"
#include <vector>

#include "utility.hpp"
#include "SolverHelper.hpp"

/**
 * @brief Checks if a point computed in the complemented nonbasic space is equivalent to the original point in the structural space
 */
void checkPoint(CoinPackedVector point, OsiSolverInterface* origSolver, const double* struct_point) {
  const int num_cols = origSolver->getNumCols();
  const int num_rows = origSolver->getNumRows();

  // Get nonbasic variables
  std::vector<int> NBVarIndex;
  for (int var = 0; var < num_cols+num_rows; var++) {
    if (!isBasicVar(origSolver, var)) {

    }
  }

  // Check objective value matches up
  double struct_obj = 0.;
  for (int i = 0; i < num_cols; i++) {
    struct_obj += struct_point[i] * origSolver->getObjCoefficients()[i];
  }

//  double nb_obj = 0.;

} /* checkPoint */

void printVector(const int n, const double* vec) {
  for (int i = 0; i < n; ++i) {
    if (vec[i] != 0.0)
      printf("\t[%d] = %e\n", i, vec[i]);
    else
      printf("\t[%d] = 0\n", i);
  }
}

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
  const int num_children = (node_ind < numNodes) ? (int) locationOfChild[node_ind].size() : 0;
  std::string str = std::string(curr_depth, ' ');
  str += std::to_string(node_ind);
  if (print_obj && num_children > 0) {
    str += "/\"" + stringValue(obj[node_ind], "%.2f") + "\"";
  }

  // Add information to print at node
  if (loc > 0 || node_ind > 0) {
    str += " [";
  }
  if (loc == 1) {
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
    str += ">\"{";
    str += std::to_string(variable[node_ind]);
    str += "}\"";
    if (child_ind == 0) {
      str += ", >swap";
    }
  }
  if (loc > 0 || node_ind > 0) {
    str += "]";
  }

  // Continue recursive process
  if (node_ind < numNodes) {
    bool addedString = false;
    for (int j = 0; j < num_children; j++) {
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
  std::vector<std::vector<int> > children(numNodes); // the children of this node
  std::vector<std::vector<int> > locationOfChild(numNodes); // 0: children, 1: feas_children, 2: pruned_children, 3: unexplored_children
  std::vector<int> variable(numNodes); // variable that was branched on to get to this node
  std::vector<double> obj(numNodes); // variable that was branched on to get to this node
  variable[0] = -1;
  obj[0] = eventHandler->getOriginalSolver()->getObjValue();

  // Assume there are two children
  // There may be more in the case of infeasibility
  for (int i = 0; i < numNodes; i++) {
    children[i].resize(2, -1);
    locationOfChild[i].resize(2, -1);
  }

  // First handle the nodes already explored
  for (int i = 0; i < (int) stats.size(); i++) {
    if (stats[i].orig_id != stats[i].id || stats[i].parent_id < 0) {
      continue;
    }
    const int node_ind = stats[i].number;
    const int parent_id = stats[i].parent_id;
    const int parent_node_ind = stats[parent_id].number;
    const int child_ind = (stats[parent_id].way <= 0) ? 0 : 1;
    setChildForTikzTreeString(children, locationOfChild, stats[i].found_integer_solution,
        parent_node_ind, child_ind, node_ind);
    variable[node_ind] = stats[parent_id].variable;
    obj[node_ind] = stats[i].obj;
  }

  // Add the leaf nodes
  int numActiveNodes = eventHandler->getNumNodesOnTree();
  for (int tmp_ind = 0; tmp_ind < numActiveNodes; tmp_ind++) {
    const int i = eventHandler->getNodeIndex(tmp_ind);
    const int branch = stats[i].branch_index;
    const int first_child_ind = (stats[i].way <= 0);
    const int parent_node_ind = stats[i].number;
    for (int b = branch; b < 2; b++) {
      const int child_ind = (b == 0) ? first_child_ind : !first_child_ind;
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
    const int child_ind = (stats[parent_id].way <= 0) ? 0 : 1;
    variable[node_ind] = stats[parent_id].variable;
    obj[node_ind] = pruned_stats[i].obj;
    if (pruned_stats[i].found_integer_solution) {
      setChildForTikzTreeString(children, locationOfChild, 1, parent_node_ind, child_ind, node_ind);
      numFeasNodes++;
    } else {
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
    std::string filename = params.get(stringParam::FILENAME) + "-Tree.alex";
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

    filename = params.get(stringParam::FILENAME) + ".tex";
    myfile = fopen(filename.c_str(), "w");
    if (!myfile) {
      error_msg(errorstring, "Failed to open %s.\n", filename.c_str());
      writeErrorToLog(errorstring, params.logfile);
      exit(1);
    }
    fprintf(myfile, "%% !TEX TS-program = lualatex\n");
    fprintf(myfile, "\\documentclass[11pt,landscape]{article}\n");
    fprintf(myfile, "\\usepackage{akazachk}\n");
    fprintf(myfile, "\\begin{document}\n");
    fprintf(myfile, "%s\n", str.c_str());
    fprintf(myfile, "\\end{document}\n");
    fclose(myfile);
  }

  return str;
} /* generateTikzTreeString */
