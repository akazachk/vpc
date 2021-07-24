/**
 * @file vpc_debug.hpp
 * @brief Useful functions for debugging that are not necessary in the rest of the code
 *
 * @author A. M. Kazachkov
 * @date 2018-12-25
 */
#include <string>
#include <vector>

namespace VPCParametersNamespace {
  struct VPCParameters;
}
class VPCEventHandler;
class PartialBBDisjunction;
class OsiSolverInterface;
class OsiCuts;

#ifdef USE_CBC
/// @brief Save tikz tree from \p orig_owner
void printTree(PartialBBDisjunction* const orig_owner,
    OsiSolverInterface* solver, OsiCuts* vpcs, OsiCuts* gmics);
#endif

/// @brief Create a string that can be fed into Mathematica to plot the tree
std::string generateTreePlotString(const VPCEventHandler* eventHandler, const VPCParametersNamespace::VPCParameters& params,
    const bool saveToFile = false);

/// @brief Create a string that can be fed into LuaLaTeX to plot the tree
std::string generateTikzTreeString(const VPCEventHandler* eventHandler, const VPCParametersNamespace::VPCParameters& params,
    const int orig_strategy, const double branching_lb, const std::string filename_stub = "");
