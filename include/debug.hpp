/**
 * debug.hpp
 * A. M. Kazachkov
 * 2018-Dec-25
 */

#include <string>

#include "params.hpp"
#include "VPCEventHandler.hpp"

std::string generateTreePlotString(const VPCEventHandler* eventHandler, const VPCParameters& params,
    const bool saveToFile = false);
std::string generateTikzTreeString(const VPCEventHandler* eventHandler, const VPCParameters& params,
    const int orig_strategy, const double branching_lb, const bool saveToFile);
