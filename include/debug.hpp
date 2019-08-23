/**
 * debug.hpp
 * A. M. Kazachkov
 * 2018-Dec-25
 */

#include <string>

struct VPCParameters;
class VPCEventHandler;
class PartialBBDisjunction;
class OsiSolverInterface;
class OsiCuts;
class CoinPackedVector;

void printVector(const CoinPackedVector& vec);
template <typename T>
void printVector(const int n, const T* vec);

void printTree(PartialBBDisjunction* const orig_owner,
    OsiSolverInterface* solver, OsiCuts* vpcs, OsiCuts* gmics);
std::string generateTreePlotString(const VPCEventHandler* eventHandler, const VPCParameters& params,
    const bool saveToFile = false);
std::string generateTikzTreeString(const VPCEventHandler* eventHandler, const VPCParameters& params,
    const int orig_strategy, const double branching_lb, const bool saveToFile);
