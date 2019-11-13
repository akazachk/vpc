/**
 * @file debug.hpp
 * @brief Useful functions for debugging that are not necessary in the rest of the code
 *
 * @author A. M. Kazachkov
 * @date 2018-12-25
 */
#include <string>
#include <vector>

struct VPCParameters;
class VPCEventHandler;
class PartialBBDisjunction;
class OsiSolverInterface;
class OsiCuts;
class CoinPackedVectorBase;
class CoinPackedVector;

void printVector(const CoinPackedVectorBase& vec, const bool use_newline = true);
void printVectors(const std::vector<CoinPackedVector>& vecs, const bool use_newline = true);
template <typename T>
void printVector(const int n, const T* vec);
template <typename T>
inline void printVector(std::vector<T>& vec) {
  printVector(vec.size(), vec.data());
}

#ifdef USE_CBC
void printTree(PartialBBDisjunction* const orig_owner,
    OsiSolverInterface* solver, OsiCuts* vpcs, OsiCuts* gmics);
#endif
std::string generateTreePlotString(const VPCEventHandler* eventHandler, const VPCParameters& params,
    const bool saveToFile = false);
std::string generateTikzTreeString(const VPCEventHandler* eventHandler, const VPCParameters& params,
    const int orig_strategy, const double branching_lb, const bool saveToFile);
