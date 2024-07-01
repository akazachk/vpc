/**
 * @file DisjunctionHelper.cpp
 * @author A. M. Kazachkov
 * @date 2018-02-25
 */
#include "DisjunctionHelper.hpp"

// Project files
#include "SolverHelper.hpp"
#include "SolverInterface.hpp"
#include "TimeStats.hpp"
#include "utility.hpp"
#include "VPCParameters.hpp"
using namespace VPCParametersNamespace;

// Disjunctions
#include "PartialBBDisjunction.hpp"
#include "SplitDisjunction.hpp"

CglVPC::ExitReason setDisjunctions(
    DisjunctionSet* const disjSet,
    const OsiSolverInterface* const si, 
    const VPCParametersNamespace::VPCParameters& params,
    const VPCParametersNamespace::VPCMode& mode) {
  if (mode == VPCParametersNamespace::VPCMode::PARTIAL_BB) {
    if (params.get(intParam::DISJ_TERMS) < 2) {
      return CglVPC::ExitReason::NO_DISJUNCTION_EXIT;
    }
    PartialBBDisjunction* disj = new PartialBBDisjunction(params);
    DisjExitReason status = disj->prepareDisjunction(si);
    disjSet->addDisjunction(disj);
    if (disj) { delete disj; }
    return matchStatus(status);
  } // PARTIAL_BB
  else if (mode == VPCParametersNamespace::VPCMode::SPLITS) {
    if (generateSplitDisjunctions(disjSet, si, params, true)) {
      return CglVPC::ExitReason::SUCCESS_EXIT;
    } else {
      return CglVPC::ExitReason::NO_DISJUNCTION_EXIT;
    }
  } else {
    error_msg(errorstring,
        "Mode that is chosen has not yet been implemented for VPC generation: %s.\n",
        VPCParametersNamespace::VPCModeName[static_cast<int>(mode)].c_str());
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }
  return CglVPC::ExitReason::UNKNOWN;
} /* setDisjunctions */
