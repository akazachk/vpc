#include "main_reporting.hpp"

#include <cstdio>

namespace main_reporting {

int wrapUp(
    int retCode,
    int argc,
    char** argv,
    MainAppState& state) {
  const int exitReasonInt = static_cast<int>(state.exitReason);

  std::time_t end_time_t;
  time(&end_time_t);
  struct tm* end_timeinfo = localtime(&end_time_t);
  char end_time_string[25];
  snprintf(end_time_string, sizeof(end_time_string) / sizeof(char), "%s", asctime(end_timeinfo));

  FILE* logfile = state.params.logfile;
  if (logfile != NULL) {
    if (state.params.get(VPCParametersNamespace::PREPROCESS) == 0) {
      printBoundAndGapInfo(state.boundInfo, state.params.logfile);
      printSummaryBBInfo({state.info_nocuts, state.info_mycuts}, state.params.logfile, state.params.get(VPCParametersNamespace::BB_RUNS) == 0);
      printOrigProbInfo(state.origSolver, state.params.logfile);
      printPostCutProbInfo(state.solver, state.cutInfoGMICs, state.cutInfo, state.params.logfile);
      printDisjInfo(state.disjInfo, state.params.logfile);
      printCutInfo(state.cutInfoGMICs, state.cutInfo, state.params.logfile);
      printFullBBInfo({state.info_nocuts, state.info_mycuts}, state.params.logfile, state.params.get(VPCParametersNamespace::BB_RUNS) == 0);
      state.timer.print(state.params.logfile, 2);
    }

#ifdef VPC_VERSION
    fprintf(logfile, "%s,", state.vpc_version_string.substr(0, 8).c_str());
#else
    fprintf(logfile, ",");
#endif
#ifdef VPC_CBC_VERSION
    fprintf(logfile, "%s,", state.cbc_version_string.substr(0, 8).c_str());
#else
    fprintf(logfile, ",");
#endif
#ifdef VPC_CLP_VERSION
    fprintf(logfile, "%s,", state.clp_version_string.substr(0, 8).c_str());
#else
    fprintf(logfile, ",");
#endif
#ifdef USE_GUROBI
    fprintf(logfile, "%s,", state.gurobi_version_string.c_str());
#else
    fprintf(logfile, ",");
#endif
#ifdef USE_CPLEX
    fprintf(logfile, "%s,", state.cplex_version_string.c_str());
#else
    fprintf(logfile, ",");
#endif
    fprintf(logfile, "%s,", state.hostname);
    fprintf(logfile, "%s,", state.cpu_model.c_str());
    fprintf(logfile, "%d,", state.cpu_id);
    fprintf(logfile, "%s,", CglVPC::ExitReasonName[exitReasonInt].c_str());
    fprintf(logfile, "%s,", end_time_string);
    fprintf(logfile, "%.f,", difftime(end_time_t, state.start_time_t));
    fprintf(logfile, "%s,", state.instname.c_str());
    fprintf(logfile, "DONE\n");
    fclose(logfile);
  }

#ifdef TRACE
  printf("\n## Parameter values ##\n");
  printParams(state.params, stdout, 7);
  printf("\n");
  printParams(state.params, stdout, 8);
  printf("\n");

  int NAME_WIDTH = 25;
  printf("\n## Time information ##\n");
  for (size_t i = 0; i < state.overallTimeStatsName.size(); i++) {
    std::string name = state.overallTimeStatsName[i];
    if ((state.params.get(VPCParametersNamespace::intParam::DISJ_TERMS) == 0)
        && (state.params.get(VPCParametersNamespace::stringParam::DISJ_OPTIONS).empty())
        && name.compare(0, 3, "VPC") == 0)
      continue;
    if (state.params.get(VPCParametersNamespace::intParam::BB_RUNS) == 0 && name.compare(0, 2, "BB") == 0)
      continue;
    printf("%-*.*s%s\n", NAME_WIDTH, NAME_WIDTH, (name).c_str(),
        stringValue(state.timer.get_time(name), "%.3f").c_str());
  }
#endif

  printf("%s", state.cut_output.c_str());
  printf("%s", state.bb_output.c_str());

  printf("\n## Exiting VPC generation with reason %s. ##\n", CglVPC::ExitReasonName[exitReasonInt].c_str());
#ifdef VPC_VERSION
  printf("VPC Version: %s\n", state.vpc_version_string.substr(0, 8).c_str());
#endif
#ifdef VPC_CBC_VERSION
  printf("Cbc Version: %s\n", state.cbc_version_string.substr(0, 8).c_str());
#endif
#ifdef VPC_CLP_VERSION
  printf("Clp Version: %s\n", state.clp_version_string.substr(0, 8).c_str());
#endif
#ifdef USE_GUROBI
  printf("Gurobi Version: %s\n", state.gurobi_version_string.c_str());
#endif
#ifdef USE_CPLEX
  printf("CPLEX Version: %s\n", state.cplex_version_string.c_str());
#endif
  fprintf(stdout, "Hostname: %s\n", state.hostname);
  fprintf(stdout, "CPU Model: %s\n", state.cpu_model.c_str());
  fprintf(stdout, "CPU ID: %d\n", state.cpu_id);
  printf("Instance: %s\n", state.instname.c_str());
  if (!state.params.get(VPCParametersNamespace::stringParam::LOGFILE).empty()) {
    printf("Log: %s\n", state.params.get(VPCParametersNamespace::stringParam::LOGFILE).c_str());
  } else {
    printf("Log: stdout\n");
  }
  printf("Start time: %s\n", state.start_time_string);
  printf("End time: %s\n", end_time_string);
  printf("Elapsed time: %.f seconds\n", difftime(end_time_t, state.start_time_t));
  {
    if (argc > 0) {
      printf("Command: ");
    }
    for (int i = 0; i < argc; i++) {
      printf("%s ", argv[i]);
    }
    if (argc > 0) {
      printf("\n");
    }
  }

  if (state.solver) {
    delete state.solver;
  }
  if (state.origSolver) {
    delete state.origSolver;
  }
  if (state.GMICSolver) {
    delete state.GMICSolver;
  }
  if (state.VPCSolver) {
    delete state.VPCSolver;
  }
  if (state.allCutsSolver) {
    delete state.allCutsSolver;
  }

  return retCode;
}

} // namespace main_reporting
