#include "main_reporting.hpp"

#include <cstdio>

namespace main_reporting {

int wrapUp(
    int retCode,
    int argc,
    char** argv,
    VPCParametersNamespace::VPCParameters& params,
    TimeStats& timer,
    const std::vector<std::string>& overallTimeStatsName,
    CglVPC::ExitReason exitReason,
    std::time_t start_time_t,
    const char* start_time_string,
    const char* hostname,
    const std::string& cpu_model,
    int cpu_id,
    const std::string& instname,
    const std::string& cut_output,
    const std::string& bb_output,
    const SummaryBoundInfo& boundInfo,
    const SummaryBBInfo& info_nocuts,
    const SummaryBBInfo& info_mycuts,
    const SummaryCutInfo& cutInfoGMICs,
    const SummaryCutInfo& cutInfo,
    const SummaryDisjunctionInfo& disjInfo,
    OsiSolverInterface*& solver,
    OsiSolverInterface*& origSolver,
    OsiSolverInterface*& GMICSolver,
    OsiSolverInterface*& VPCSolver,
    OsiSolverInterface*& allCutsSolver,
    const std::string& vpc_version_string,
    const std::string& cbc_version_string,
    const std::string& clp_version_string,
    const std::string& gurobi_version_string,
    const std::string& cplex_version_string) {
  const int exitReasonInt = static_cast<int>(exitReason);

  std::time_t end_time_t;
  time(&end_time_t);
  struct tm* end_timeinfo = localtime(&end_time_t);
  char end_time_string[25];
  snprintf(end_time_string, sizeof(end_time_string) / sizeof(char), "%s", asctime(end_timeinfo));

  FILE* logfile = params.logfile;
  if (logfile != NULL) {
    if (params.get(VPCParametersNamespace::PREPROCESS) == 0) {
      printBoundAndGapInfo(boundInfo, params.logfile);
      printSummaryBBInfo({info_nocuts, info_mycuts}, params.logfile, params.get(VPCParametersNamespace::BB_RUNS) == 0);
      printOrigProbInfo(origSolver, params.logfile);
      printPostCutProbInfo(solver, cutInfoGMICs, cutInfo, params.logfile);
      printDisjInfo(disjInfo, params.logfile);
      printCutInfo(cutInfoGMICs, cutInfo, params.logfile);
      printFullBBInfo({info_nocuts, info_mycuts}, params.logfile, params.get(VPCParametersNamespace::BB_RUNS) == 0);
      timer.print(params.logfile, 2);
    }

#ifdef VPC_VERSION
    fprintf(logfile, "%s,", vpc_version_string.substr(0, 8).c_str());
#else
    fprintf(logfile, ",");
#endif
#ifdef VPC_CBC_VERSION
    fprintf(logfile, "%s,", cbc_version_string.substr(0, 8).c_str());
#else
    fprintf(logfile, ",");
#endif
#ifdef VPC_CLP_VERSION
    fprintf(logfile, "%s,", clp_version_string.substr(0, 8).c_str());
#else
    fprintf(logfile, ",");
#endif
#ifdef USE_GUROBI
    fprintf(logfile, "%s,", gurobi_version_string.c_str());
#else
    fprintf(logfile, ",");
#endif
#ifdef USE_CPLEX
    fprintf(logfile, "%s,", cplex_version_string.c_str());
#else
    fprintf(logfile, ",");
#endif
    fprintf(logfile, "%s,", hostname);
    fprintf(logfile, "%s,", cpu_model.c_str());
    fprintf(logfile, "%d,", cpu_id);
    fprintf(logfile, "%s,", CglVPC::ExitReasonName[exitReasonInt].c_str());
    fprintf(logfile, "%s,", end_time_string);
    fprintf(logfile, "%.f,", difftime(end_time_t, start_time_t));
    fprintf(logfile, "%s,", instname.c_str());
    fprintf(logfile, "DONE\n");
    fclose(logfile);
  }

#ifdef TRACE
  printf("\n## Parameter values ##\n");
  printParams(params, stdout, 7);
  printf("\n");
  printParams(params, stdout, 8);
  printf("\n");

  int NAME_WIDTH = 25;
  printf("\n## Time information ##\n");
  for (size_t i = 0; i < overallTimeStatsName.size(); i++) {
    std::string name = overallTimeStatsName[i];
    if ((params.get(VPCParametersNamespace::intParam::DISJ_TERMS) == 0)
        && (params.get(VPCParametersNamespace::stringParam::DISJ_OPTIONS).empty())
        && name.compare(0, 3, "VPC") == 0)
      continue;
    if (params.get(VPCParametersNamespace::intParam::BB_RUNS) == 0 && name.compare(0, 2, "BB") == 0)
      continue;
    printf("%-*.*s%s\n", NAME_WIDTH, NAME_WIDTH, (name).c_str(),
        stringValue(timer.get_time(name), "%.3f").c_str());
  }
#endif

  printf("%s", cut_output.c_str());
  printf("%s", bb_output.c_str());

  printf("\n## Exiting VPC generation with reason %s. ##\n", CglVPC::ExitReasonName[exitReasonInt].c_str());
#ifdef VPC_VERSION
  printf("VPC Version: %s\n", vpc_version_string.substr(0, 8).c_str());
#endif
#ifdef VPC_CBC_VERSION
  printf("Cbc Version: %s\n", cbc_version_string.substr(0, 8).c_str());
#endif
#ifdef VPC_CLP_VERSION
  printf("Clp Version: %s\n", clp_version_string.substr(0, 8).c_str());
#endif
#ifdef USE_GUROBI
  printf("Gurobi Version: %s\n", gurobi_version_string.c_str());
#endif
#ifdef USE_CPLEX
  printf("CPLEX Version: %s\n", cplex_version_string.c_str());
#endif
  fprintf(stdout, "Hostname: %s\n", hostname);
  fprintf(stdout, "CPU Model: %s\n", cpu_model.c_str());
  fprintf(stdout, "CPU ID: %d\n", cpu_id);
  printf("Instance: %s\n", instname.c_str());
  if (!params.get(VPCParametersNamespace::stringParam::LOGFILE).empty()) {
    printf("Log: %s\n", params.get(VPCParametersNamespace::stringParam::LOGFILE).c_str());
  } else {
    printf("Log: stdout\n");
  }
  printf("Start time: %s\n", start_time_string);
  printf("End time: %s\n", end_time_string);
  printf("Elapsed time: %.f seconds\n", difftime(end_time_t, start_time_t));
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

  if (solver) {
    delete solver;
  }
  if (origSolver) {
    delete origSolver;
  }
  if (GMICSolver) {
    delete GMICSolver;
  }
  if (VPCSolver) {
    delete VPCSolver;
  }
  if (allCutsSolver) {
    delete allCutsSolver;
  }

  return retCode;
}

} // namespace main_reporting
