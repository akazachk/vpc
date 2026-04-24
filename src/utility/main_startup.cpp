#include "main_startup.hpp"

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <unistd.h>

#include "main_input.hpp"
#include "utility.hpp"

#ifdef __linux__
#include <sched.h>
#endif

using namespace VPCParametersNamespace;

namespace main_startup {

int startUp(
    int argc,
    char** argv,
    MainAppState& state) {
  int status = 0;
  state.exitReason = CglVPC::ExitReason::UNKNOWN;

  printf("## V-Polyhedral Disjunctive Cuts ##\n");
  printf("# Aleksandr M. Kazachkov\n");
  printf("# Based on joint work with Egon Balas\n");
  for (int i = 0; i < argc; i++) {
    std::cout << argv[i] << " ";
  }
  std::cout << std::endl;

  time(&state.start_time_t);
  struct tm* start_timeinfo = localtime(&state.start_time_t);
  snprintf(state.start_time_string, state.start_time_string_size / sizeof(char), "%s", asctime(start_timeinfo));
  printf("Start time: %s\n", state.start_time_string);

  {
    const int host_status = gethostname(state.hostname, state.hostname_max);
    if (host_status) {
      perror("gethostname");
      return EXIT_FAILURE;
    }
  }

  {
#ifdef __linux__
    std::ifstream cpuInfoFile("/proc/cpuinfo");

    if (!cpuInfoFile.is_open()) {
      std::cerr << "Error opening /proc/cpuinfo" << std::endl;
    }

    std::string line;
    while (std::getline(cpuInfoFile, line)) {
      size_t pos = line.find("model name");
      if (pos != std::string::npos) {
        state.cpu_model = line.substr(pos + 13);
        break;
      }
    }

    cpuInfoFile.close();
#else
    state.cpu_model = "";
#endif
  }

  {
#ifdef __linux__
    state.cpu_id = sched_getcpu();
#else
    state.cpu_id = -1;
#endif
  }

  status = main_input::processArgs(argc, argv, state.params);
  if (status) {
    return status;
  }

  if (state.params.get(stringParam::FILENAME).empty()) {
    error_msg(errorstring,
        "No instance file provided. Use -f /path/to/instance.mps to specify the instance.\n");
    exit(1);
  }
  printf("Instance file: %s\n", state.params.get(stringParam::FILENAME).c_str());

  if (parseFilename(state.dir, state.instname, state.in_file_ext, state.params.get(stringParam::FILENAME), state.params.logfile) != 0) {
    error_msg(errorstring,
        "Unable to parse filename: %s. Found: dir=\"%s\", instname=\"%s\",ext=\"%s\".\n",
        state.params.get(stringParam::FILENAME).c_str(), state.dir.c_str(),
        state.instname.c_str(), state.in_file_ext.c_str());
    exit(1);
  }
  state.filename_stub = state.dir + "/" + state.instname;

  const std::string logname = state.params.get(stringParam::LOGFILE);
  bool logexists = false;
  if (!logname.empty()) {
    logexists = fexists(logname.c_str());
    state.params.logfile = fopen(logname.c_str(), "a");
  }

  if (!isInfinity(state.params.get(doubleParam::IP_OBJ))) {
    state.boundInfo.ip_obj = state.params.get(doubleParam::IP_OBJ);
  }
  if (isInfinity(state.boundInfo.ip_obj) && !state.params.get(stringParam::OPTFILE).empty()) {
    std::string optfile = state.params.get(stringParam::OPTFILE);
    std::string csvext = ".csv";
    if (optfile.size() > csvext.size() && optfile.compare(optfile.size() - csvext.size(), csvext.size(), csvext) == 0) {
#ifdef TRACE
      fprintf(stdout, "Reading objective information from \"%s\".\n", optfile.c_str());
#endif
      state.boundInfo.ip_obj = getObjValueFromFile(optfile, state.params.get(stringParam::FILENAME), state.params.logfile);
      state.params.set(doubleParam::IP_OBJ, state.boundInfo.ip_obj);
      fprintf(stdout, "Best known IP objective value is %s.\n", stringValue(state.boundInfo.ip_obj, "%f").c_str());
      if (isInfinity(state.boundInfo.ip_obj)) {
        warning_msg(warnstring, "Did not find objective value.\n");
      }
    }
  }
  if (state.params.logfile != NULL) {
    if (!logexists) {
      if (state.params.get(PREPROCESS) == 0) {
        printHeader(state.params, state.overallTimeStatsName);
      } else {
        printPreprocessingHeader(state.params);
      }
    }
    fprintf(state.params.logfile, "%s,", state.instname.c_str());
    if (state.params.get(PREPROCESS) == 0) {
      printParams(state.params, state.params.logfile, 8);
    }
    fflush(state.params.logfile);
  }

  return status;
}

} // namespace main_startup
