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
    VPCParameters& params,
    CglVPC::ExitReason& exitReason,
    SummaryBoundInfo& boundInfo,
    std::time_t& start_time_t,
    char* start_time_string,
    size_t start_time_string_size,
    char* hostname,
    int hostname_max,
    std::string& cpu_model,
    int& cpu_id,
    std::string& dir,
    std::string& filename_stub,
    std::string& instname,
    std::string& in_file_ext,
    const std::vector<std::string>& overallTimeStatsName) {
  int status = 0;
  exitReason = CglVPC::ExitReason::UNKNOWN;

  printf("## V-Polyhedral Disjunctive Cuts ##\n");
  printf("# Aleksandr M. Kazachkov\n");
  printf("# Based on joint work with Egon Balas\n");
  for (int i = 0; i < argc; i++) {
    std::cout << argv[i] << " ";
  }
  std::cout << std::endl;

  time(&start_time_t);
  struct tm* start_timeinfo = localtime(&start_time_t);
  snprintf(start_time_string, start_time_string_size / sizeof(char), "%s", asctime(start_timeinfo));
  printf("Start time: %s\n", start_time_string);

  {
    const int host_status = gethostname(hostname, hostname_max);
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
        cpu_model = line.substr(pos + 13);
        break;
      }
    }

    cpuInfoFile.close();
#else
    cpu_model = "";
#endif
  }

  {
#ifdef __linux__
    cpu_id = sched_getcpu();
#else
    cpu_id = -1;
#endif
  }

  status = main_input::processArgs(argc, argv, params);
  if (status) {
    return status;
  }

  if (params.get(stringParam::FILENAME).empty()) {
    error_msg(errorstring,
        "No instance file provided. Use -f /path/to/instance.mps to specify the instance.\n");
    exit(1);
  }
  printf("Instance file: %s\n", params.get(stringParam::FILENAME).c_str());

  if (parseFilename(dir, instname, in_file_ext, params.get(stringParam::FILENAME), params.logfile) != 0) {
    error_msg(errorstring,
        "Unable to parse filename: %s. Found: dir=\"%s\", instname=\"%s\",ext=\"%s\".\n",
        params.get(stringParam::FILENAME).c_str(), dir.c_str(),
        instname.c_str(), in_file_ext.c_str());
    exit(1);
  }
  filename_stub = dir + "/" + instname;

  const std::string logname = params.get(stringParam::LOGFILE);
  bool logexists = false;
  if (!logname.empty()) {
    logexists = fexists(logname.c_str());
    params.logfile = fopen(logname.c_str(), "a");
  }

  if (!isInfinity(params.get(doubleParam::IP_OBJ))) {
    boundInfo.ip_obj = params.get(doubleParam::IP_OBJ);
  }
  if (isInfinity(boundInfo.ip_obj) && !params.get(stringParam::OPTFILE).empty()) {
    std::string optfile = params.get(stringParam::OPTFILE);
    std::string csvext = ".csv";
    if (optfile.size() > csvext.size() && optfile.compare(optfile.size() - csvext.size(), csvext.size(), csvext) == 0) {
#ifdef TRACE
      fprintf(stdout, "Reading objective information from \"%s\".\n", optfile.c_str());
#endif
      boundInfo.ip_obj = getObjValueFromFile(optfile, params.get(stringParam::FILENAME), params.logfile);
      params.set(doubleParam::IP_OBJ, boundInfo.ip_obj);
      fprintf(stdout, "Best known IP objective value is %s.\n", stringValue(boundInfo.ip_obj, "%f").c_str());
      if (isInfinity(boundInfo.ip_obj)) {
        warning_msg(warnstring, "Did not find objective value.\n");
      }
    }
  }
  if (params.logfile != NULL) {
    if (!logexists) {
      if (params.get(PREPROCESS) == 0) {
        printHeader(params, overallTimeStatsName);
      } else {
        printPreprocessingHeader(params);
      }
    }
    fprintf(params.logfile, "%s,", instname.c_str());
    if (params.get(PREPROCESS) == 0) {
      printParams(params, params.logfile, 8);
    }
    fflush(params.logfile);
  }

  return status;
}

} // namespace main_startup
