//============================================================================
// Name        : timestats.hpp
// Author      : A. M. Kazachkov and S. Nadadarajah
// Version     : 2018-Dec-25
// Description : For various statistics, such as timer
//
// This code is licensed under the terms of the Eclipse Public License (EPL).
//============================================================================

#pragma once

#include <assert.h>
#include <cstdlib>
#include <string>
#include <ctime>
#include <map>
#include <vector>

#define UNINITIALIZED_STAT    -1    // initial value stat field receives
#define INITIAL_STAT_SIZE  100    // initial number of statistics considered

/**
 * Struct for C-string comparison
 */
struct ltstr {
  bool operator()(const std::string &s1, const std::string &s2) const {
    return (strcmp(s1.c_str(), s2.c_str()) < 0);
  }
};

/**
 * Struct for statistic value
 */
struct data_t {
  int id;
  data_t() :
      id(UNINITIALIZED_STAT) {
  }
};

/**
 * Class to collect general statistics on the code
 */
class TimeStats {
public:
  /** Constructor: Just reserve memory */
  TimeStats() {
    timer_start.reserve(INITIAL_STAT_SIZE);
    value.reserve(INITIAL_STAT_SIZE);
    timer_running.reserve(INITIAL_STAT_SIZE);
  }

  /** Register a new statistic. The initial value is 0 by default */
  int register_name(const std::string &name, clock_t initial_value = 0);

  /** Start timer for a time statistic */
  void start_timer(const std::string &name);
  /** Start timer (by id) for a time statistic */
  void start_timer(int id);

  /** End timer for a time statistic, accumulating the result (in seconds) */
  void end_timer(const std::string &name);
  /** End timer (by id) for a time statistic, accumulating the result (in seconds) */
  void end_timer(int id);
  /** End all timers */
  void end_all();

  /** Add value for a numerical statistic by name. Default value to add is 1. */
  void add_value(const std::string &name, clock_t val = 1);
  /** Add value for a numerical statistic by id. Default value to add is 1. */
  void add_value(int id, clock_t val = 1);

  /** Add time for a numerical statistic by name. Default value to add is 1. */
  void add_time(const std::string &name, clock_t val = 1);
  /** Add time for a numerical statistic by id. Default value to add is 1. */
  void add_time(int id, clock_t val = 1);

  /** Get value for a statistic */
  clock_t get_value(const std::string &name) const;
  /** Get value for a statistic */
  clock_t get_value(const int id) const;

  /** Return value (by name) interpreted as time */
  double get_time(const std::string &name) const;
  /** Return value (by id) interpreted as time */
  double get_time(int id) const;

  /** Return value (by name) interpreted as time, taking current time clock for measure */
  double get_current_time(const std::string &name) const;
  /** Return value (by id) interpreted as time, taking current time clock for measure */
  double get_current_time(const int id) const;

  /** Get value for a statistic interpreted as time, using total time as measure */
  double get_total_time(const std::string &name) const;
  /** Get value for a statistic interpreted as time, using total time as measure */
  double get_total_time(const int id) const;

  /** Return id of name */
  int get_id(const std::string &name) const;

private:
  std::map<std::string, data_t, ltstr> name_to_id; /** map from name to stat identifier */ // string key should probably be const
  std::vector<clock_t> timer_start; /** timer start for statistic */
  std::vector<clock_t> value; /** statistic value */
  std::vector<bool> timer_running; /** is the timer running? */
}; /* TimeStats */

/**
 * -------------------------------------------------------------
 * Inline Implementations
 * -------------------------------------------------------------
 */

/**
 * Register a new statistic
 * */
inline int TimeStats::register_name(const std::string &name,
    clock_t initial_value) {
  name_to_id[name].id = value.size();
  value.push_back(initial_value);
  timer_start.push_back(clock());
  timer_running.push_back(false);
  return value.size() - 1;
}

/**
 * Start the timer for a statistic
 */
inline void TimeStats::start_timer(const std::string &name) {
  start_timer(get_id(name));
}

/**
 * Start the timer (by id) for a statistic
 */
inline void TimeStats::start_timer(int id) {
  //assert(id >= 0 && id < (int )value.size());
  if (id >= 0 && id < (int )value.size()) {
    timer_start[id] = clock();
    timer_running[id] = true;
  }
}

/**
 * End the timer for a statistic, accumulating the result (in seconds)
 */
inline void TimeStats::end_timer(const std::string &name) {
  end_timer(get_id(name));
}

/**
 * End the timer for a statistic (by id), accumulating the result (in seconds)
 */
inline void TimeStats::end_timer(int id) {
  //assert(id >= 0 && id < (int )value.size());
  if (id >= 0 && id < (int )value.size()) {
    value[id] += clock() - timer_start[id];
    timer_running[id] = false;
  }
}

/**
 * End all timers
 */
inline void TimeStats::end_all() {
  for (int id = 0; id < (int) value.size(); id++)
    end_timer(id);
} /* end_all */

/**
 * Add value for a numerical statistic. Default value to add is 1.
 */
inline void TimeStats::add_value(const std::string &name, clock_t val) {
  add_value(get_id(name), val);
}

/**
 * Add value by id
 */
inline void TimeStats::add_value(int id, clock_t val) {
  //assert(id >= 0 && id < (int )value.size());
  if (id >= 0 && id < (int )value.size())
    value[id] += val;
}

/**
 * Add time for a numerical statistic. Default value to add is 1.
 */
inline void TimeStats::add_time(const std::string &name, clock_t val) {
  add_time(get_id(name), val);
}

/**
 * Add time by id
 */
inline void TimeStats::add_time(int id, clock_t val) {
  //assert(id >= 0 && id < (int )value.size());
  if (id >= 0 && id < (int )value.size())
    value[id] += val * CLOCKS_PER_SEC;
}

/**
 * Get value for a statistic
 */
inline clock_t TimeStats::get_value(const std::string &name) const {
  return (get_value(get_id(name)));
}

/**
 * Get value for a statistic
 */
inline clock_t TimeStats::get_value(const int id) const {
  //assert(id >= 0 && id < (int )value.size());
  if (id >= 0 && id < (int )value.size())
    return (value[id]);
  else
    return -1;
}

/**
 * Get value for a statistic interpreted as time
 */
inline double TimeStats::get_time(const std::string &name) const {
  return get_time(get_id(name));
}

/**
 * Get value for a statistic interpreted as time
 */
inline double TimeStats::get_time(int id) const {
  //assert(id >= 0 && id < (int )value.size());
  if (id >= 0 && id < (int )value.size())
    return ((double) (value[id])) / (double) CLOCKS_PER_SEC;
  else
    return -1.0;
} /* get_time (by id) */

/**
 * Get value for a statistic interpreted as time, using current time as measure
 */
inline double TimeStats::get_current_time(const std::string &name) const {
  return get_current_time(get_id(name));
}

/**
 * Get value for a statistic interpreted as time, using current time as measure
 */
inline double TimeStats::get_current_time(const int id) const {
  //assert(id >= 0 && id < (int )value.size());
  if (id >= 0 && id < (int )value.size())
    return ((double) (clock() - timer_start[id])) / (double) CLOCKS_PER_SEC;
  else
    return -1.0;
} /* get_current_time (by id) */

/**
 * Get value for a statistic interpreted as time, using total time as measure
 */
inline double TimeStats::get_total_time(const std::string &name) const {
  return get_total_time(get_id(name));
}

/**
 * Get value for a statistic interpreted as time, using total time as measure
 */
inline double TimeStats::get_total_time(const int id) const {
  //assert(id >= 0 && id < (int )value.size());
  if (id >= 0 && id < (int )value.size()) {
    if (timer_running[id]) {
      return ((double) (value[id] + clock() - timer_start[id])) / (double) CLOCKS_PER_SEC;
    } else {
      return get_time(id);
    }
  } else {
    return -1.0;
  }
} /* get_total_time (by id) */

/**
 * Check and return id
 */
inline int TimeStats::get_id(const std::string &name) const {
  std::map<std::string, data_t, ltstr>::const_iterator it = name_to_id.find(name);
  if (it != name_to_id.end()) {
    if (it->second.id != UNINITIALIZED_STAT)
      return (it->second.id);
    else
      return -1;
  } else {
    return -1;
  }
} /* get_id */
