#ifndef ROO_TIMER
#define ROO_TIMER

#include <chrono>
#include <ctime>
#include <string>
#include <vector>
#include "RooJsonListFile.h"

class RooTimer {
public:
  virtual void start() = 0;
  virtual void stop() = 0;
  double timing_s();
  void set_timing_s(double timing_s);
  void store_timing_in_RooTrace(const std::string &name);

  static std::vector<RooJsonListFile> timing_outfiles;
  static std::vector<double> timings;

private:
  double _timing_s;
};

class RooWallTimer: public RooTimer {
public:
  RooWallTimer();
  virtual void start();
  virtual void stop();

private:
  std::chrono::time_point<std::chrono::high_resolution_clock> _timing_begin, _timing_end;
};

/// @class RooCPUTimer
/// Measures the CPU time on the local process. Note that for multi-process runs,
/// e.g. when using RooRealMPFE, the child process CPU times are not included!
/// Use a separate timer in child processes to measure their CPU timing.
class RooCPUTimer: public RooTimer {
public:
  RooCPUTimer();
  virtual void start();
  virtual void stop();

private:
  struct timespec _timing_begin, _timing_end;
};

#endif
