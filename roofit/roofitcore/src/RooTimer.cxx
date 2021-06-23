#include "RooTimer.h"
#include "RooTrace.h"


double RooTimer::timing_s() {
  return _timing_s;
}

void RooTimer::set_timing_s(double timing_s) {
  _timing_s = timing_s;
}

void RooTimer::store_timing_in_RooTrace(const std::string &name) {
  RooTrace::objectTiming[name] = _timing_s;  // subscript operator overwrites existing values, insert does not
}

std::vector<RooJsonListFile> RooTimer::timing_outfiles;
std::vector<double> RooTimer::timings;


RooWallTimer::RooWallTimer() {
  start();
}

void RooWallTimer::start() {
  _timing_begin = std::chrono::high_resolution_clock::now();
}

void RooWallTimer::stop() {
  _timing_end = std::chrono::high_resolution_clock::now();
  set_timing_s(std::chrono::duration_cast<std::chrono::nanoseconds>(_timing_end - _timing_begin).count() / 1.e9);
}


RooCPUTimer::RooCPUTimer() {
  start();
}

void RooCPUTimer::start() {
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &_timing_begin);
}

void RooCPUTimer::stop() {
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &_timing_end);
  long long timing_end_nsec = _timing_end.tv_nsec + 1000000000 * _timing_end.tv_sec; // 1'000'000'000 in c++14
  long long timing_begin_nsec = _timing_begin.tv_nsec + 1000000000 * _timing_begin.tv_sec; // 1'000'000'000 in c++14
  set_timing_s((timing_end_nsec - timing_begin_nsec) / 1.e9);
}