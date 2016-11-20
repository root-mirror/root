#include "ROOT/TThreadExecutor.hxx"
#include "tbb/tbb.h"

namespace ROOT{
  TThreadExecutor::TThreadExecutor():fInitTBB(new tbb::task_scheduler_init()){
  }

  TThreadExecutor::TThreadExecutor(size_t nThreads):fInitTBB(new tbb::task_scheduler_init(nThreads)){
  }

  TThreadExecutor::~TThreadExecutor() {
    fInitTBB->terminate();
  }

  void TThreadExecutor::ParallelFor(unsigned int start, unsigned int end, unsigned step, const std::function<void(unsigned int i)> &f){
    tbb::parallel_for(start, end, step, f);
  }

  double TThreadExecutor::ParallelReduce(const std::vector<double> &objs, const std::function<double(double a, double b)> &redfunc){
   return tbb::parallel_reduce(tbb::blocked_range<decltype(objs.begin())>(objs.begin(), objs.end()), double{},
                              [redfunc](tbb::blocked_range<decltype(objs.begin())> const & range, double init) {
                              return std::accumulate(range.begin(), range.end(), init, redfunc);
                              }, redfunc);
  }

  float TThreadExecutor::ParallelReduce(const std::vector<float> &objs, const std::function<float(float a, float b)> &redfunc){
   return tbb::parallel_reduce(tbb::blocked_range<decltype(objs.begin())>(objs.begin(), objs.end()), float{},
                              [redfunc](tbb::blocked_range<decltype(objs.begin())> const & range, float init) {
                              return std::accumulate(range.begin(), range.end(), init, redfunc);
                              }, redfunc);
  }
}