#include "ROOT/TThreadExecutor.hxx"
#include "ROOT/TTaskGroup.hxx"

#include "TROOT.h"

#if !defined(_MSC_VER)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#endif

#include "tbb/tbb.h"

#if !defined(_MSC_VER)
#pragma GCC diagnostic pop
#endif

//////////////////////////////////////////////////////////////////////////
///
/// \class ROOT::TThreadExecutor
/// \ingroup Parallelism
/// \brief This class provides a simple interface to execute the same task
/// multiple times in parallel, possibly with different arguments every
/// time. This mimics the behaviour of python's pool.Map method.
///
/// ### ROOT::TThreadExecutor::Map
/// This class inherits its interfaces from ROOT::TExecutor\n.
/// The two possible usages of the Map method are:\n
/// * Map(F func, unsigned nTimes): func is executed nTimes with no arguments
/// * Map(F func, T& args): func is executed on each element of the collection of arguments args
///
/// For either signature, func is executed as many times as needed by a pool of
/// nThreads threads; It defaults to the number of cores.\n
/// A collection containing the result of each execution is returned.\n
/// **Note:** the user is responsible for the deletion of any object that might
/// be created upon execution of func, returned objects included: ROOT::TThreadExecutor never
/// deletes what it returns, it simply forgets it.\n
///
/// \param func
/// \parblock
/// a lambda expression, an std::function, a loaded macro, a
/// functor class or a function that takes zero arguments (for the first signature)
/// or one (for the second signature).
/// \endparblock
/// \param args
/// \parblock
/// a standard vector, a ROOT::TSeq of integer type or an initializer list for the second signature.
/// An integer only for the first.
/// \endparblock
/// **Note:** in cases where the function to be executed takes more than
/// zero/one argument but all are fixed except zero/one, the function can be wrapped
/// in a lambda or via std::bind to give it the right signature.\n
///
/// #### Return value:
/// An std::vector. The elements in the container
/// will be the objects returned by func.
///
///
/// #### Examples:
///
/// ~~~{.cpp}
/// root[] ROOT::TThreadExecutor pool; auto hists = pool.Map(CreateHisto, 10);
/// root[] ROOT::TThreadExecutor pool(2); auto squares = pool.Map([](int a) { return a*a; }, {1,2,3});
/// ~~~
///
/// ### ROOT::TThreadExecutor::MapReduce
/// This set of methods behaves exactly like Map, but takes an additional
/// function as a third argument. This function is applied to the set of
/// objects returned by the corresponding Map execution to "squash" them
/// to a single object. This function should be independent of the size of
/// the vector returned by Map due to optimization of the number of chunks.
///
/// If this function is a binary operator, the "squashing" will be performed in parallel.
/// This is exclusive to ROOT::TThreadExecutor and not any other ROOT::TExecutor-derived classes.\n
/// An integer can be passed as the fourth argument indicating the number of chunks we want to divide our work in.
/// This may be useful to avoid the overhead introduced when running really short tasks.
///
/// #### Examples:
/// ~~~{.cpp}
/// root[] ROOT::TThreadExecutor pool; auto ten = pool.MapReduce([]() { return 1; }, 10, [](std::vector<int> v) { return std::accumulate(v.begin(), v.end(), 0); })
/// root[] ROOT::TThreadExecutor pool; auto hist = pool.MapReduce(CreateAndFillHists, 10, PoolUtils::ReduceObjects);
/// ~~~
///
//////////////////////////////////////////////////////////////////////////

/*
VERY IMPORTANT NOTE ABOUT WORK ISOLATION

We enclose the parallel_for and parallel_reduce invocations in a
task_arena::isolate because we want to prevent a thread to start executing an
outer task when the task it's running spawned subtasks, e.g. with a parallel_for,
and is waiting on inner tasks to be completed.

While this change has a negligible performance impact, it has benefits for
several applications, for example big parallelised HEP frameworks and
RDataFrame analyses.
- For HEP Frameworks, without work isolation, it can happen that a huge
framework task is pulled by a yielding ROOT task.
This causes to delay the processing of the event which is interrupted by the
long task.
For example, work isolation avoids that during the wait due to the parallel
flushing of baskets, a very long simulation task is pulled in by the idle task.
- For RDataFrame analyses we want to guarantee that each entry is processed from
the beginning to the end without TBB interrupting it to pull in other work items.
As a corollary, the usage of ROOT (or TBB in work isolation mode) in actions
and transformations guarantee that each entry is processed from the beginning to
the end without being interrupted by the processing of outer tasks.
*/

namespace ROOT {
namespace Internal {

/// A helper function to implement the TThreadExecutor::ParallelReduce methods
template<typename T>
static T ParallelReduceHelper(const std::vector<T> &objs, const std::function<T(T a, T b)> &redfunc)
{
   using BRange_t = tbb::blocked_range<decltype(objs.begin())>;

   auto pred = [redfunc](BRange_t const & range, T init) {
      return std::accumulate(range.begin(), range.end(), init, redfunc);
   };

   BRange_t objRange(objs.begin(), objs.end());

   return tbb::this_task_arena::isolate([&]{
      return tbb::parallel_reduce(objRange, T{}, pred, redfunc);
   });

}

} // End NS Internal
} // End NS ROOT

namespace ROOT {

   //////////////////////////////////////////////////////////////////////////
   /// Class constructor.
   /// If the scheduler is active (e.g. because another TThreadExecutor is in flight, or ROOT::EnableImplicitMT() was
   /// called), work with the current pool of threads.
   /// If not, initialize the pool of threads, spawning nThreads. nThreads' default value, 0, initializes the
   /// pool with as many logical threads as are available in the system (see NLogicalCores in TPoolManager.cxx).
   ///
   /// At construction time, TThreadExecutor automatically enables ROOT's thread-safety locks as per calling
   /// ROOT::EnableThreadSafety().
   TThreadExecutor::TThreadExecutor(UInt_t nThreads)
   {
      ROOT::EnableThreadSafety();

      auto current = ROOT::Internal::TPoolManager::GetPoolSize();
      if (nThreads && current && (current != nThreads))
      {
         Warning("TThreadExecutor", "There's already an active pool of threads. Proceeding with the current %d threads", current);
      }
      fSched = ROOT::Internal::GetPoolManager(nThreads);
   }

   void TThreadExecutor::ParallelFor(unsigned int start, unsigned int end, unsigned step, const std::function<void(unsigned int i)> &f)
   {
      tbb::this_task_arena::isolate([&]{
         tbb::parallel_for(start, end, step, f);
      });
   }

   double TThreadExecutor::ParallelReduce(const std::vector<double> &objs, const std::function<double(double a, double b)> &redfunc)
   {
      return ROOT::Internal::ParallelReduceHelper<double>(objs, redfunc);
   }

   float TThreadExecutor::ParallelReduce(const std::vector<float> &objs, const std::function<float(float a, float b)> &redfunc)
   {
      return ROOT::Internal::ParallelReduceHelper<float>(objs, redfunc);
   }

   unsigned TThreadExecutor::GetPoolSize(){
      return ROOT::Internal::TPoolManager::GetPoolSize();
   }

}
