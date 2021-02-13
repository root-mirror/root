// @(#)root/thread:$Id$
// Author: Xavier Valls September 2020

/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "ROOT/TExecutor.hxx"

//////////////////////////////////////////////////////////////////////////
///
/// \class ROOT::Internal::TExecutor
/// \brief This class implements the interface to execute the same task
/// multiple times, sequentially or in parallel depending on the execution policy passed
/// as a first parameter on construction, and possibly with different arguments every time.
///
/// ###ROOT::Internal::TExecutor::Map
/// The two possible usages of the Map method are:\n
/// * `Map(F func, unsigned nTimes)`: func is executed nTimes with no arguments
/// * `Map(F func, T& args)`: func is executed on each element of the collection of arguments args
///
/// For either signature, func is executed as many times as needed by a pool of
/// n workers; where n tipically defaults to the number of cores.\n
/// A collection containing the result of each execution is returned.\n
/// **Note:** the user is responsible for the deletion of any object that might
/// be created upon execution of func, returned objects included. ROOT::::Internal::TExecutor never
/// deletes what it returns, it simply forgets it.\n
///
/// \param func
/// \parblock
/// a callable object, such as a lambda expression, an std::function, a
/// functor object or a function that takes zero arguments (for the first signature)
/// or one (for the second signature).
/// \endparblock
/// \param args
/// \parblock
/// a standard vector, a ROOT::TSeq of integer type or an initializer list for the second signature.
/// An integer only for the first.\n
/// \endparblock
///
/// **Note:** in cases where the function to be executed takes more than
/// zero/one argument but all are fixed except zero/one, the function can be wrapped
/// in a lambda or via std::bind to give it the right signature.\n
///
/// #### Return value:
/// An std::vector. The elements in the container
/// will be the objects returned by func.
///
/// ### ROOT::Internal::TExecutor::MapReduce
/// This set of methods behaves exactly like Map, but takes an additional
/// function as a third argument. This function is applied to the set of
/// objects returned by the corresponding Map execution to "squash" them
/// into a single object.
///
/// An integer can be passed as the fourth argument indicating the number of chunks we want to divide our work in.
/// <b>(Note: Please be aware that chunking is only available when the policy is kMultiThread, ignoring this argument in other cases)</b>
/// This may be useful to avoid the overhead introduced when running really short tasks. In this case, the reduction
/// function should be independent of the size of the vector returned by Map due to optimization of the number of
/// chunks.
///
/// #### Examples:
/// ~~~{.cpp}
/// root[] ROOT::Internal::TExecutor pool; auto ten = pool.MapReduce([]() { return 1; }, 10, [](const std::vector<int> &v) { return std::accumulate(v.begin(), v.end(), 0); })
/// root[] ROOT::Internal::TExecutor pool(ROOT::EExecutionPolicy::kMultiProcess); auto hist = pool.MapReduce(CreateAndFillHists, 10, PoolUtils::ReduceObjects);
/// ~~~
///
//////////////////////////////////////////////////////////////////////////


namespace ROOT {
namespace Internal {
TExecutor::TExecutor(ROOT::EExecutionPolicy execPolicy, unsigned nWorkers): fExecPolicy(execPolicy) {
   switch(fExecPolicy) {
      case ROOT::EExecutionPolicy::kSequential:
         fSequentialExecutor = std::make_unique<ROOT::TSequentialExecutor>();
         break;
#ifdef R__USE_IMT
      case ROOT::EExecutionPolicy::kMultiThread:
         fThreadExecutor = std::make_unique<ROOT::TThreadExecutor>(nWorkers);
         break;
#endif
#ifndef R__WIN32
      case ROOT::EExecutionPolicy::kMultiProcess:
         fProcessExecutor = std::make_unique<ROOT::TProcessExecutor>(nWorkers);
         break;
#endif
      default:
         throw std::invalid_argument(
            "Invalid execution policy. Potential issues:\n* kMultiThread policy not available when ROOT is compiled with IMT=OFF.\n* kMultiprocess policy not available on Windows");
   }
}
}
}
