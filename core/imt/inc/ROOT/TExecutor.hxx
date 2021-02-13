// @(#)root/thread:$Id$
// Author: Xavier Valls September 2020

/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/
#ifndef ROOT_TExecutor
#define ROOT_TExecutor

#include "ROOT/RConfig.hxx"
#include "ROOT/RMakeUnique.hxx"
#include "ROOT/TExecutorCRTP.hxx"
#include "ROOT/TSeq.hxx"
#include "ROOT/TSequentialExecutor.hxx"
#ifdef R__USE_IMT
#include "ROOT/TThreadExecutor.hxx"
#endif
#ifndef R__WIN32
#include "ROOT/TProcessExecutor.hxx"
#endif
#include "TROOT.h"
#include "ROOT/EExecutionPolicy.hxx"

#include <initializer_list>
#include <memory>
#include <thread>
#include <type_traits> //std::enable_if, std::result_of
#include <stdexcept> //std::invalid_argument
#include <utility> //std::move

namespace ROOT{

namespace Internal{
class TExecutor: public TExecutorCRTP<TExecutor> {
   friend TExecutorCRTP;
public:

   /// \brief Class constructor. Sets the default execution policy and initializes the corresponding executor.
   /// Defaults to multithreaded execution policy if ROOT is compiled with IMT=ON and IsImplicitMTEnabled. Otherwise it defaults to a serial execution policy
   /// \param nWorkers [optional] Number of parallel workers, only taken into account if the execution policy is kMultiThread
   explicit TExecutor(unsigned nWorkers = 0) :
      TExecutor(ROOT::IsImplicitMTEnabled() ? ROOT::EExecutionPolicy::kMultiThread : ROOT::EExecutionPolicy::kSequential, nWorkers) {}

   /// \brief Class constructor. Sets the execution policy and initializes the corresponding executor.
   /// \param execPolicy Execution policy(kMultiThread, kMultiprocess, kSerial) to process the data
   /// \param nWorkers [optional] Number of parallel workers, only taken into account if the execution policy is kMultiThread
   explicit TExecutor(ROOT::EExecutionPolicy execPolicy, unsigned nWorkers = 0);

   TExecutor(const TExecutor &) = delete;
   TExecutor &operator=(const TExecutor &) = delete;

   /// Return the execution policy the executor is set to
   ROOT::EExecutionPolicy Policy() const { return fExecPolicy; }

   // Map
   //
   using TExecutorCRTP<TExecutor>::Map;

   // MapReduce
   // the late return types also check at compile-time whether redfunc is compatible with func,
   // other than checking that func is compatible with the type of arguments.
   // a static_assert check in TExecutor::Reduce is used to check that redfunc is compatible with the type returned by func
   using TExecutorCRTP<TExecutor>::MapReduce;
   template<class F, class R, class Cond = noReferenceCond<F>>
   auto MapReduce(F func, unsigned nTimes, R redfunc, unsigned nChunks) -> typename std::result_of<F()>::type;
   template<class F, class INTEGER, class R, class Cond = noReferenceCond<F, INTEGER>>
   auto MapReduce(F func, ROOT::TSeq<INTEGER> args, R redfunc, unsigned nChunks) -> typename std::result_of<F(INTEGER)>::type;
   template<class F, class T, class R, class Cond = noReferenceCond<F, T>>
   auto MapReduce(F func, std::initializer_list<T> args, R redfunc, unsigned nChunks) -> typename std::result_of<F(T)>::type;
   template<class F, class T, class R, class Cond = noReferenceCond<F, T>>
   auto MapReduce(F func, std::vector<T> &args, R redfunc, unsigned nChunks) -> typename std::result_of<F(T)>::type;
   template<class F, class T, class R, class Cond = noReferenceCond<F, T>>
   auto MapReduce(F func, const std::vector<T> &args, R redfunc, unsigned nChunks) -> typename std::result_of<F(T)>::type;

   // Reduce
   //
   using TExecutorCRTP<TExecutor>::Reduce;

   unsigned GetPoolSize() const;

private:
   // Implementation of the Map functions declared in the parent class (TExecutorCRTP)
   //
   template<class F, class Cond = noReferenceCond<F>>
   auto MapImpl(F func, unsigned nTimes) -> std::vector<typename std::result_of<F()>::type>;
   template<class F, class INTEGER, class Cond = noReferenceCond<F, INTEGER>>
   auto MapImpl(F func, ROOT::TSeq<INTEGER> args) -> std::vector<typename std::result_of<F(INTEGER)>::type>;
   template<class F, class T, class Cond = noReferenceCond<F, T>>
   auto MapImpl(F func, std::vector<T> &args) -> std::vector<typename std::result_of<F(T)>::type>;
   template<class F, class T, class Cond = noReferenceCond<F, T>>
   auto MapImpl(F func, const std::vector<T> &args) -> std::vector<typename std::result_of<F(T)>::type>;

   ROOT::EExecutionPolicy fExecPolicy;

   // When they are not available, we use a placeholder type instead of TThreadExecutor or TProcessExecutor.
   // The corresponding data members will not be used.
   using Unused_t = ROOT::TSequentialExecutor;
#ifdef R__USE_IMT
# define R__EXECUTOR_THREAD ROOT::TThreadExecutor
#else
# define R__EXECUTOR_THREAD Unused_t
#endif
#ifndef R__WIN32
# define R__EXECUTOR_PROCESS ROOT::TProcessExecutor
#else
# define R__EXECUTOR_PROCESS Unused_t
#endif

   std::unique_ptr<R__EXECUTOR_THREAD> fThreadExecutor;
   std::unique_ptr<R__EXECUTOR_PROCESS> fProcessExecutor;
   std::unique_ptr<ROOT::TSequentialExecutor> fSequentialExecutor;

#undef R__EXECUTOR_THREAD
#undef R__EXECUTOR_PROCESS

   /// \brief Helper class to get the correct return type from the Map function,
   /// necessary to infer the ResolveExecutorAndMap function type
   template<class F, class CONTAINER>
   struct MapRetType {
      using type = typename std::result_of<F(typename CONTAINER::value_type)>::type;
   };

   template<class F>
   struct MapRetType<F, unsigned> {
      using type = typename std::result_of<F()>::type;
   };


   /// \brief Function called from Map to select and execute the correct Executor
   /// according to the set Execution Policy.
   template<class F, class T>
   auto ResolveExecutorAndMap(F func, T&& args) -> std::vector<typename MapRetType<F, typename std::decay<T>::type>::type> {
      std::vector<typename MapRetType<F, typename std::decay<T>::type>::type> res;
      switch(fExecPolicy) {
         case ROOT::EExecutionPolicy::kSequential:
            res = fSequentialExecutor->Map(func, std::forward<T>(args));
            break;
         case ROOT::EExecutionPolicy::kMultiThread:
            res = fThreadExecutor->Map(func, std::forward<T>(args));
            break;
         case ROOT::EExecutionPolicy::kMultiProcess:
            res = fProcessExecutor->Map(func, std::forward<T>(args));
            break;
         default:
            break;
      }
      return res;
   }
};


//////////////////////////////////////////////////////////////////////////
/// \brief Execute a function without arguments several times.
/// Implementation of the Map method.
///
/// \copydetails TExecutorCRTP::Map(F func,unsigned nTimes)
template<class F, class Cond>
auto TExecutor::MapImpl(F func, unsigned nTimes) -> std::vector<typename std::result_of<F()>::type> {
   return ResolveExecutorAndMap(func, nTimes);
}

//////////////////////////////////////////////////////////////////////////
/// \brief Execute a function over a sequence of indexes.
/// Implementation of the Map method.
///
/// \copydetails TExecutorCRTP::Map(F func,ROOT::TSeq<INTEGER> args)
template<class F, class INTEGER, class Cond>
auto TExecutor::MapImpl(F func, ROOT::TSeq<INTEGER> args) -> std::vector<typename std::result_of<F(INTEGER)>::type> {
   return ResolveExecutorAndMap(func, args);
}

//////////////////////////////////////////////////////////////////////////
/// \brief Execute a function over the elements of a vector.
/// Implementation of the Map method.
///
/// \copydetails TExecutorCRTP::Map(F func,std::vector<T> &args)
template<class F, class T, class Cond>
auto TExecutor::MapImpl(F func, std::vector<T> &args) -> std::vector<typename std::result_of<F(T)>::type> {
   return ResolveExecutorAndMap(func, args);
}

//////////////////////////////////////////////////////////////////////////
/// \brief Execute a function over the elements of an immutable vector.
/// Implementation of the Map method.
///
/// \copydetails TExecutorCRTP::Map(F func,const std::vector<T> &args)
template<class F, class T, class Cond>
auto TExecutor::MapImpl(F func, const std::vector<T> &args) -> std::vector<typename std::result_of<F(T)>::type> {
   return ResolveExecutorAndMap(func, args);
}

//////////////////////////////////////////////////////////////////////////
/// \brief Execute a function `nTimes` (Map) and accumulate the results into a single value (Reduce).
/// Benefits from partial reduction into `nChunks` intermediate results if the execution policy is multithreaded.
/// Otherwise, <b>it ignores the nChunks argument</b> and performs a normal MapReduce operation.
///
/// \param func Function to be executed. Must take an element of the sequence passed as second argument as a parameter.
/// \param nTimes Number of times function should be called.
/// \param redfunc Reduction function to combine the results of the calls to `func` into partial results, and these
/// into a final result. Must return the same type as `func`.
/// \param nChunks Number of chunks to split the input data for processing.
/// \return A value result of "reducing" the vector returned by the Map operation into a single object.
template<class F, class R, class Cond>
auto TExecutor::MapReduce(F func, unsigned nTimes, R redfunc, unsigned nChunks) -> typename std::result_of<F()>::type {
   if (fExecPolicy == ROOT::EExecutionPolicy::kMultiThread) {
      return fThreadExecutor->MapReduce(func, nTimes, redfunc, nChunks);
   }
   return Reduce(Map(func, nTimes), redfunc);
}

//////////////////////////////////////////////////////////////////////////
/// \brief Execute a function over a sequence of indexes (Map) and accumulate the results into a single value (Reduce).
/// Benefits from partial reduction into `nChunks` intermediate results if the execution policy is multithreaded.
/// Otherwise, <b>it ignores the nChunks argument</b> and performs a normal MapReduce operation.
///
/// \param func Function to be executed. Must take an element of the sequence passed assecond argument as a parameter.
/// \param args Sequence of indexes to execute `func` on.
/// \param redfunc Reduction function to combine the results of the calls to `func` into partial results, and these
/// into a final result. Must return the same type as `func`.
/// \param nChunks Number of chunks to split the input data for processing.
/// \return A value result of "reducing" the vector returned by the Map operation into a single object.
template<class F, class INTEGER, class R, class Cond>
auto TExecutor::MapReduce(F func, ROOT::TSeq<INTEGER> args, R redfunc, unsigned nChunks) -> typename std::result_of<F(INTEGER)>::type {
   if (fExecPolicy == ROOT::EExecutionPolicy::kMultiThread) {
      return fThreadExecutor->MapReduce(func, args, redfunc, nChunks);
   }
   return Reduce(Map(func, args), redfunc);
}

//////////////////////////////////////////////////////////////////////////
/// \brief Execute a function over the elements of an initializer_list (Map) and accumulate the results into a single value (Reduce).
/// Benefits from partial reduction into `nChunks` intermediate results if the execution policy is multithreaded.
/// Otherwise, <b>it ignores the nChunks argument</b> and performs a normal MapReduce operation.
///
/// \param func Function to be executed. Must take an element of the sequence passed assecond argument as a parameter.
/// \param args initializer_list for a vector to apply `func` on.
/// \param redfunc Reduction function to combine the results of the calls to `func` into partial results, and these
/// into a final result. Must return the same type as `func`.
/// \param nChunks Number of chunks to split the input data for processing.
/// \return A value result of "reducing" the vector returned by the Map operation into a single object.
template<class F, class T, class R, class Cond>
auto TExecutor::MapReduce(F func, std::initializer_list<T> args, R redfunc, unsigned nChunks) -> typename std::result_of<F(T)>::type {
   if (fExecPolicy == ROOT::EExecutionPolicy::kMultiThread) {
      return fThreadExecutor->MapReduce(func, args, redfunc, nChunks);
   }
   return Reduce(Map(func, args), redfunc);
}

//////////////////////////////////////////////////////////////////////////
/// \brief Execute a function over the elements of a vector (Map) and accumulate the results into a single value (Reduce).
/// Benefits from partial reduction into `nChunks` intermediate results if the execution policy is multithreaded.
/// Otherwise, <b>it ignores the nChunks argument</b> and performs a normal MapReduce operation.
///
/// \param func Function to be executed. Must take an element of the sequence passed assecond argument as a parameter.
/// \param args Vector of elements passed as an argument to `func`.
/// \param redfunc Reduction function to combine the results of the calls to `func` into partial results, and these
/// into a final result. Must return the same type as `func`.
/// \param nChunks Number of chunks to split the input data for processing.
/// \return A value result of "reducing" the vector returned by the Map operation into a single object.
template<class F, class T, class R, class Cond>
auto TExecutor::MapReduce(F func, std::vector<T> &args, R redfunc, unsigned nChunks) -> typename std::result_of<F(T)>::type {
   if (fExecPolicy == ROOT::EExecutionPolicy::kMultiThread) {
      return fThreadExecutor->MapReduce(func, args, redfunc, nChunks);
   }
   return Reduce(Map(func, args), redfunc);
}

//////////////////////////////////////////////////////////////////////////
/// \brief Execute a function over the elements of an immutable vector (Map) and accumulate the results into a single value (Reduce).
/// Benefits from partial reduction into `nChunks` intermediate results if the execution policy is multithreaded.
/// Otherwise, <b>it ignores the nChunks argument</b> and performs a normal MapReduce operation.
///
/// \param func Function to be executed. Must take an element of the sequence passed assecond argument as a parameter.
/// \param args Immutable vector, whose elements are passed as an argument to `func`.
/// \param redfunc Reduction function to combine the results of the calls to `func` into partial results, and these
/// into a final result. Must return the same type as `func`.
/// \param nChunks Number of chunks to split the input data for processing.
/// \return A value result of "reducing" the vector returned by the Map operation into a single object.
template<class F, class T, class R, class Cond>
auto TExecutor::MapReduce(F func, const std::vector<T> &args, R redfunc, unsigned nChunks) -> typename std::result_of<F(T)>::type {
   if (fExecPolicy == ROOT::EExecutionPolicy::kMultiThread) {
      return fThreadExecutor->MapReduce(func, args, redfunc, nChunks);
   }
   return Reduce(Map(func, args), redfunc);
}

//////////////////////////////////////////////////////////////////////////
/// \brief Return the number of pooled workers.
///
/// \return The number of workers in the pool in the executor used as a backend.

inline unsigned TExecutor::GetPoolSize() const
{
   unsigned poolSize{0u};
   switch(fExecPolicy){
      case ROOT::EExecutionPolicy::kSequential:
         poolSize = fSequentialExecutor->GetPoolSize();
         break;
      case ROOT::EExecutionPolicy::kMultiThread:
         poolSize = fThreadExecutor->GetPoolSize();
         break;
      case ROOT::EExecutionPolicy::kMultiProcess:
         poolSize = fProcessExecutor->GetPoolSize();
         break;
      default:
         break;
   }
   return poolSize;
}

} // namespace Internal
} // namespace ROOT

#endif
