/* @(#)root/multiproc:$Id$ */
// Author: Enrico Guiraud July 2015
// Modified: G Ganis Jan 2017

/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TProcessExecutor
#define ROOT_TProcessExecutor

#include "MPCode.h"
#include "MPSendRecv.h"
#include "PoolUtils.h"
#include "ROOT/RConfig.hxx" //R__DEPRECATED
#include "ROOT/TExecutorCRTP.hxx"
#include "ROOT/TSeq.hxx"
#include "TError.h"
#include "TFileCollection.h"
#include "TFileInfo.h"
#include "THashList.h"
#include "TMPClient.h"
#include "TMPWorkerExecutor.h"

#include <algorithm> //std::generate
#include <numeric> //std::iota
#include <string>
#include <type_traits> //std::result_of, std::enable_if
#include <functional> //std::reference_wrapper
#include <vector>

namespace ROOT {

class TProcessExecutor : public TExecutorCRTP<TProcessExecutor>, private TMPClient {
   friend TExecutorCRTP;
public:
   explicit TProcessExecutor(unsigned nWorkers = 0); //default number of workers is the number of processors
   ~TProcessExecutor() = default;
   //it doesn't make sense for a TProcessExecutor to be copied
   TProcessExecutor(const TProcessExecutor &) = delete;
   TProcessExecutor &operator=(const TProcessExecutor &) = delete;

   // Map
   //
   using TExecutorCRTP<TProcessExecutor>::Map;

   // MapReduce
   // Redefinition of the MapReduce classes of the base class, to adapt them to
   // TProcessExecutor's logic
   using TExecutorCRTP<TProcessExecutor>::MapReduce;
   template<class F, class R, class Cond = noReferenceCond<F>>
   auto MapReduce(F func, unsigned nTimes, R redfunc) -> typename std::result_of<F()>::type;
   template<class F, class T, class R, class Cond = noReferenceCond<F, T>>
   auto MapReduce(F func, std::vector<T> &args, R redfunc) -> typename std::result_of<F(T)>::type;
   template<class F, class T, class R, class Cond = noReferenceCond<F, T>>
   auto MapReduce(F func, const std::vector<T> &args, R redfunc) -> typename std::result_of<F(T)>::type;

   // Reduce
   //
   using TExecutorCRTP<TProcessExecutor>::Reduce;

   void SetNWorkers(unsigned n) { TMPClient::SetNWorkers(n); }

   /// \deprecated Use GetPoolSize()
   unsigned GetNWorkers() const _R__DEPRECATED_626("Please use GetPoolSize instead") { return TMPClient::GetNWorkers(); }

   //////////////////////////////////////////////////////////////////////////
   /// \brief Return the number of pooled parallel workers.
   ///
   /// \return The number of workers in the pool.
   unsigned GetPoolSize() const { return TMPClient::GetNWorkers(); }

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

   template<class T> void Collect(std::vector<T> &reslist);
   template<class T> void HandlePoolCode(MPCodeBufPair &msg, TSocket *sender, std::vector<T> &reslist);

   void Reset();
   void ReplyToFuncResult(TSocket *s);
   void ReplyToIdle(TSocket *s);

   unsigned fNProcessed; ///< number of arguments already passed to the workers
   unsigned fNToProcess; ///< total number of arguments to pass to the workers

   /// A collection of the types of tasks that TProcessExecutor can execute.
   /// It is used to interpret in the right way and properly reply to the
   /// messages received (see, for example, TProcessExecutor::HandleInput)
   enum class ETask : unsigned char {
      kNoTask,       ///< no task is being executed
      kMap,          ///< a Map method with no arguments is being executed
      kMapWithArg,   ///< a Map method with arguments is being executed
      kMapRed,       ///< a MapReduce method with no arguments is being executed
      kMapRedWithArg ///< a MapReduce method with arguments is being executed
   };

   ETask fTaskType = ETask::kNoTask; ///< the kind of task that is being executed, if any
};


/************ TEMPLATE METHODS IMPLEMENTATION ******************/

//////////////////////////////////////////////////////////////////////////
/// \brief Execute a function without arguments several times in parallel.
/// Implementation of the Map method.
///
/// \copydetails TExecutorCRTP::Map(F func,unsigned nTimes)
template<class F, class Cond>
auto TProcessExecutor::MapImpl(F func, unsigned nTimes) -> std::vector<typename std::result_of<F()>::type>
{
   using retType = decltype(func());
   //prepare environment
   Reset();
   fTaskType = ETask::kMap;

   //fork max(nTimes, fNWorkers) times
   unsigned oldNWorkers = GetPoolSize();
   if (nTimes < oldNWorkers)
      SetNWorkers(nTimes);
   TMPWorkerExecutor<F> worker(func);
   bool ok = Fork(worker);
   SetNWorkers(oldNWorkers);
   if (!ok)
   {
      Error("TProcessExecutor::Map", "[E][C] Could not fork. Aborting operation.");
      return std::vector<retType>();
   }

   //give out tasks
   fNToProcess = nTimes;
   std::vector<retType> reslist;
   reslist.reserve(fNToProcess);
   fNProcessed = Broadcast(MPCode::kExecFunc, fNToProcess);

   //collect results, give out other tasks if needed
   Collect(reslist);

   //clean-up and return
   ReapWorkers();
   fTaskType = ETask::kNoTask;
   return reslist;
}

//////////////////////////////////////////////////////////////////////////
/// \brief Execute a function over the elements of a vector in parallel
/// Implementation of the Map method.
///
/// \copydetails TExecutorCRTP::Map(F func,std::vector<T> &args)
template<class F, class T, class Cond>
auto TProcessExecutor::MapImpl(F func, std::vector<T> &args) -> std::vector<typename std::result_of<F(T)>::type>
{
   //check whether func is callable
   using retType = decltype(func(args.front()));
   //prepare environment
   Reset();
   fTaskType = ETask::kMapWithArg;

   //fork max(args.size(), fNWorkers) times
   //N.B. from this point onwards, args is filled with undefined (but valid) values, since TMPWorkerExecutor moved its content away
   unsigned oldNWorkers = GetPoolSize();
   if (args.size() < oldNWorkers)
      SetNWorkers(args.size());
   TMPWorkerExecutor<F, T> worker(func, args);
   bool ok = Fork(worker);
   SetNWorkers(oldNWorkers);
   if (!ok)
   {
      Error("TProcessExecutor::Map", "[E][C] Could not fork. Aborting operation.");
      return std::vector<retType>();
   }

   //give out tasks
   fNToProcess = args.size();
   std::vector<retType> reslist;
   reslist.reserve(fNToProcess);
   std::vector<unsigned> range(fNToProcess);
   std::iota(range.begin(), range.end(), 0);
   fNProcessed = Broadcast(MPCode::kExecFuncWithArg, range);

   //collect results, give out other tasks if needed
   Collect(reslist);

   //clean-up and return
   ReapWorkers();
   fTaskType = ETask::kNoTask;
   return reslist;
}

//////////////////////////////////////////////////////////////////////////
/// \brief Execute a function over the elements of an immutable vector in parallel
/// Implementation of the Map method.
///
/// \copydetails TExecutorCRTP::Map(F func,const std::vector<T> &args)
template<class F, class T, class Cond>
auto TProcessExecutor::MapImpl(F func, const std::vector<T> &args) -> std::vector<typename std::result_of<F(T)>::type>
{
   //check whether func is callable
   using retType = decltype(func(args.front()));
   //prepare environment
   Reset();
   fTaskType = ETask::kMapWithArg;

   //fork max(args.size(), fNWorkers) times
   //N.B. from this point onwards, args is filled with undefined (but valid) values, since TMPWorkerExecutor moved its content away
   unsigned oldNWorkers = GetPoolSize();
   if (args.size() < oldNWorkers)
      SetNWorkers(args.size());
   TMPWorkerExecutor<F, T> worker(func, args);
   bool ok = Fork(worker);
   SetNWorkers(oldNWorkers);
   if (!ok)
   {
      Error("TProcessExecutor::Map", "[E][C] Could not fork. Aborting operation.");
      return std::vector<retType>();
   }

   //give out tasks
   fNToProcess = args.size();
   std::vector<retType> reslist;
   reslist.reserve(fNToProcess);
   std::vector<unsigned> range(fNToProcess);
   std::iota(range.begin(), range.end(), 0);
   fNProcessed = Broadcast(MPCode::kExecFuncWithArg, range);

   //collect results, give out other tasks if needed
   Collect(reslist);

   //clean-up and return
   ReapWorkers();
   fTaskType = ETask::kNoTask;
   return reslist;
}

//////////////////////////////////////////////////////////////////////////
/// \brief Execute a function over a sequence of indexes in parallel.
/// Implementation of the Map method.
///
/// \copydetails TExecutorCRTP::Map(F func,ROOT::TSeq<INTEGER> args)
template<class F, class INTEGER, class Cond>
auto TProcessExecutor::MapImpl(F func, ROOT::TSeq<INTEGER> args) -> std::vector<typename std::result_of<F(INTEGER)>::type>
{
   std::vector<INTEGER> vargs(args.size());
   std::copy(args.begin(), args.end(), vargs.begin());
   const auto &reslist = Map(func, vargs);
   return reslist;
}

//////////////////////////////////////////////////////////////////////////
/// \brief Execute a function `nTimes` in parallel (Map) and accumulate the results into a single value (Reduce).
/// \copydetails  ROOT::Internal::TExecutor::MapReduce(F func,unsigned nTimes,R redfunc)
template<class F, class R, class Cond>
auto TProcessExecutor::MapReduce(F func, unsigned nTimes, R redfunc) -> typename std::result_of<F()>::type
{
   using retType = decltype(func());
   //prepare environment
   Reset();
   fTaskType= ETask::kMapRed;

   //fork max(nTimes, fNWorkers) times
   unsigned oldNWorkers = GetPoolSize();
   if (nTimes < oldNWorkers)
      SetNWorkers(nTimes);
   TMPWorkerExecutor<F, void, R> worker(func, redfunc);
   bool ok = Fork(worker);
   SetNWorkers(oldNWorkers);
   if (!ok) {
      std::cerr << "[E][C] Could not fork. Aborting operation\n";
      return retType();
   }

   //give workers their first task
   fNToProcess = nTimes;
   std::vector<retType> reslist;
   reslist.reserve(fNToProcess);
   fNProcessed = Broadcast(MPCode::kExecFunc, fNToProcess);

   //collect results/give workers their next task
   Collect(reslist);

   //clean-up and return
   ReapWorkers();
   fTaskType= ETask::kNoTask;
   return redfunc(reslist);
}

//////////////////////////////////////////////////////////////////////////
/// \brief Execute a function in parallel over the elements of a vector (Map) and accumulate the results into a single value (Reduce).
/// Benefits from partial reduction into `nChunks` intermediate results.
///
/// \copydetails ROOT::Internal::TExecutor::MapReduce(F func,std::vector<T> &args,R redfunc,unsigned nChunks).
template<class F, class T, class R, class Cond>
auto TProcessExecutor::MapReduce(F func, std::vector<T> &args, R redfunc) -> typename std::result_of<F(T)>::type
{

   using retType = decltype(func(args.front()));
   //prepare environment
   Reset();
   fTaskType= ETask::kMapRedWithArg;

   //fork max(args.size(), fNWorkers) times
   unsigned oldNWorkers = GetPoolSize();
   if (args.size() < oldNWorkers)
      SetNWorkers(args.size());
   TMPWorkerExecutor<F, T, R> worker(func, args, redfunc);
   bool ok = Fork(worker);
   SetNWorkers(oldNWorkers);
   if (!ok) {
      std::cerr << "[E][C] Could not fork. Aborting operation\n";
      return decltype(func(args.front()))();
   }

   //give workers their first task
   fNToProcess = args.size();
   std::vector<retType> reslist;
   reslist.reserve(fNToProcess);
   std::vector<unsigned> range(fNToProcess);
   std::iota(range.begin(), range.end(), 0);
   fNProcessed = Broadcast(MPCode::kExecFuncWithArg, range);

   //collect results/give workers their next task
   Collect(reslist);

   ReapWorkers();
   fTaskType= ETask::kNoTask;
   return Reduce(reslist, redfunc);
}

//////////////////////////////////////////////////////////////////////////
/// \brief Execute a function in parallel over the elements of an immutable vector (Map) and accumulate the results into a single value (Reduce).
/// Benefits from partial reduction into `nChunks` intermediate results.
///
/// \copydetails ROOT::Internal::TExecutor::MapReduce(F func,const std::vector<T> &args,R redfunc,unsigned nChunks).
template<class F, class T, class R, class Cond>
auto TProcessExecutor::MapReduce(F func, const std::vector<T> &args, R redfunc) -> typename std::result_of<F(T)>::type
{

   using retType = decltype(func(args.front()));
   //prepare environment
   Reset();
   fTaskType= ETask::kMapRedWithArg;

   //fork max(args.size(), fNWorkers) times
   unsigned oldNWorkers = GetPoolSize();
   if (args.size() < oldNWorkers)
      SetNWorkers(args.size());
   TMPWorkerExecutor<F, T, R> worker(func, args, redfunc);
   bool ok = Fork(worker);
   SetNWorkers(oldNWorkers);
   if (!ok) {
      std::cerr << "[E][C] Could not fork. Aborting operation\n";
      return decltype(func(args.front()))();
   }

   //give workers their first task
   fNToProcess = args.size();
   std::vector<retType> reslist;
   reslist.reserve(fNToProcess);
   std::vector<unsigned> range(fNToProcess);
   std::iota(range.begin(), range.end(), 0);
   fNProcessed = Broadcast(MPCode::kExecFuncWithArg, range);

   //collect results/give workers their next task
   Collect(reslist);

   ReapWorkers();
   fTaskType= ETask::kNoTask;
   return Reduce(reslist, redfunc);
}

//////////////////////////////////////////////////////////////////////////
/// Handle message and reply to the worker
template<class T>
void TProcessExecutor::HandlePoolCode(MPCodeBufPair &msg, TSocket *s, std::vector<T> &reslist)
{
   unsigned code = msg.first;
   if (code == MPCode::kFuncResult) {
      reslist.push_back(std::move(ReadBuffer<T>(msg.second.get())));
      ReplyToFuncResult(s);
   } else if (code == MPCode::kIdling) {
      ReplyToIdle(s);
   } else if(code == MPCode::kProcResult) {
      if(msg.second != nullptr)
         reslist.push_back(std::move(ReadBuffer<T>(msg.second.get())));
      MPSend(s, MPCode::kShutdownOrder);
   } else if(code == MPCode::kProcError) {
      const char *str = ReadBuffer<const char*>(msg.second.get());
      Error("TProcessExecutor::HandlePoolCode", "[E][C] a worker encountered an error: %s\n"
                                         "Continuing execution ignoring these entries.", str);
      ReplyToIdle(s);
      delete [] str;
   } else {
      // UNKNOWN CODE
      Error("TProcessExecutor::HandlePoolCode", "[W][C] unknown code received from server. code=%d", code);
   }
}

//////////////////////////////////////////////////////////////////////////
/// Listen for messages sent by the workers and call the appropriate handler function.
/// TProcessExecutor::HandlePoolCode is called on messages with a code < 1000 and
/// TMPClient::HandleMPCode is called on messages with a code >= 1000.
template<class T>
void TProcessExecutor::Collect(std::vector<T> &reslist)
{
   TMonitor &mon = GetMonitor();
   mon.ActivateAll();
   while (mon.GetActive() > 0) {
      TSocket *s = mon.Select();
      MPCodeBufPair msg = MPRecv(s);
      if (msg.first == MPCode::kRecvError) {
         Error("TProcessExecutor::Collect", "[E][C] Lost connection to a worker");
         Remove(s);
      } else if (msg.first < 1000)
         HandlePoolCode(msg, s, reslist);
      else
         HandleMPCode(msg, s);
   }
}

} // ROOT namespace

#endif
