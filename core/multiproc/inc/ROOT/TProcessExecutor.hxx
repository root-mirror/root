/* @(#)root/multiproc:$Id$ */
// Author: Enrico Guiraud July 2015
// Modified: G Ganis Jan 2017

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
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
#include "TChain.h"
#include "TChainElement.h"
#include "TError.h"
#include "TFileCollection.h"
#include "TFileInfo.h"
#include "THashList.h"
#include "TMPClient.h"
#include "ROOT/TExecutor.hxx"
#include "TMPWorkerExecutor.h"
#include <algorithm> //std::generate
#include <numeric> //std::iota
#include <string>
#include <type_traits> //std::result_of, std::enable_if
#include <functional> //std::reference_wrapper
#include <vector>

namespace ROOT {

   class TProcessExecutor : public TExecutor<TProcessExecutor>, private TMPClient {
   public:
      explicit TProcessExecutor(unsigned nWorkers = 0); //default number of workers is the number of processors
      ~TProcessExecutor() = default;
      //it doesn't make sense for a TProcessExecutor to be copied
      TProcessExecutor(const TProcessExecutor &) = delete;
      TProcessExecutor &operator=(const TProcessExecutor &) = delete;

      // Map
      template<class F, class Cond = noReferenceCond<F>>
      auto Map(F func, unsigned nTimes) -> std::vector<typename std::result_of<F()>::type>;
      /// \cond
      template<class F, class INTEGER, class Cond = noReferenceCond<F, INTEGER>>
      auto Map(F func, ROOT::TSeq<INTEGER> args) -> std::vector<typename std::result_of<F(INTEGER)>::type>;
      template<class F, class T, class Cond = noReferenceCond<F, T>>
      auto Map(F func, std::vector<T> &args) -> std::vector<typename std::result_of<F(T)>::type>;
      /// \endcond
      using TExecutor<TProcessExecutor>::Map;

      template<class F, class R, class Cond = noReferenceCond<F>>
      auto MapReduce(F func, unsigned nTimes, R redfunc, unsigned nChunks) -> typename std::result_of<F()>::type;
      template<class F, class INTEGER, class R, class Cond = noReferenceCond<F, INTEGER>>
      auto MapReduce(F func, ROOT::TSeq<INTEGER> args, R redfunc, unsigned nChunks) -> typename std::result_of<F(INTEGER)>::type;
      // /// \cond doxygen should ignore these methods
      template<class F, class T, class R, class Cond = noReferenceCond<F, T>>
      auto MapReduce(F func, std::initializer_list<T> args, R redfunc, unsigned nChunks) -> typename std::result_of<F(T)>::type;
      template<class F, class T, class R, class Cond = noReferenceCond<F, T>>
      auto MapReduce(F func, std::vector<T> &args, R redfunc, unsigned nChunks) -> typename std::result_of<F(T)>::type;
      using TExecutor<TProcessExecutor>::MapReduce;

      void SetNWorkers(unsigned n)
      {
         TMPClient::SetNWorkers(n);
      }
      unsigned GetNWorkers() const
      {
         return TMPClient::GetNWorkers();
      }

      template<class T, class R> T Reduce(const std::vector<T> &objs, R redfunc);
      using TExecutor<TProcessExecutor>::Reduce;

   protected:
      template<class F, class R, class Cond = noReferenceCond<F>>
      auto Map(F func, unsigned nTimes, R redfunc, unsigned nChunks) -> std::vector<typename std::result_of<F()>::type>;
      template<class F, class INTEGER, class R, class Cond = noReferenceCond<F, INTEGER>>
      auto Map(F func, ROOT::TSeq<INTEGER> args, R redfunc, unsigned nChunks) -> std::vector<typename std::result_of<F(INTEGER)>::type>;
      template<class F, class T, class R, class Cond = noReferenceCond<F, T>>
      auto Map(F func, std::vector<T> &args, R redfunc, unsigned nChunks) -> std::vector<typename std::result_of<F(T)>::type>;
      template<class F, class T, class R, class Cond = noReferenceCond<F, T>>
      auto Map(F func, std::initializer_list<T> args, R redfunc, unsigned nChunks) -> std::vector<typename std::result_of<F(T)>::type>;

   private:
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
         kMapWithArg    ///< a Map method with arguments is being executed
      };

      ETask fTaskType = ETask::kNoTask; ///< the kind of task that is being executed, if any
   };


   /************ TEMPLATE METHODS IMPLEMENTATION ******************/

//////////////////////////////////////////////////////////////////////////
/// Execute func (with no arguments) nTimes in parallel.
/// A vector containg executions' results is returned.
/// Functions that take more than zero arguments can be executed (with
/// fixed arguments) by wrapping them in a lambda or with std::bind.
   template<class F, class Cond>
   auto TProcessExecutor::Map(F func, unsigned nTimes) -> std::vector<typename std::result_of<F()>::type>
   {
      using retType = decltype(func());
      //prepare environment
      Reset();
      fTaskType = ETask::kMap;

      //fork max(nTimes, fNWorkers) times
      unsigned oldNWorkers = GetNWorkers();
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

    template<class F, class R, class Cond>
   auto TProcessExecutor::Map(F func, unsigned nTimes, R redfunc, unsigned nChunks) -> std::vector<typename std::result_of<F()>::type>
   {
      if (nChunks == 0)
      {
         return Map(func, nTimes);
      }

      using retType = decltype(func());
      //prepare environment
      Reset();
      fTaskType = ETask::kMapWithArg;
      unsigned step = (nTimes + nChunks - 1) / nChunks; //ceiling the division
      auto lambda = [&](unsigned int i)
      {
         std::vector<retType> partialResults(step);
         for (unsigned j = 0; j < step && (i + j) < nTimes; j++) {
            partialResults[j] = func();
         }
         return redfunc(partialResults);
      };

      unsigned oldNWorkers = GetNWorkers();
      if (nTimes < oldNWorkers)
         SetNWorkers(nTimes);
      fNToProcess = nChunks;
      std::vector<retType> reslist;
      reslist.reserve(fNToProcess);
      std::vector<unsigned> range(fNToProcess);
      std::iota(range.begin(), range.end(), 0);
      TMPWorkerExecutor<decltype(lambda), unsigned> worker(lambda, range);
      bool ok = Fork(worker);
      SetNWorkers(oldNWorkers);
      if (!ok)
      {
         Error("TProcessExecutor::Map", "[E][C] Could not fork. Aborting operation.");
         return std::vector<retType>();
      }

      //give out tasks
      fNProcessed = Broadcast(MPCode::kExecFuncWithArg, range);

      //collect results, give out other tasks if needed
      Collect(reslist);

      //clean-up and return
      ReapWorkers();
      fTaskType = ETask::kNoTask;
      return reslist;
   }

// tell doxygen to ignore this (\endcond closes the statement)
/// \cond

// actual implementation of the Map method. all other calls with arguments eventually
// call this one
   template<class F, class T, class Cond>
   auto TProcessExecutor::Map(F func, std::vector<T> &args) -> std::vector<typename std::result_of<F(T)>::type>
   {
      //check whether func is callable
      using retType = decltype(func(args.front()));
      //prepare environment
      Reset();
      fTaskType = ETask::kMapWithArg;

      //fork max(args.size(), fNWorkers) times
      //N.B. from this point onwards, args is filled with undefined (but valid) values, since TMPWorkerExecutor moved its content away
      unsigned oldNWorkers = GetNWorkers();
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

   template<class F, class T, class R, class Cond>
   auto TProcessExecutor::Map(F func, std::vector<T> &args, R redfunc, unsigned nChunks) -> std::vector<typename std::result_of<F(T)>::type>
   {

      if (nChunks == 0) {
         return Map(func, args);
      }
      // //check whether func is callable
      using retType = decltype(func(args.front()));
       //prepare environment
      Reset();
      fTaskType = ETask::kMapWithArg;

      fNToProcess = nChunks;
      std::vector<retType> reslist;
      reslist.reserve(fNToProcess);
      unsigned step = (args.size() + nChunks - 1) / nChunks; //ceiling the division

      auto lambda = [&](int i) {
         std::vector<retType> partialResults(step);
         for (unsigned j = 0; j < step && (i + j) < args.size(); j++) {
            partialResults[j] = func(args[i + j]);
         }
         return redfunc(partialResults);
      };

      //fork max(args.size(), fNWorkers) times
      //N.B. from this point onwards, args is filled with undefined (but valid) values, since TMPWorkerExecutor moved its content away
      unsigned oldNWorkers = GetNWorkers();
      if (nChunks < oldNWorkers)
         SetNWorkers(nChunks);
      TMPWorkerExecutor<decltype(lambda), T> worker(lambda, args);
      bool ok = Fork(worker);
      SetNWorkers(oldNWorkers);
      if (!ok) {
         Error("TProcessExecutor::Map", "[E][C] Could not fork. Aborting operation.");
         return std::vector<retType>();
      }

      ROOT::TSeq<unsigned> chunkIt(0, args.size(), step);
      std::vector<unsigned> vargs(args.size());
      fNProcessed = Broadcast(MPCode::kExecFuncWithArg, vargs);

      //collect results, give out other tasks if needed
      Collect(reslist);
      //clean-up and return
      ReapWorkers();
      fTaskType = ETask::kNoTask;
      return reslist;
   }


   template<class F, class INTEGER, class Cond>
   auto TProcessExecutor::Map(F func, ROOT::TSeq<INTEGER> args) -> std::vector<typename std::result_of<F(INTEGER)>::type>
   {
      std::vector<INTEGER> vargs(args.size());
      std::copy(args.begin(), args.end(), vargs.begin());
      const auto &reslist = Map(func, vargs);
      return reslist;
   }

   template<class F, class INTEGER, class R, class Cond>
   auto TProcessExecutor::Map(F func, ROOT::TSeq<INTEGER> args, R redfunc, unsigned nChunks) -> std::vector<typename std::result_of<F(INTEGER)>::type> {
      std::vector<INTEGER> vargs(args.size());
      std::copy(args.begin(), args.end(), vargs.begin());
      const auto &reslist = Map(func, vargs, redfunc, nChunks);
      return reslist;
   }

   template<class F, class T, class R, class Cond>
   auto TProcessExecutor::Map(F func, std::initializer_list<T> args, R redfunc, unsigned nChunks) -> std::vector<typename std::result_of<F(T)>::type> {
      std::vector<T> vargs(std::move(args));
      const auto &reslist = Map(func, vargs, redfunc, nChunks);
      return reslist;
   }

// tell doxygen to stop ignoring code
/// \endcond
// //////////////////////////////////////////////////////////////////////////
// /// This method behaves just like Map, but an additional redfunc function
// /// must be provided. redfunc is applied to the vector Map would return and
// /// must return the same type as func. In practice, redfunc can be used to
// /// "squash" the vector returned by Map into a single object by merging,
// /// adding, mixing the elements of the vector.
   template<class F, class R, class Cond>
   auto TProcessExecutor::MapReduce(F func, unsigned nTimes, R redfunc, unsigned nChunks) -> typename std::result_of<F()>::type {
      return Reduce(Map(func, nTimes, redfunc, nChunks), redfunc);
   }

/// \cond doxygen should ignore these methods
   template<class F, class INTEGER, class R, class Cond>
   auto TProcessExecutor::MapReduce(F func, ROOT::TSeq<INTEGER> args, R redfunc, unsigned nChunks) -> typename std::result_of<F(INTEGER)>::type {
      return Reduce(Map(func, args, redfunc, nChunks), redfunc);
   }

   template<class F, class T, class R, class Cond>
   auto TProcessExecutor::MapReduce(F func, std::initializer_list<T> args, R redfunc, unsigned nChunks) -> typename std::result_of<F(T)>::type {
      return Reduce(Map(func, args, redfunc, nChunks), redfunc);
   }

   template<class F, class T, class R, class Cond>
   auto TProcessExecutor::MapReduce(F func, std::vector<T> &args, R redfunc, unsigned nChunks) -> typename std::result_of<F(T)>::type {
      return Reduce(Map(func, args, redfunc, nChunks), redfunc);
   }

   template<class T, class R>
   T TProcessExecutor::Reduce(const std::vector<T> &objs, R redfunc)
   {
      // check we can apply reduce to objs
      static_assert(std::is_same<decltype(redfunc(objs)), T>::value, "redfunc does not have the correct signature");
      return redfunc(objs);
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
      } else if (code == MPCode::kProcResult) {
         if (msg.second != nullptr)
            reslist.push_back(std::move(ReadBuffer<T>(msg.second.get())));
         MPSend(s, MPCode::kShutdownOrder);
      } else if (code == MPCode::kProcError) {
         const char *str = ReadBuffer<const char *>(msg.second.get());
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
