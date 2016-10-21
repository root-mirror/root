// @(#)root/thread:$Id$
// Author: Xavier Valls March 2016

/*************************************************************************
 * Copyright (C) 1995-2006, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TThreadExecutor
#define ROOT_TThreadExecutor

#include "RConfigure.h"
// exclude in case ROOT does not have IMT support
#ifdef R__USE_IMT

// exclude in case ROOT does not have IMT support
#ifndef R__USE_IMT
// No need to error out for dictionaries.
# if !defined(__ROOTCLING__) && !defined(G__DICTIONARY)
#  error "Cannot use ROOT::TThreadExecutor without defining R__USE_IMT."
# endif
#else
#include "ROOT/TExecutor.hxx"
#include <numeric>
#include <functional>

namespace tbb { class task_scheduler_init;}

namespace ROOT {

class TThreadExecutor: public TExecutor<TThreadExecutor> {
public:
   explicit TThreadExecutor();

   explicit TThreadExecutor(size_t nThreads);
   
   TThreadExecutor(TThreadExecutor &) = delete;
   TThreadExecutor & operator=(TThreadExecutor &) = delete;

   ~TThreadExecutor();

   template<class F, class Cond = noReferenceCond<F>>
   auto Map(F func, unsigned nTimes) -> std::vector<typename std::result_of<F()>::type>;
   /// \cond
   template<class F, class INTEGER, class Cond = noReferenceCond<F, INTEGER>>
   auto Map(F func, ROOT::TSeq<INTEGER> args) -> std::vector<typename std::result_of<F(INTEGER)>::type>;
   template<class F, class T, class Cond = noReferenceCond<F, T>>
   auto Map(F func, std::vector<T> &args) -> std::vector<typename std::result_of<F(T)>::type>;
   // / \endcond
   using TExecutor<TThreadExecutor>::Map;
   
   template<class T, class BINARYOP> auto Reduce(const std::vector<T> &objs, BINARYOP redfunc) -> decltype(redfunc(objs.front(), objs.front()));
   using TExecutor<TThreadExecutor>::Reduce;

private:
    void _parallelFor(unsigned start, unsigned end, const std::function<void(unsigned int i)> &f);
    tbb::task_scheduler_init *fInitTBB;
};

/************ TEMPLATE METHODS IMPLEMENTATION ******************/

//////////////////////////////////////////////////////////////////////////
/// Execute func (with no arguments) nTimes in parallel.
/// A vector containg executions' results is returned.
/// Functions that take more than zero arguments can be executed (with
/// fixed arguments) by wrapping them in a lambda or with std::bind.
template<class F, class Cond>
auto TThreadExecutor::Map(F func, unsigned nTimes) -> std::vector<typename std::result_of<F()>::type>
{
   using retType = decltype(func());
   std::vector<retType> reslist(nTimes);

   auto lambda = [&](unsigned int i){reslist[i] = func();};
   _parallelFor(0U, nTimes, lambda);

   return reslist;
}

template<class F, class INTEGER, class Cond>
auto TThreadExecutor::Map(F func, ROOT::TSeq<INTEGER> args) -> std::vector<typename std::result_of<F(INTEGER)>::type>
{
   unsigned start = *args.begin();
   unsigned end = *args.end();
   using retType = decltype(func(start));
   std::vector<retType> reslist(end-start);

   auto lambda = [&](unsigned int i){reslist[i] = func(i);};
   _parallelFor(start, end, lambda);

   return reslist;
}

// tell doxygen to ignore this (\endcond closes the statement)
/// \cond

// actual implementation of the Map method. all other calls with arguments eventually
// call this one
template<class F, class T, class Cond>
auto TThreadExecutor::Map(F func, std::vector<T> &args) -> std::vector<typename std::result_of<F(T)>::type>
{
   // //check whether func is callable
   using retType = decltype(func(args.front()));

   unsigned int fNToProcess = args.size();
   std::vector<retType> reslist(fNToProcess);

   auto lambda = [&](unsigned int i){reslist[i] = func(args[i]);};
   _parallelFor(0U, fNToProcess, lambda);

   return reslist;
}

// // tell doxygen to stop ignoring code
// /// \endcond

template<class T, class BINARYOP>
auto TThreadExecutor::Reduce(const std::vector<T> &objs, BINARYOP redfunc) -> decltype(redfunc(objs.front(), objs.front()))
{
   // check we can apply reduce to objs
   static_assert(std::is_same<decltype(redfunc(objs.front(), objs.front())), T>::value, "redfunc does not have the correct signature");   
   return std::accumulate(objs.begin(), objs.end(), T{}, redfunc);
}

} // namespace ROOT

#endif   // R__USE_IMT
#endif
