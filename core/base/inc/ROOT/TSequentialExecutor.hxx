// @(#)root/thread:$Id$
// Author: Xavier Valls November 2017

/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/
#ifndef ROOT_TSequentialExecutor
#define ROOT_TSequentialExecutor

#include "ROOT/EExecutionPolicy.hxx"
#include "ROOT/TExecutorCRTP.hxx"
#include "ROOT/TSeq.hxx"

#include <initializer_list>
#include <numeric> //std::accumulate
#include <type_traits> //std::enable_if, std::result_of
#include <utility> //std::move
#include <vector>

namespace ROOT {

   class TSequentialExecutor: public TExecutorCRTP<TSequentialExecutor> {
      friend TExecutorCRTP;
   public:

      TSequentialExecutor() = default;
      TSequentialExecutor(const TSequentialExecutor &) = delete;
      TSequentialExecutor &operator=(const TSequentialExecutor &) = delete;

      // Foreach
      //
      template<class F>
      void Foreach(F func, unsigned nTimes);
      template<class F, class INTEGER>
      void Foreach(F func, ROOT::TSeq<INTEGER> args);
      template<class F, class T>
      void Foreach(F func, std::initializer_list<T> args);
      template<class F, class T>
      void Foreach(F func, std::vector<T> &args);
      template<class F, class T>
      void Foreach(F func, const std::vector<T> &args);

      // Map
      //
      using TExecutorCRTP<TSequentialExecutor>::Map;

      // MapReduce
      // the late return types also check at compile-time whether redfunc is compatible with func,
      // other than checking that func is compatible with the type of arguments.
      // a static_assert check in TSequentialExecutor::Reduce is used to check that redfunc is compatible with the type returned by func
      using TExecutorCRTP<TSequentialExecutor>::MapReduce;

      // Reduce
      //
      using TExecutorCRTP<TSequentialExecutor>::Reduce;

      //////////////////////////////////////////////////////////////////////////
      /// \brief Return the number of workers in the sequential executor: a single one.
      ///
      /// \return The number of workers in the pool, one.
      unsigned GetPoolSize() const { return 1u; }

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
   };

   /************ TEMPLATE METHODS IMPLEMENTATION ******************/

   //////////////////////////////////////////////////////////////////////////
   /// \brief Execute a function without arguments several times, dividing the execution in nChunks.
   ///
   /// \param func Function to be executed.
   /// \param nTimes Number of times function should be called.
   template<class F>
   void TSequentialExecutor::Foreach(F func, unsigned nTimes) {
      for (auto i = 0U; i < nTimes; ++i) func();
   }

   //////////////////////////////////////////////////////////////////////////
   /// \brief Execute a function over a sequence of indexes, dividing the execution in nChunks.
   ///
   /// \param func Function to be executed. Must take an element of the sequence passed assecond argument as a parameter.
   /// \param args Sequence of indexes to execute `func` on.
   template<class F, class INTEGER>
   void TSequentialExecutor::Foreach(F func, ROOT::TSeq<INTEGER> args) {
      for(auto i : args) func(i);
   }

   //////////////////////////////////////////////////////////////////////////
   /// \brief Execute a function over the elements of an initializer_list, dividing the execution in nChunks.
   ///
   /// \param func Function to be executed on the elements of the initializer_list passed as second parameter.
   /// \param args initializer_list for a vector to apply `func` on.
   template<class F, class T>
   void TSequentialExecutor::Foreach(F func, std::initializer_list<T> args) {
      std::vector<T> vargs(std::move(args));
      Foreach(func, vargs);
   }

   //////////////////////////////////////////////////////////////////////////
   /// \brief Execute a function over the elements of a vector, dividing the execution in nChunks.
   ///
   /// \param func Function to be executed on the elements of the vector passed as second parameter.
   /// \param args Vector of elements passed as an argument to `func`.
   template<class F, class T>
   void TSequentialExecutor::Foreach(F func, std::vector<T> &args) {
      for(auto &&arg: args) {
         func(arg);
      }
   }

   //////////////////////////////////////////////////////////////////////////
   /// \brief Execute a function over the elements of an immutable vector, dividing the execution in nChunks.
   ///
   /// \param func Function to be executed on the elements of the immutable vector passed as second parameter.
   /// \param args Immutable vector of elements passed as an argument to `func`.
   template<class F, class T>
   void TSequentialExecutor::Foreach(F func, const std::vector<T> &args) {
      for(auto &&arg: args) {
         func(arg);
      }
   }

   //////////////////////////////////////////////////////////////////////////
   /// \brief Execute a function without arguments several times.
   /// Implementation of the Map method.
   ///
   /// \copydetails TExecutorCRTP::Map(F func,unsigned nTimes)
   template<class F, class Cond>
   auto TSequentialExecutor::MapImpl(F func, unsigned nTimes) -> std::vector<typename std::result_of<F()>::type> {
      using retType = decltype(func());
      std::vector<retType> reslist;
      reslist.reserve(nTimes);
      while(reslist.size() < nTimes) {
         reslist.emplace_back(func());
      }
      return reslist;
   }

   //////////////////////////////////////////////////////////////////////////
   /// \brief Execute a function over a sequence of indexes.
   /// Implementation of the Map method.
   ///
   /// \copydetails TExecutorCRTP::Map(F func,ROOT::TSeq<INTEGER> args)
   template<class F, class INTEGER, class Cond>
   auto TSequentialExecutor::MapImpl(F func, ROOT::TSeq<INTEGER> args) -> std::vector<typename std::result_of<F(INTEGER)>::type> {
      using retType = decltype(func(*args.begin()));
      std::vector<retType> reslist;
      reslist.reserve(args.size());
      for(auto i: args)
         reslist.emplace_back(func(i));
      return reslist;
   }

   //////////////////////////////////////////////////////////////////////////
   /// \brief Execute a function over the elements of a vector in parallel
   /// Implementation of the Map method.
   ///
   /// \copydetails TExecutorCRTP::Map(F func,std::vector<T> &args)
   template<class F, class T, class Cond>
   auto TSequentialExecutor::MapImpl(F func, std::vector<T> &args) -> std::vector<typename std::result_of<F(T)>::type> {
      // //check whether func is callable
      using retType = decltype(func(args.front()));
      std::vector<retType> reslist;
      reslist.reserve(args.size());
      for(auto &&arg: args) {
         reslist.emplace_back(func(arg));
      }
      return reslist;
   }

   //////////////////////////////////////////////////////////////////////////
   /// \brief Execute a function over the elements of an immutable vector.
   /// Implementation of the Map method.
   ///
   /// \copydetails TExecutorCRTP::Map(F func,const std::vector<T> &args)
   template<class F, class T, class Cond>
   auto TSequentialExecutor::MapImpl(F func, const std::vector<T> &args) -> std::vector<typename std::result_of<F(T)>::type> {
      // //check whether func is callable
      using retType = decltype(func(args.front()));
      std::vector<retType> reslist;
      reslist.reserve(args.size());
      for(auto &&arg: args) {
         reslist.emplace_back(func(arg));
      }
      return reslist;
   }

} // namespace ROOT
#endif
