// @(#)root/thread:$Id$
// Author: Xavier Valls March 2016

/*************************************************************************
 * Copyright (C) 1995-2006, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TPOOL
#define ROOT_TPOOL

#include "TCollection.h"
#include "ROOT/TSeq.h"


template<class subc>
class TPool {
public:
   explicit TPool() = default;
   explicit TPool(size_t /* nThreads */ ){};

   template< class F, class... T>
   using noReferenceCond = typename std::enable_if<"Function can't return a reference" && !(std::is_reference<typename std::result_of<F(T...)>::type>::value)>::type;

   // // Map
   // //these late return types allow for a compile-time check of compatibility between function signatures and args,
   // //and a compile-time check that the argument list implements a front() method (all STL sequence containers have it)
   template<class F, class Cond = noReferenceCond<F>>
   auto Map(F func, unsigned nTimes) -> std::vector<typename std::result_of<F()>::type>;
   template<class F, class T, class Cond = noReferenceCond<F, typename T::value_type>>
   auto Map(F func, T &args) -> std::vector < decltype(++(args.begin()), args.end(), func(args.front()))>;
   // /// \cond doxygen should ignore these methods
   template<class F, class INTEGER, class Cond = noReferenceCond<F, INTEGER>>
   auto Map(F func, ROOT::TSeq<INTEGER> args) -> std::vector<typename std::result_of<F(INTEGER)>::type>;
   template<class F, class Cond = noReferenceCond<F, TObject*>>
   auto Map(F func, TCollection &args) -> std::vector<typename std::result_of<F(TObject*)>::type>;
   template<class F, class T, class Cond = noReferenceCond<F, T>>
   auto Map(F func, std::initializer_list<T> args) -> std::vector<typename std::result_of<F(T)>::type>;
   template<class F, class T, class Cond = noReferenceCond<F, T>>
   auto Map(F func, std::vector<T> &args) -> std::vector<typename std::result_of<F(T)>::type>;
   // // // / \endcond

   // // MapReduce
   // // the late return types also check at compile-time whether redfunc is compatible with func,
   // // other than checking that func is compatible with the type of arguments.
   // // a static_assert check in TPool<subc>::Reduce is used to check that redfunc is compatible with the type returned by func
   template<class F, class R, class Cond = noReferenceCond<F>>
   auto MapReduce(F func, unsigned nTimes, R redfunc) -> typename std::result_of<F()>::type;
   template<class F, class INTEGER, class R, class Cond = noReferenceCond<F, INTEGER>>
   auto MapReduce(F func, ROOT::TSeq<INTEGER> args, R redfunc) -> typename std::result_of<F(INTEGER)>::type;
   // /// \cond doxygen should ignore these methods
   template<class F, class T, class R, class Cond = noReferenceCond<F, T>>
   auto MapReduce(F func, std::initializer_list<T> args, R redfunc) -> typename std::result_of<F(T)>::type;
   template<class F, class T, class R, class Cond = noReferenceCond<F, T>>
   auto MapReduce(F func, std::vector<T> &args, R redfunc) -> typename std::result_of<F(T)>::type;
   // /// \endcond

protected:
   template<class T, class R> auto Reduce(const std::vector<T> &objs, R redfunc) -> decltype(redfunc(objs));
   template<class T, class BINARYOP> auto Reduce(const std::vector<T> &objs, BINARYOP redfunc) -> decltype(redfunc(objs.front(), objs.front()));

private:
  inline subc & Derived()
  {
    return *static_cast<subc*>(this);
  }
};

//////////////////////////////////////////////////////////////////////////
/// Execute func (with no arguments) nTimes in parallel.
/// A vector containg executions' results is returned.
/// Functions that take more than zero arguments can be executed (with
/// fixed arguments) by wrapping them in a lambda or with std::bind.
template<class subc> template<class F, class Cond>
auto TPool<subc>::Map(F func, unsigned nTimes) -> std::vector<typename std::result_of<F()>::type>
{
   return Derived().Map(func, nTimes);
}

// //////////////////////////////////////////////////////////////////////////
// /// Execute func in parallel distributing the elements of the args collection between the workers.
// /// See class description for the valid types of collections and containers that can be used.
// /// A vector containing each execution's result is returned. The user is responsible of deleting
// /// objects that might be created upon the execution of func, returned objects included.
// /// **Note:** the collection of arguments is modified by Map and should be considered empty or otherwise
// /// invalidated after Map's execution (std::move might be applied to it).

// tell doxygen to ignore this (\endcond closes the statement)
/// \cond
template<class subc> template<class F, class INTEGER, class Cond>
auto TPool<subc>::Map(F func, ROOT::TSeq<INTEGER> args) -> std::vector<typename std::result_of<F(INTEGER)>::type>
{
  return Derived().Map(func, args);
}

template<class subc> template<class F, class T, class Cond>
auto TPool<subc>::Map(F func, std::initializer_list<T> args) -> std::vector<typename std::result_of<F(T)>::type>
{
   std::vector<T> vargs(std::move(args));
   const auto &reslist = Map(func, vargs);
   return reslist;
}

// actual implementation of the Map method. all other calls with arguments eventually
// call this one

template<class subc> template<class F, class T, class Cond>
auto TPool<subc>::Map(F func, std::vector<T> &args) -> std::vector<typename std::result_of<F(T)>::type>
{
   return Derived().Map(func, args);
}

// //////////////////////////////////////////////////////////////////////////
// /// This method behaves just like Map, but an additional redfunc function
// /// must be provided. redfunc is applied to the vector Map would return and
// /// must return the same type as func. In practice, redfunc can be used to
// /// "squash" the vector returned by Map into a single object by merging,
// /// adding, mixing the elements of the vector.
template<class subc> template<class F, class R, class Cond>
auto TPool<subc>::MapReduce(F func, unsigned nTimes, R redfunc) -> typename std::result_of<F()>::type
{
   return Reduce(Map(func, nTimes), redfunc);
}

//////////////////////////////////////////////////////////////////////////
/// This method behaves just like Map, but an additional redfunc function
/// must be provided. redfunc is applied to the vector Map would return and
/// must return the same type as func. In practice, redfunc can be used to
/// "squash" the vector returned by Map into a single object by merging,
/// adding, mixing the elements of the vector.

/// \cond doxygen should ignore these methods
template<class subc> template<class F, class INTEGER, class R, class Cond>
auto TPool<subc>::MapReduce(F func, ROOT::TSeq<INTEGER> args, R redfunc) -> typename std::result_of<F(INTEGER)>::type
{
  return Reduce(Map(func, args), redfunc);
}

template<class subc> template<class F, class T, class R, class Cond>
auto TPool<subc>::MapReduce(F func, std::initializer_list<T> args, R redfunc) -> typename std::result_of<F(T)>::type
{
   return Reduce(Map(func, args), redfunc);
}

template<class subc> template<class F, class T, class R, class Cond>
auto TPool<subc>::MapReduce(F func, std::vector<T> &args, R redfunc) -> typename std::result_of<F(T)>::type
{
   return Reduce(Derived().Map(func, args), redfunc);
}

/// \endcond

/// Check that redfunc has the right signature and call it on objs
template<class subc> template<class T, class BINARYOP>
auto TPool<subc>::Reduce(const std::vector<T> &objs, BINARYOP redfunc) -> decltype(redfunc(objs.front(), objs.front()))
{
   // check we can apply reduce to objs
   static_assert(std::is_same<decltype(redfunc(objs.front(), objs.front())), T>::value, "redfunc does not have the correct signature");
   return Derived().Reduce(objs, redfunc);
}

template<class subc> template<class T, class R>
auto TPool<subc>::Reduce(const std::vector<T> &objs, R redfunc) -> decltype(redfunc(objs))
{
   // check we can apply reduce to objs
   static_assert(std::is_same<decltype(redfunc(objs)), T>::value, "redfunc does not have the correct signature");
   return redfunc(objs);
}

#endif
