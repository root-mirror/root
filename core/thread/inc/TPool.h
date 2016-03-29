#ifndef TPOOL_H_
#define TPOOL_H_

#include "TCollection.h"
#include "TObjArray.h"
#include "ROOT/TSeq.h"

template<class subc>
class TPool {
public:
   explicit TPool(){};
   explicit TPool(size_t nThreads){};

   // // Map
   // //these late return types allow for a compile-time check of compatibility between function signatures and args,
   // //and a compile-time check that the argument list implements a front() method (all STL sequence containers have it)
   template<class F> auto Map(F func, unsigned nTimes) -> std::vector<decltype(func())>;
   template<class F, class T> auto Map(F func, T &args) -> std::vector < decltype(++(args.begin()), args.end(), func(args.front())) >;
   // /// \cond doxygen should ignore these methods
   template<class F, class INTEGER> auto Map(F func, ROOT::TSeq<INTEGER> args) -> std::vector<decltype(func(args.front()))>;
   template<class F> TObjArray Map(F func, TCollection &args);
   template<class F, class T> auto Map(F func, std::initializer_list<T> args) -> std::vector<decltype(func(*args.begin()))>;
   template<class F, class T> auto Map(F func, std::vector<T> &args) -> std::vector<decltype(func(args.front()))>;
   // // // / \endcond

   // // MapReduce
   // // the late return types also check at compile-time whether redfunc is compatible with func,
   // // other than checking that func is compatible with the type of arguments.
   // // a static_assert check in TPool<subc>::Reduce is used to check that redfunc is compatible with the type returned by func
   template<class F, class T, class R> auto MapReduce(F func, T*args, int nEvents, R redfunc) -> decltype(func(args, 0));
   template<class F, class R> auto MapReduce(F func, unsigned nTimes, R redfunc) -> decltype(func());
   template<class F, class INTEGER, class R> auto MapReduce(F func, ROOT::TSeq<INTEGER> args, R redfunc) -> decltype(func(args.front()));
   template<class F, class T, class R> auto MapReduce(F func, T &args, R redfunc) -> decltype(++(args.begin()), args.end(), func(args.front()));
   // /// \cond doxygen should ignore these methods
   template<class F, class R> auto MapReduce(F func, TCollection &args, R redfunc) -> decltype(func(nullptr));
   template<class F, class T, class R> auto MapReduce(F func, std::initializer_list<T> args, R redfunc) -> decltype(func(*args.begin()));
   template<class F, class T, class R> auto MapReduce(F func, std::vector<T> &args, R redfunc) -> decltype(func(args.front()));
   // /// \endcond

protected:
   template<class T, class R> auto Reduce(const std::vector<T> &objs, R redfunc) -> decltype(redfunc(objs));
   template<class T, class BINARYOP> auto Reduce(const std::vector<T> &objs, BINARYOP redfunc) -> decltype(redfunc(objs.front(), objs.front()));
};

//////////////////////////////////////////////////////////////////////////
/// Execute func (with no arguments) nTimes in parallel.
/// A vector containg executions' results is returned.
/// Functions that take more than zero arguments can be executed (with
/// fixed arguments) by wrapping them in a lambda or with std::bind.
template<class subc> template<class F>
auto TPool<subc>::Map(F func, unsigned nTimes) -> std::vector<decltype(func())>
{
   return static_cast<subc*>(this)->Map(func, nTimes);
}

// //////////////////////////////////////////////////////////////////////////
// /// Execute func in parallel distributing the elements of the args collection between the workers.
// /// See class description for the valid types of collections and containers that can be used.
// /// A vector containing each execution's result is returned. The user is responsible of deleting
// /// objects that might be created upon the execution of func, returned objects included.
// /// **Note:** the collection of arguments is modified by Map and should be considered empty or otherwise
// /// invalidated after Map's execution (std::move might be applied to it).
template<class subc> template<class F, class T>
auto TPool<subc>::Map(F func, T &args) -> std::vector < decltype(++(args.begin()), args.end(), func(args.front())) >
{
   std::vector<typename T::value_type> vargs(
      std::make_move_iterator(std::begin(args)),
      std::make_move_iterator(std::end(args))
   );
   const auto &reslist = static_cast<subc*>(this)->Map(func, vargs);
   return reslist;
}


// tell doxygen to ignore this (\endcond closes the statement)
/// \cond
template<class subc> template<class F, class INTEGER>
auto TPool<subc>::Map(F func, ROOT::TSeq<INTEGER> args) -> std::vector<decltype(func(args.front()))>
{
  return static_cast<subc*>(this)->Map(func, args);
}

template<class subc> template<class F>
TObjArray TPool<subc>::Map(F func, TCollection &args)
{
   // check the function returns something from which we can build a TObject*
   static_assert(std::is_constructible<TObject *, typename std::result_of<F(TObject *)>::type>::value,
                 "func should return a pointer to TObject or derived classes");

   //build vector with same elements as args
   std::vector<TObject *> vargs(args.GetSize());
   auto it = vargs.begin();
   for (auto o : args) {
      *it = o;
      ++it;
   }

   //call Map
   const auto &reslist = static_cast<subc*>(this)->Map(func, vargs);

   //build TObjArray with same elements as reslist
   TObjArray resarray;
   for (const auto &res : reslist)
      resarray.Add(res);
   return resarray;
}

template<class subc> template<class F, class T>
auto TPool<subc>::Map(F func, std::initializer_list<T> args) -> std::vector<decltype(func(*args.begin()))>
{
   std::vector<T> vargs(std::move(args));
   const auto &reslist = Map(func, vargs);
   return reslist;
}

// actual implementation of the Map method. all other calls with arguments eventually
// call this one

template<class subc> template<class F, class T>
auto TPool<subc>::Map(F func, std::vector<T> &args) -> std::vector<decltype(func(args.front()))>
{
   return static_cast<subc*>(this)->Map(func, args);
}

// //////////////////////////////////////////////////////////////////////////
// /// This method behaves just like Map, but an additional redfunc function
// /// must be provided. redfunc is applied to the vector Map would return and
// /// must return the same type as func. In practice, redfunc can be used to
// /// "squash" the vector returned by Map into a single object by merging,
// /// adding, mixing the elements of the vector.
template<class subc> template<class F, class R>
auto TPool<subc>::MapReduce(F func, unsigned nTimes, R redfunc) -> decltype(func())
{
   return Reduce(Map(func, nTimes), redfunc);
}

//////////////////////////////////////////////////////////////////////////
/// This method behaves just like Map, but an additional redfunc function
/// must be provided. redfunc is applied to the vector Map would return and
/// must return the same type as func. In practice, redfunc can be used to
/// "squash" the vector returned by Map into a single object by merging,
/// adding, mixing the elements of the vector.
template<class subc> template<class F, class T, class R>
auto TPool<subc>::MapReduce(F func, T &args, R redfunc) -> decltype(++(args.begin()), args.end(), func(args.front()))
{
   return Reduce(static_cast<subc*>(this)->Map(func, args), redfunc);
}

/// \cond doxygen should ignore these methods

template<class subc> template<class F, class INTEGER, class R>
auto TPool<subc>::MapReduce(F func, ROOT::TSeq<INTEGER> args, R redfunc) -> decltype(func(args.front()))
{
  return Reduce(Map(func, args), redfunc);
}

template<class subc> template<class F, class R>
auto TPool<subc>::MapReduce(F func, TCollection &args, R redfunc) -> decltype(func(nullptr))
{
   //build vector with same elements as args
   std::vector<TObject *> vargs(args.GetSize());
   auto it = vargs.begin();
   for (auto o : args) {
      *it = o;
      ++it;
   }

   //call MapReduce
   auto res = MapReduce(func, vargs, redfunc); //TODO useless copying by value here, but it looks like the return type of this MapReduce is a reference otherwise
   return res;
}

template<class subc> template<class F, class T, class R>
auto TPool<subc>::MapReduce(F func, std::initializer_list<T> args, R redfunc) -> decltype(func(*args.begin()))
{
   return Reduce(Map(func, args), redfunc);
}

template<class subc> template<class F, class T, class R>
auto TPool<subc>::MapReduce(F func, std::vector<T> &args, R redfunc) -> decltype(func(args.front()))
{
   return Reduce(static_cast<subc*>(this)->Map(func, args), redfunc);
}

/// \endcond

/// Check that redfunc has the right signature and call it on objs
template<class subc> template<class T, class BINARYOP>
auto TPool<subc>::Reduce(const std::vector<T> &objs, BINARYOP redfunc) -> decltype(redfunc(objs.front(), objs.front()))
{
   // check we can apply reduce to objs
   static_assert(std::is_same<decltype(redfunc(objs.front(), objs.front())), T>::value, "redfunc does not have the correct signature");
   return static_cast<subc*>(this)->Reduce(objs, redfunc);
}

template<class subc> template<class T, class R>
auto TPool<subc>::Reduce(const std::vector<T> &objs, R redfunc) -> decltype(redfunc(objs))
{
   // check we can apply reduce to objs
   static_assert(std::is_same<decltype(redfunc(objs)), T>::value, "redfunc does not have the correct signature");
   return redfunc(objs);
}

#endif
