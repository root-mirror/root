#include "ROOT/TypeTraits.hxx"

#include <functional>
#include <memory>
#include <type_traits>
#include <vector>
#include <string>

#include "gtest/gtest.h"

using namespace ROOT::TypeTraits;

TEST(TypeTraits, TypeList)
{
   static_assert(TypeList<>::list_size == 0, "");
   static_assert(TypeList<void>::list_size == 1, "");
   static_assert((TypeList<int, double *>::list_size) == 2, "");
}

TEST(TypeTraits, TakeFirstType)
{
   ::testing::StaticAssertTypeEq<TakeFirstType_t<int, void, double *>, int>();
   ::testing::StaticAssertTypeEq<TakeFirstType_t<TypeList<int, int>, TypeList<void>>, TypeList<int, int>>();
}

TEST(TypeTraits, RemoveFirst)
{
   ::testing::StaticAssertTypeEq<RemoveFirst_t<int, void, double *>, TypeList<void, double *>>();
   ::testing::StaticAssertTypeEq<RemoveFirst_t<TypeList<int, int>, TypeList<void>>, TypeList<TypeList<void>>>();
}

TEST(TypeTraits, TakeFirstParameter)
{
   ::testing::StaticAssertTypeEq<TakeFirstParameter_t<TypeList<int, void>>, int>();
   ::testing::StaticAssertTypeEq<TakeFirstParameter_t<std::tuple<void>>, void>();
}

TEST(TypeTraits, RemoveFirstParameter)
{
   ::testing::StaticAssertTypeEq<RemoveFirstParameter_t<TypeList<int, void>>, TypeList<void>>();
   ::testing::StaticAssertTypeEq<RemoveFirstParameter_t<std::tuple<void>>, std::tuple<>>();
}

/******** helper objects ***********/
struct Dummy {
};

int freeFun1(int)
{
   return 0;
}

Dummy freeFun2(Dummy &, int *, const std::vector<Dummy> &)
{
   return Dummy();
}

struct Functor1 {
   void operator()(int) { return; }
};

struct Functor2 {
   Dummy operator()(Dummy &, int *, const std::vector<Dummy> &) { return Dummy(); }
};
/***********************************/

TEST(TypeTraits, CallableTraits)
{
   // check CallableTraits with one arg, two args, return type, nodecaytypes
   // free function
   ::testing::StaticAssertTypeEq<TypeList<int>, CallableTraits<decltype(freeFun1)>::arg_types>();
   ::testing::StaticAssertTypeEq<TypeList<int>, CallableTraits<decltype(freeFun1)>::arg_types_nodecay>();
   ::testing::StaticAssertTypeEq<int, CallableTraits<decltype(freeFun1)>::ret_type>();
   ::testing::StaticAssertTypeEq<TypeList<Dummy, int *, std::vector<Dummy>>,
                                 CallableTraits<decltype(freeFun2)>::arg_types>();
   ::testing::StaticAssertTypeEq<TypeList<Dummy &, int *, const std::vector<Dummy> &>,
                                 CallableTraits<decltype(freeFun2)>::arg_types_nodecay>();
   ::testing::StaticAssertTypeEq<Dummy, CallableTraits<decltype(freeFun2)>::ret_type>();

   // function pointer
   auto freeFun1Ptr = &freeFun1;
   auto freeFun2Ptr = &freeFun2;
   ::testing::StaticAssertTypeEq<TypeList<int>, CallableTraits<decltype(freeFun1Ptr)>::arg_types>();
   ::testing::StaticAssertTypeEq<TypeList<int>, CallableTraits<decltype(freeFun1Ptr)>::arg_types_nodecay>();
   ::testing::StaticAssertTypeEq<int, CallableTraits<decltype(freeFun1Ptr)>::ret_type>();
   ::testing::StaticAssertTypeEq<TypeList<Dummy, int *, std::vector<Dummy>>,
                                 CallableTraits<decltype(freeFun2Ptr)>::arg_types>();
   ::testing::StaticAssertTypeEq<TypeList<Dummy &, int *, const std::vector<Dummy> &>,
                                 CallableTraits<decltype(freeFun2Ptr)>::arg_types_nodecay>();
   ::testing::StaticAssertTypeEq<Dummy, CallableTraits<decltype(freeFun2Ptr)>::ret_type>();

   // functor class
   ::testing::StaticAssertTypeEq<TypeList<int>, CallableTraits<Functor1>::arg_types>();
   ::testing::StaticAssertTypeEq<TypeList<int>, CallableTraits<Functor1>::arg_types_nodecay>();
   ::testing::StaticAssertTypeEq<void, CallableTraits<Functor1>::ret_type>();
   ::testing::StaticAssertTypeEq<TypeList<Dummy, int *, std::vector<Dummy>>, CallableTraits<Functor2>::arg_types>();
   ::testing::StaticAssertTypeEq<TypeList<Dummy &, int *, const std::vector<Dummy> &>,
                                 CallableTraits<Functor2>::arg_types_nodecay>();
   ::testing::StaticAssertTypeEq<Dummy, CallableTraits<Functor2>::ret_type>();

   // mutable lambda
   auto boolvec = std::make_shared<std::vector<bool>>(2, false);
   auto lambda1 = [boolvec](const bool b) mutable -> std::vector<bool> & {
      boolvec->at(0) = b;
      return *boolvec;
   };

   // lambda
   auto lambda2 = [](Dummy &, int *, const std::vector<Dummy> &) { return Dummy(); };
   ::testing::StaticAssertTypeEq<TypeList<bool>, CallableTraits<decltype(lambda1)>::arg_types>();
   ::testing::StaticAssertTypeEq<TypeList<bool>, CallableTraits<decltype(lambda1)>::arg_types_nodecay>();
   ::testing::StaticAssertTypeEq<std::vector<bool> &, CallableTraits<decltype(lambda1)>::ret_type>();
   ::testing::StaticAssertTypeEq<TypeList<Dummy, int *, std::vector<Dummy>>,
                                 CallableTraits<decltype(lambda2)>::arg_types>();
   ::testing::StaticAssertTypeEq<TypeList<Dummy &, int *, const std::vector<Dummy> &>,
                                 CallableTraits<decltype(lambda2)>::arg_types_nodecay>();
   ::testing::StaticAssertTypeEq<Dummy, CallableTraits<decltype(lambda2)>::ret_type>();

   // std::function
   // masking signature int(int) of freeFunc1
   std::function<int(double)> stdFun1(freeFun1);
   // masking signature of lambda2
   std::function<Dummy(Dummy &, int *, const std::vector<Dummy> &)> stdFun2(lambda2);
   ::testing::StaticAssertTypeEq<TypeList<double>, CallableTraits<decltype(stdFun1)>::arg_types>();
   ::testing::StaticAssertTypeEq<TypeList<double>, CallableTraits<decltype(stdFun1)>::arg_types_nodecay>();
   ::testing::StaticAssertTypeEq<int, CallableTraits<decltype(stdFun1)>::ret_type>();
   ::testing::StaticAssertTypeEq<TypeList<Dummy, int *, std::vector<Dummy>>,
                                 CallableTraits<decltype(stdFun2)>::arg_types>();
   ::testing::StaticAssertTypeEq<TypeList<Dummy &, int *, const std::vector<Dummy> &>,
                                 CallableTraits<decltype(stdFun2)>::arg_types_nodecay>();
   ::testing::StaticAssertTypeEq<Dummy, CallableTraits<decltype(stdFun2)>::ret_type>();
}

template<typename F, typename R = typename CallableTraits<F>::ret_type>
constexpr bool HasRetType(F) { return true; }
template<typename F, typename...Args>
constexpr bool HasRetType(F, Args...) { return false; }

TEST(TypeTraits, SFINAEOnCallableTraits)
{
   EXPECT_FALSE(HasRetType(int(42)));
   EXPECT_FALSE(HasRetType(std::vector<int>()));
   EXPECT_TRUE(HasRetType(freeFun1));
   EXPECT_TRUE(HasRetType(freeFun2));
   EXPECT_TRUE(HasRetType(Functor1()));
   EXPECT_TRUE(HasRetType(Functor2()));
   EXPECT_TRUE(HasRetType(std::function<void(int)>(Functor1())));
   EXPECT_TRUE(HasRetType([]() { return 42; }));
   EXPECT_TRUE(HasRetType([]() { return Dummy(); }));
}
