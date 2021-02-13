#ifndef _DICT_CHECK_TEST_
#define _DICT_CHECK_TEST_

#include "stlDictCheckAux.h"

#include <array>
#include <memory>
#include <utility>

#ifdef __ROOTCLING__
#pragma link C++ class A+;
#pragma link C++ class B+;
#endif


struct A {
};

struct D {
};

struct B {
   std::unique_ptr<A> fA;
   std::tuple<A, int, double*> fT;
   std::array<int, 5> fArr;
   std::array<const std::unique_ptr<C**const*>, 5> fArrPtr;
   std::array<const std::unique_ptr<D**const*>, 5> fArrPtr2;
};

#endif

