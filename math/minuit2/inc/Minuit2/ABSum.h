// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Minuit2_ABSum
#define ROOT_Minuit2_ABSum

#include "Minuit2/ABObj.h"

namespace ROOT {

namespace Minuit2 {

template <class M1, class M2>
class ABSum {

private:
   ABSum() : fA(M1()), fB(M2()) {}

   ABSum &operator=(const ABSum &) { return *this; }

   template <class MI1, class MI2>
   ABSum &operator=(const ABSum<MI1, MI2> &)
   {
      return *this;
   }

public:
   ABSum(const M1 &a, const M2 &b) : fA(a), fB(b) {}

   ~ABSum() {}

   ABSum(const ABSum &sum) : fA(sum.fA), fB(sum.fB) {}

   template <class MI1, class MI2>
   ABSum(const ABSum<MI1, MI2> &sum) : fA(M1(sum.A())), fB(M2(sum.B()))
   {
   }

   const M1 &A() const { return fA; }
   const M2 &B() const { return fB; }

private:
   M1 fA;
   M2 fB;
};

// ABObj + ABObj
template <class atype, class A, class btype, class B, class T>
inline ABObj<typename AlgebraicSumType<atype, btype>::Type, ABSum<ABObj<atype, A, T>, ABObj<btype, B, T>>, T>
operator+(const ABObj<atype, A, T> &a, const ABObj<btype, B, T> &b)
{

   return ABObj<typename AlgebraicSumType<atype, btype>::Type, ABSum<ABObj<atype, A, T>, ABObj<btype, B, T>>, T>(
      ABSum<ABObj<atype, A, T>, ABObj<btype, B, T>>(a, b));
}

// ABObj - ABObj
template <class atype, class A, class btype, class B, class T>
inline ABObj<typename AlgebraicSumType<atype, btype>::Type, ABSum<ABObj<atype, A, T>, ABObj<btype, B, T>>, T>
operator-(const ABObj<atype, A, T> &a, const ABObj<btype, B, T> &b)
{

   return ABObj<typename AlgebraicSumType<atype, btype>::Type, ABSum<ABObj<atype, A, T>, ABObj<btype, B, T>>, T>(
      ABSum<ABObj<atype, A, T>, ABObj<btype, B, T>>(a, ABObj<btype, B, T>(b.Obj(), T(-1.) * b.f())));
}

} // namespace Minuit2

} // namespace ROOT

#endif // ROOT_Minuit2_ABSum
