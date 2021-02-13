// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef MA_LaOuterProd_H_
#define MA_LaOuterProd_H_

#include "Minuit2/VectorOuterProduct.h"
#include "Minuit2/ABSum.h"
#include "Minuit2/LAVector.h"
#include "Minuit2/LASymMatrix.h"

namespace ROOT {

namespace Minuit2 {

///    LAPACK Algebra function
///    specialize the Outer_product function for LAVector;

inline ABObj<sym, VectorOuterProduct<ABObj<vec, LAVector, double>, double>, double>
Outer_product(const ABObj<vec, LAVector, double> &obj)
{
   //   std::cout<<"ABObj<sym, VectorOuterProduct<ABObj<vec, LAVector, double>, double>, double> Outer_product(const
   //   ABObj<vec, LAVector, double>& obj)"<<std::endl;
   return ABObj<sym, VectorOuterProduct<ABObj<vec, LAVector, double>, double>, double>(
      VectorOuterProduct<ABObj<vec, LAVector, double>, double>(obj));
}

// f*outer
template <class T>
inline ABObj<sym, VectorOuterProduct<ABObj<vec, LAVector, T>, T>, T>
operator*(T f, const ABObj<sym, VectorOuterProduct<ABObj<vec, LAVector, T>, T>, T> &obj)
{
   //   std::cout<<"ABObj<sym, VectorOuterProduct<ABObj<vec, LAVector, T>, T>, T> operator*(T f, const ABObj<sym,
   //   VectorOuterProduct<ABObj<vec, LAVector, T>, T>, T>& obj)"<<std::endl;
   return ABObj<sym, VectorOuterProduct<ABObj<vec, LAVector, T>, T>, T>(obj.Obj(), obj.f() * f);
}

// outer/f
template <class T>
inline ABObj<sym, VectorOuterProduct<ABObj<vec, LAVector, T>, T>, T>
operator/(const ABObj<sym, VectorOuterProduct<ABObj<vec, LAVector, T>, T>, T> &obj, T f)
{
   //   std::cout<<"ABObj<sym, VectorOuterProduct<ABObj<vec, LAVector, T>, T>, T> operator/(const ABObj<sym,
   //   VectorOuterProduct<ABObj<vec, LAVector, T>, T>, T>& obj, T f)"<<std::endl;
   return ABObj<sym, VectorOuterProduct<ABObj<vec, LAVector, T>, T>, T>(obj.Obj(), obj.f() / f);
}

// -outer
template <class T>
inline ABObj<sym, VectorOuterProduct<ABObj<vec, LAVector, T>, T>, T>
operator-(const ABObj<sym, VectorOuterProduct<ABObj<vec, LAVector, T>, T>, T> &obj)
{
   //   std::cout<<"ABObj<sym, VectorOuterProduct<ABObj<vec, LAVector, T>, T>, T> operator/(const ABObj<sym,
   //   VectorOuterProduct<ABObj<vec, LAVector, T>, T>, T>& obj, T f)"<<std::endl;
   return ABObj<sym, VectorOuterProduct<ABObj<vec, LAVector, T>, T>, T>(obj.Obj(), T(-1.) * obj.f());
}

void Outer_prod(LASymMatrix &, const LAVector &, double f = 1.);

} // namespace Minuit2

} // namespace ROOT

#endif // MA_LaOuterProd_H_
