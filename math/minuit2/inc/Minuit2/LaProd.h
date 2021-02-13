// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Minuit2_LaProd
#define ROOT_Minuit2_LaProd

#include "Minuit2/ABProd.h"
#include "Minuit2/LAVector.h"
#include "Minuit2/LASymMatrix.h"

namespace ROOT {

namespace Minuit2 {

/*
  LM" remove this for fixing alpha poblem
#define OP_MULT1(MT1,MT2,MAT1,MAT2,T) \
template<class B>                     \
inline ABObj<AlgebraicProdType<MT1,MT2>::Type,ABProd<ABObj<MT1,MAT1,T>, ABObj<MT2,B,T> >,T>  operator*(const
ABObj<MT1,MAT1,T>& a, const ABObj<MT2,B,T>& b) { return ABObj<AlgebraicProdType<MT1,MT2>::Type,ABProd<ABObj<MT1,MAT1,T>,
ABObj<MT2,B,T> >,T>(ABProd<ABObj<MT1,MAT1,T>, ABObj<MT2,B,T> >(a, b)); }   \
template<class A>             \
inline ABObj<AlgebraicProdType<MT1,MT2>::Type,ABProd<ABObj<MT1,A,T>, ABObj<MT2,MAT2,T> >,T>  operator*(const
ABObj<MT1,A,T>& a, const ABObj<MT2,MAT2,T>& b) { \
  return ABObj<AlgebraicProdType<MT1,MT2>::Type,ABProd<ABObj<MT1,A,T>, ABObj<MT2,MAT2,T> >,T>(ABProd<ABObj<MT1,A,T>,
ABObj<MT2,MAT2,T> >(a, b));    \
} \
  \
*/

#define OP_MULT1(MT1, MT2, MAT1, MAT2, T)                                                                          \
   inline ABObj<AlgebraicProdType<MT1, MT2>::Type, ABProd<ABObj<MT1, MAT1, T>, ABObj<MT2, MAT2, T>>, T> operator*( \
      const ABObj<MT1, MAT1, T> &a, const ABObj<MT2, MAT2, T> &b)                                                  \
   {                                                                                                               \
      return ABObj<AlgebraicProdType<MT1, MT2>::Type, ABProd<ABObj<MT1, MAT1, T>, ABObj<MT2, MAT2, T>>, T>(        \
         ABProd<ABObj<MT1, MAT1, T>, ABObj<MT2, MAT2, T>>(a, b));                                                  \
   }

OP_MULT1(sym, vec, LASymMatrix, LAVector, double)
// OP_MULT1(sym,gen,LASymMatrix,LAGenMatrix,double)
// OP_MULT1(sym,sym,LASymMatrix,LASymMatrix,double)
// OP_MULT1(gen,vec,LAGenMatrix,LAVector,double)
// OP_MULT1(gen,sym,LAGenMatrix,LASymMatrix,double)
// OP_MULT1(gen,gen,LAGenMatrix,LAGenMatrix,double)

} // namespace Minuit2

} // namespace ROOT

#endif // ROOT_Minuit2_LaProd
