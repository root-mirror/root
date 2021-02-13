// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Minuit2_LaSum
#define ROOT_Minuit2_LaSum

#include "Minuit2/ABSum.h"
#include "Minuit2/LAVector.h"
#include "Minuit2/LASymMatrix.h"

namespace ROOT {

namespace Minuit2 {

#define OP_ADD1(MT, MAT1, T)                                                                                 \
   inline ABObj<MT, MAT1, T> operator-(const ABObj<MT, MAT1, T> &m)                                          \
   {                                                                                                         \
      return ABObj<MT, MAT1, T>(m.Obj(), T(-1.) * m.f());                                                    \
   }                                                                                                         \
                                                                                                             \
   inline ABObj<MT, ABSum<ABObj<MT, MAT1, T>, ABObj<MT, MAT1, T>>, T> operator+(const ABObj<MT, MAT1, T> &a, \
                                                                                const ABObj<MT, MAT1, T> &b) \
   {                                                                                                         \
      return ABObj<MT, ABSum<ABObj<MT, MAT1, T>, ABObj<MT, MAT1, T>>, T>(                                    \
         ABSum<ABObj<MT, MAT1, T>, ABObj<MT, MAT1, T>>(a, b));                                               \
   }                                                                                                         \
   inline ABObj<MT, ABSum<ABObj<MT, MAT1, T>, ABObj<MT, MAT1, T>>, T> operator-(const ABObj<MT, MAT1, T> &a, \
                                                                                const ABObj<MT, MAT1, T> &b) \
   {                                                                                                         \
      return ABObj<MT, ABSum<ABObj<MT, MAT1, T>, ABObj<MT, MAT1, T>>, T>(                                    \
         ABSum<ABObj<MT, MAT1, T>, ABObj<MT, MAT1, T>>(a, ABObj<MT, MAT1, T>(b.Obj(), T(-1.) * b.f())));     \
   }

OP_ADD1(vec, LAVector, double)
OP_ADD1(sym, LASymMatrix, double)

#define OP_SCALE(MT, MAT1, T) \
   inline ABObj<MT, MAT1, T> operator*(T f, const MAT1 &obj) { return ABObj<MT, MAT1, T>(obj, f); }

OP_SCALE(sym, LASymMatrix, double)
OP_SCALE(vec, LAVector, double)

#define OP_SCALE1(MT, MAT1, T) \
   inline ABObj<MT, MAT1, T> operator/(const MAT1 &obj, T f) { return ABObj<MT, MAT1, T>(obj, 1. / f); }

OP_SCALE1(sym, LASymMatrix, double)
OP_SCALE1(vec, LAVector, double)

#define OP_MIN(MT, MAT1, T) \
   inline ABObj<MT, MAT1, T> operator-(const MAT1 &obj) { return ABObj<MT, MAT1, T>(obj, T(-1.)); }

OP_MIN(sym, LASymMatrix, double)
OP_MIN(vec, LAVector, double)

} // namespace Minuit2

} // namespace ROOT

#endif // ROOT_Minuit2_LaSum
