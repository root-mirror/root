// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Minuit2_VectorOuterProduct
#define ROOT_Minuit2_VectorOuterProduct

#include "Minuit2/ABTypes.h"
#include "Minuit2/ABObj.h"

namespace ROOT {

namespace Minuit2 {

template <class M, class T>
class VectorOuterProduct {

public:
   VectorOuterProduct(const M &obj) : fObject(obj) {}

   ~VectorOuterProduct() {}

   typedef sym Type;

   const M &Obj() const { return fObject; }

private:
   M fObject;
};

template <class M, class T>
inline ABObj<sym, VectorOuterProduct<ABObj<vec, M, T>, T>, T> Outer_product(const ABObj<vec, M, T> &obj)
{
   return ABObj<sym, VectorOuterProduct<ABObj<vec, M, T>, T>, T>(VectorOuterProduct<ABObj<vec, M, T>, T>(obj));
}

} // namespace Minuit2

} // namespace ROOT

#endif // ROOT_Minuit2_VectorOuterProduct
