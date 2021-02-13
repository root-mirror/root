// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#include "Minuit2/LaInverse.h"
#include "Minuit2/LASymMatrix.h"

namespace ROOT {

namespace Minuit2 {

int mnvert(LASymMatrix &t);

// symmetric matrix (positive definite only)

int Invert(LASymMatrix &t)
{
   // function for inversion of symmetric matrices using  mnvert function
   // (from Fortran Minuit)

   int ifail = 0;

   if (t.size() == 1) {
      double tmp = t.Data()[0];
      if (!(tmp > 0.))
         ifail = 1;
      else
         t.Data()[0] = 1. / tmp;
   } else {
      ifail = mnvert(t);
   }

   return ifail;
}

} // namespace Minuit2

} // namespace ROOT
