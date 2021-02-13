// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#include "Minuit2/MnMatrix.h"

#include <cmath>

namespace ROOT {

namespace Minuit2 {

/** Inverts a symmetric matrix. Matrix is first scaled to have all ones on
    the diagonal (equivalent to change of units) but no pivoting is done
    since matrix is positive-definite.
 */

int mnvert(MnAlgebraicSymMatrix &a)
{

   unsigned int nrow = a.Nrow();
   MnAlgebraicVector s(nrow);
   MnAlgebraicVector q(nrow);
   MnAlgebraicVector pp(nrow);

   for (unsigned int i = 0; i < nrow; i++) {
      double si = a(i, i);
      if (si < 0.)
         return 1;
      s(i) = 1. / std::sqrt(si);
   }

   for (unsigned int i = 0; i < nrow; i++)
      for (unsigned int j = i; j < nrow; j++)
         a(i, j) *= (s(i) * s(j));

   for (unsigned i = 0; i < nrow; i++) {
      unsigned int k = i;
      if (a(k, k) == 0.)
         return 1;
      q(k) = 1. / a(k, k);
      pp(k) = 1.;
      a(k, k) = 0.;
      unsigned int kp1 = k + 1;
      if (k != 0) {
         for (unsigned int j = 0; j < k; j++) {
            pp(j) = a(j, k);
            q(j) = a(j, k) * q(k);
            a(j, k) = 0.;
         }
      }
      if (k != nrow - 1) {
         for (unsigned int j = kp1; j < nrow; j++) {
            pp(j) = a(k, j);
            q(j) = -a(k, j) * q(k);
            a(k, j) = 0.;
         }
      }
      for (unsigned int j = 0; j < nrow; j++)
         for (k = j; k < nrow; k++)
            a(j, k) += (pp(j) * q(k));
   }

   for (unsigned int j = 0; j < nrow; j++)
      for (unsigned int k = j; k < nrow; k++)
         a(j, k) *= (s(j) * s(k));

   return 0;
}

} // namespace Minuit2

} // namespace ROOT
