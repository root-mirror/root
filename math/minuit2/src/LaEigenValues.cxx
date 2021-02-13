// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#include "Minuit2/LAVector.h"
#include "Minuit2/LASymMatrix.h"

namespace ROOT {

namespace Minuit2 {

int mneigen(double *, unsigned int, unsigned int, unsigned int, double *, double);

LAVector eigenvalues(const LASymMatrix &mat)
{
   // calculate eigenvalues of symmetric matrices using mneigen function (transalte from fortran Minuit)
   unsigned int nrow = mat.Nrow();

   LAVector tmp(nrow * nrow);
   LAVector work(2 * nrow);

   for (unsigned int i = 0; i < nrow; i++)
      for (unsigned int j = 0; j <= i; j++) {
         tmp(i + j * nrow) = mat(i, j);
         tmp(i * nrow + j) = mat(i, j);
      }

   int info = mneigen(tmp.Data(), nrow, nrow, work.size(), work.Data(), 1.e-6);
   (void)info;
   assert(info == 0);

   LAVector result(nrow);
   for (unsigned int i = 0; i < nrow; i++)
      result(i) = work(i);

   return result;
}

} // namespace Minuit2

} // namespace ROOT
