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

double mndasum(unsigned int, const double *, int);

double sum_of_elements(const LAVector &v)
{
   // calculate the absolute sum of the vector elements using mndasum
   // which is a translation from dasum from BLAS
   return mndasum(v.size(), v.Data(), 1);
}

double sum_of_elements(const LASymMatrix &m)
{
   // calculate the absolute sum of all the matrix elements using mndasum
   // which is a translation of dasum from BLAS
   return mndasum(m.size(), m.Data(), 1);
}

} // namespace Minuit2

} // namespace ROOT
