// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#include "Minuit2/LAVector.h"

namespace ROOT {

namespace Minuit2 {

double mnddot(unsigned int, const double *, int, const double *, int);

double inner_product(const LAVector &v1, const LAVector &v2)
{
   // calculate inner (dot) product of two vectors  using mnddot function
   return mnddot(v1.size(), v1.Data(), 1, v2.Data(), 1);
}

} // namespace Minuit2

} // namespace ROOT
