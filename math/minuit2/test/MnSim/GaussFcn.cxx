// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#include "GaussFcn.h"
#include "GaussFunction.h"

#include <cassert>

namespace ROOT {

namespace Minuit2 {

double GaussFcn::operator()(const std::vector<double> &par) const
{

   assert(par.size() == 3);
   GaussFunction gauss(par[0], par[1], par[2]);

   double chi2 = 0.;
   for (unsigned int n = 0; n < fMeasurements.size(); n++) {
      chi2 += ((gauss(fPositions[n]) - fMeasurements[n]) * (gauss(fPositions[n]) - fMeasurements[n]) / fMVariances[n]);
   }

   return chi2;
}

} // namespace Minuit2

} // namespace ROOT
