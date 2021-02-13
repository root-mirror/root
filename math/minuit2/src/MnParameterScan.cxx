// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#include "Minuit2/MnParameterScan.h"
#include "Minuit2/FCNBase.h"

namespace ROOT {

namespace Minuit2 {

MnParameterScan::MnParameterScan(const FCNBase &fcn, const MnUserParameters &par)
   : fFCN(fcn), fParameters(par), fAmin(fcn(par.Params()))
{
}

MnParameterScan::MnParameterScan(const FCNBase &fcn, const MnUserParameters &par, double fval)
   : fFCN(fcn), fParameters(par), fAmin(fval)
{
}

std::vector<std::pair<double, double>> MnParameterScan::
operator()(unsigned int par, unsigned int maxsteps, double low, double high)
{
   // do the scan for parameter par between low and high values

   // if(maxsteps > 101) maxsteps = 101;
   std::vector<std::pair<double, double>> result;
   result.reserve(maxsteps + 1);
   std::vector<double> params = fParameters.Params();
   result.push_back(std::pair<double, double>(params[par], fAmin));

   if (low > high)
      return result;
   if (maxsteps < 2)
      return result;

   if (low == 0. && high == 0.) {
      low = params[par] - 2. * fParameters.Error(par);
      high = params[par] + 2. * fParameters.Error(par);
   }

   if (low == 0. && high == 0. && fParameters.Parameter(par).HasLimits()) {
      if (fParameters.Parameter(par).HasLowerLimit())
         low = fParameters.Parameter(par).LowerLimit();
      if (fParameters.Parameter(par).HasUpperLimit())
         high = fParameters.Parameter(par).UpperLimit();
   }

   if (fParameters.Parameter(par).HasLimits()) {
      if (fParameters.Parameter(par).HasLowerLimit())
         low = std::max(low, fParameters.Parameter(par).LowerLimit());
      if (fParameters.Parameter(par).HasUpperLimit())
         high = std::min(high, fParameters.Parameter(par).UpperLimit());
   }

   double x0 = low;
   double stp = (high - low) / double(maxsteps - 1);
   for (unsigned int i = 0; i < maxsteps; i++) {
      params[par] = x0 + double(i) * stp;
      double fval = fFCN(params);
      if (fval < fAmin) {
         fParameters.SetValue(par, params[par]);
         fAmin = fval;
      }
      result.push_back(std::pair<double, double>(params[par], fval));
   }

   return result;
}

} // namespace Minuit2

} // namespace ROOT
