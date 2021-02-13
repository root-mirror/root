// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#include "Minuit2/MnScan.h"
#include "Minuit2/MnParameterScan.h"

namespace ROOT {

namespace Minuit2 {

std::vector<std::pair<double, double>> MnScan::Scan(unsigned int par, unsigned int maxsteps, double low, double high)
{
   // perform a scan of the function in the parameter par
   MnParameterScan scan(fFCN, fState.Parameters());
   double amin = scan.Fval();

   std::vector<std::pair<double, double>> result = scan(par, maxsteps, low, high);
   if (scan.Fval() < amin) {
      fState.SetValue(par, scan.Parameters().Value(par));
      amin = scan.Fval();
   }

   return result;
}

} // namespace Minuit2

} // namespace ROOT
