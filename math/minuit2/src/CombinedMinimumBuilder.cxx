// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#include "Minuit2/CombinedMinimumBuilder.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnStrategy.h"
#include "Minuit2/MnPrint.h"

namespace ROOT {

namespace Minuit2 {

FunctionMinimum CombinedMinimumBuilder::Minimum(const MnFcn &fcn, const GradientCalculator &gc, const MinimumSeed &seed,
                                                const MnStrategy &strategy, unsigned int maxfcn, double edmval) const
{
   // find minimum using combined method
   // (Migrad then if fails try Simplex and then Migrad again)

   MnPrint print("CombinedMinimumBuilder");

   FunctionMinimum min = fVMMinimizer.Builder().Minimum(fcn, gc, seed, strategy, maxfcn, edmval);

   if (!min.IsValid()) {
      print.Warn("Migrad method fails, will try with simplex method first");

      MnStrategy str(2);
      FunctionMinimum min1 = fSimplexMinimizer.Builder().Minimum(fcn, gc, seed, str, maxfcn, edmval);
      if (!min1.IsValid()) {
         print.Warn("Both Migrad and Simplex methods failed");

         return min1;
      }
      MinimumSeed seed1 = fVMMinimizer.SeedGenerator()(fcn, gc, min1.UserState(), str);

      FunctionMinimum min2 = fVMMinimizer.Builder().Minimum(fcn, gc, seed1, str, maxfcn, edmval);
      if (!min2.IsValid()) {

         print.Warn("Both migrad and method failed also at 2nd attempt; return simplex Minimum");
         return min1;
      }

      return min2;
   }

   return min;
}

} // namespace Minuit2

} // namespace ROOT
