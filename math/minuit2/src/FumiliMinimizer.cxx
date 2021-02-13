// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#include "Minuit2/MnConfig.h"
#include "Minuit2/FumiliMinimizer.h"
#include "Minuit2/MinimumSeedGenerator.h"
#include "Minuit2/FumiliGradientCalculator.h"
#include "Minuit2/Numerical2PGradientCalculator.h"
#include "Minuit2/AnalyticalGradientCalculator.h"
#include "Minuit2/MinimumBuilder.h"
#include "Minuit2/MinimumSeed.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnUserTransformation.h"
#include "Minuit2/MnUserFcn.h"
#include "Minuit2/FumiliFCNBase.h"
#include "Minuit2/FCNGradientBase.h"
#include "Minuit2/MnStrategy.h"
#include "Minuit2/MnPrint.h"

namespace ROOT {

namespace Minuit2 {

// for Fumili implement Minimize here because need downcast

FunctionMinimum FumiliMinimizer::Minimize(const FCNBase &fcn, const MnUserParameterState &st,
                                          const MnStrategy &strategy, unsigned int maxfcn, double toler) const
{
   // Minimize using Fumili. Create seed and Fumili gradient calculator.
   // The FCNBase passed must be a FumiliFCNBase type otherwise method will fail !

   MnPrint print("FumiliMinimizer");

   MnUserFcn mfcn(fcn, st.Trafo());
   Numerical2PGradientCalculator gc(mfcn, st.Trafo(), strategy);

   unsigned int npar = st.VariableParameters();
   if (maxfcn == 0)
      maxfcn = 200 + 100 * npar + 5 * npar * npar;
   // FUMILI needs much less function calls
   maxfcn = int(0.1 * maxfcn);

   MinimumSeed mnseeds = SeedGenerator()(mfcn, gc, st, strategy);

   // downcast fcn

   // std::cout << "FCN type " << typeid(&fcn).Name() << std::endl;

   FumiliFCNBase *fumiliFcn = dynamic_cast<FumiliFCNBase *>(const_cast<FCNBase *>(&fcn));
   if (!fumiliFcn) {
      print.Error("Wrong FCN type; try to use default minimizer");
      return FunctionMinimum(mnseeds, fcn.Up());
   }

   FumiliGradientCalculator fgc(*fumiliFcn, st.Trafo(), npar);
   print.Debug("Using FumiliMinimizer");

   return ModularFunctionMinimizer::Minimize(mfcn, fgc, mnseeds, strategy, maxfcn, toler);
}

FunctionMinimum FumiliMinimizer::Minimize(const FCNGradientBase &fcn, const MnUserParameterState &st,
                                          const MnStrategy &strategy, unsigned int maxfcn, double toler) const
{

   MnPrint print("FumiliMinimizer::Minimize");

   // Minimize using Fumili. Case of interface is a FCNGradientBase.
   // Normally other method is used  - probably this could be removed (t.b.i.)

   // need MnUserFcn
   MnUserFcn mfcn(fcn, st.Trafo());
   AnalyticalGradientCalculator gc(fcn, st.Trafo());

   unsigned int npar = st.VariableParameters();
   if (maxfcn == 0)
      maxfcn = 200 + 100 * npar + 5 * npar * npar;

   MinimumSeed mnseeds = SeedGenerator()(mfcn, gc, st, strategy);

   // downcast fcn

   FumiliFCNBase *fumiliFcn = dynamic_cast<FumiliFCNBase *>(const_cast<FCNGradientBase *>(&fcn));
   if (!fumiliFcn) {
      print.Error("Wrong FCN type; try to use default minimizer");
      return FunctionMinimum(mnseeds, fcn.Up());
   }

   FumiliGradientCalculator fgc(*fumiliFcn, st.Trafo(), npar);
   print.Debug("Using FumiliMinimizer");

   return ModularFunctionMinimizer::Minimize(mfcn, fgc, mnseeds, strategy, maxfcn, toler);
}

} // namespace Minuit2

} // namespace ROOT
