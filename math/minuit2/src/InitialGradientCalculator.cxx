// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#include "Minuit2/InitialGradientCalculator.h"
#include "Minuit2/MnFcn.h"
#include "Minuit2/MnUserTransformation.h"
#include "Minuit2/MnMachinePrecision.h"
#include "Minuit2/MinimumParameters.h"
#include "Minuit2/FunctionGradient.h"
#include "Minuit2/MnStrategy.h"
#include "Minuit2/MnPrint.h"

#include <cmath>

namespace ROOT {

namespace Minuit2 {

FunctionGradient InitialGradientCalculator::operator()(const MinimumParameters &par) const
{
   // initial rough  estimate of the gradient using the parameter step size

//   std::cout << "########### InitialGradientCalculator::operator()" <<std::endl;

   assert(par.IsValid());

   unsigned int n = Trafo().VariableParameters();
   assert(n == par.Vec().size());

   MnPrint print("InitialGradientCalculator");

   print.Debug("Calculating initial gradient at point", par.Vec());

   MnAlgebraicVector gr(n), gr2(n), gst(n);

   for (unsigned int i = 0; i < n; i++) {
      unsigned int exOfIn = Trafo().ExtOfInt(i);

      double var = par.Vec()(i);

      double werr = Trafo().Parameter(exOfIn).Error();
      double sav = Trafo().Int2ext(i, var);

      double sav2 = sav + werr;
      if (Trafo().Parameter(exOfIn).HasLimits()) {
         if (Trafo().Parameter(exOfIn).HasUpperLimit() && sav2 > Trafo().Parameter(exOfIn).UpperLimit())
            sav2 = Trafo().Parameter(exOfIn).UpperLimit();
      }
      double var2 = Trafo().Ext2int(exOfIn, sav2);
      double vplu = var2 - var;
      sav2 = sav - werr;
      if (Trafo().Parameter(exOfIn).HasLimits()) {
         if (Trafo().Parameter(exOfIn).HasLowerLimit() && sav2 < Trafo().Parameter(exOfIn).LowerLimit())
            sav2 = Trafo().Parameter(exOfIn).LowerLimit();
      }
      var2 = Trafo().Ext2int(exOfIn, sav2);
      double vmin = var2 - var;
      double gsmin = 8. * Precision().Eps2() * (std::fabs(var) + Precision().Eps2());
      // protect against very small step sizes which can cause dirin to zero and then nan values in grd
      double dirin = std::max(0.5 * (std::fabs(vplu) + std::fabs(vmin)), gsmin);
      double g2 = 2.0 * fFcn.ErrorDef() / (dirin * dirin);
      double gstep = std::max(gsmin, 0.1 * dirin);
      double grd = g2 * dirin;
      if (Trafo().Parameter(exOfIn).HasLimits()) {
         if (gstep > 0.5)
            gstep = 0.5;
      }
      gr(i) = grd;
      gr2(i) = g2;
      gst(i) = gstep;

      print.Debug("Computed initial gradient for parameter", Trafo().Name(exOfIn), "value", var, "[", vmin, ",", vplu,
                  "]", "dirin", dirin, "grd", grd, "g2", g2);
   }

   return FunctionGradient(gr, gr2, gst);
}

FunctionGradient InitialGradientCalculator::operator()(const MinimumParameters &par, const FunctionGradient &) const
{
   // Base class interface
   return (*this)(par);
}

const MnMachinePrecision &InitialGradientCalculator::Precision() const
{
   // return precision (is set in trasformation class)
   return fTransformation.Precision();
}

unsigned int InitialGradientCalculator::Ncycle() const
{
   // return ncyles (from Strategy)
   return Strategy().GradientNCycles();
}

double InitialGradientCalculator::StepTolerance() const
{
   // return Gradient step tolerance (from Strategy)
   return Strategy().GradientStepTolerance();
}

double InitialGradientCalculator::GradTolerance() const
{
   // return Gradient tolerance
   return Strategy().GradientTolerance();
}

} // namespace Minuit2

} // namespace ROOT
