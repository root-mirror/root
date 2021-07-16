// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei, E.G.P. Bos   2003-2017

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#include <vector>
#include "Minuit2/ExternalInternalGradientCalculator.h"
#include "Minuit2/FCNGradientBase.h"
#include "Minuit2/MnUserTransformation.h"
#include "Minuit2/FunctionGradient.h"
#include "Minuit2/MinimumParameters.h"
#include "Minuit2/MnPrint.h"

namespace ROOT {
namespace Minuit2 {

FunctionGradient ExternalInternalGradientCalculator::operator()(const MinimumParameters &par) const
{
   std::vector<double> par_vec;
   par_vec.resize(par.Vec().size());
   for (std::size_t ix = 0; ix < par.Vec().size(); ++ix) {
      par_vec[ix] = par.Vec()(ix);
   }

   std::vector<double> grad = fGradCalc.Gradient(par_vec);
   assert(grad.size() == fTransformation.Parameters().size());

   MnAlgebraicVector v(par.Vec().size());
   for (unsigned int i = 0; i < par.Vec().size(); i++) {
      unsigned int ext = fTransformation.ExtOfInt(i);
      v(i) = grad[ext];
   }

   MnPrint print("ExternalInternalGradientCalculator");
   print.Debug("User given gradient in Minuit2", v);

   return FunctionGradient(v);
}

FunctionGradient
ExternalInternalGradientCalculator::operator()(const MinimumParameters &par, const FunctionGradient &functionGradient) const
{
   std::vector<double> par_vec;
   par_vec.resize(par.Vec().size());
   for (std::size_t ix = 0; ix < par.Vec().size(); ++ix) {
      par_vec[ix] = par.Vec()(ix);
   }

   std::vector<double> previous_grad(functionGradient.Grad().Data(), functionGradient.Grad().Data() + functionGradient.Grad().size());
   std::vector<double> previous_g2(functionGradient.G2().Data(), functionGradient.G2().Data() + functionGradient.G2().size());
   std::vector<double> previous_gstep(functionGradient.Gstep().Data(), functionGradient.Gstep().Data() + functionGradient.Gstep().size());

   std::vector<double> grad = fGradCalc.Gradient(par_vec, previous_grad.data(), previous_g2.data(), previous_gstep.data());
   assert(grad.size() == fTransformation.Parameters().size());

   MnAlgebraicVector v(par.Vec().size());
   MnAlgebraicVector v_g2(par.Vec().size());
   MnAlgebraicVector v_gstep(par.Vec().size());
   for (unsigned int i = 0; i < par.Vec().size(); i++) {
      unsigned int ext = fTransformation.ExtOfInt(i);
      v(i) = grad[ext];
      v_g2(i) = previous_g2[ext];
      v_gstep(i) = previous_gstep[ext];
   }

   MnPrint print("ExternalInternalGradientCalculator");
   print.Debug("User given gradient in Minuit2", v, "g2", v_g2, "step size", v_gstep);

   return FunctionGradient(v, v_g2, v_gstep);
}

} // namespace Minuit2
} // namespace ROOT
