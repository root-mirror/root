// @(#)root/mathcore:$Id$
// Authors: L. Moneta, J.T. Offermann, E.G.P. Bos    2013-2018
//
/**********************************************************************
 *                                                                    *
 * Copyright (c) 2013 , LCG ROOT MathLib Team                         *
 *                                                                    *
 **********************************************************************/
/**
 * \class NumericalDerivator
 *
 *  Original version created on: Aug 14, 2013
 *      Authors: L. Moneta, J. T. Offermann
 *  Modified version created on: Sep 27, 2017
 *      Author: E. G. P. Bos
 *
 *      NumericalDerivator was essentially a slightly modified copy of code
 *      written by M. Winkler, F. James, L. Moneta, and A. Zsenei for Minuit2,
 *      Copyright (c) 2005 LCG ROOT Math team, CERN/PH-SFT. Original version:
 *      https://github.com/lmoneta/root/blob/lvmini/math/mathcore/src/NumericalDerivator.cxx
 *
 *      This class attempts to more closely follow the Minuit2 implementation.
 *      Modified things (w.r.t. original) are indicated by MODIFIED.
 */

#include "Minuit2/NumericalDerivator.h"
#include <cmath>
#include <algorithm>
#include <Math/IFunction.h>
#include <iostream>
#include <TMath.h>
#include <cassert>
#include "Fit/ParameterSettings.h"

#include <Math/Minimizer.h> // needed here because in Fitter is only a forward declaration

namespace ROOT {
namespace Minuit2 {

NumericalDerivator::NumericalDerivator(bool always_exactly_mimic_minuit2)
   : fAlwaysExactlyMimicMinuit2(always_exactly_mimic_minuit2)
{
}

NumericalDerivator::NumericalDerivator(double step_tolerance, double grad_tolerance, unsigned int ncycles,
                                       double error_level, bool always_exactly_mimic_minuit2)
   : fStepTolerance(step_tolerance), fGradTolerance(grad_tolerance), fUp(error_level), fNCycles(ncycles),
     fAlwaysExactlyMimicMinuit2(always_exactly_mimic_minuit2)
{
}

NumericalDerivator::NumericalDerivator(const NumericalDerivator &/*other*/) = default;


/// This function sets internal state based on input parameters. This state
/// setup is used in the actual (partial) derivative calculations.
void NumericalDerivator::SetupDifferentiate(const ROOT::Math::IBaseFunctionMultiDim *function, const double *cx,
                                            const std::vector<ROOT::Fit::ParameterSettings> &parameters)
{
   assert(function != nullptr && "function is a nullptr");

   fVx.resize(function->NDim());
   fVxExternal.resize(function->NDim());
   fVxFValCache.resize(function->NDim());
   std::copy(cx, cx + function->NDim(), fVx.data());

   // convert to Minuit external parameters
   for (unsigned i = 0; i < function->NDim(); i++) {
      fVxExternal[i] = Int2ext(parameters[i], fVx[i]);
   }

   if (fVx != fVxFValCache) {
      fVxFValCache = fVx;
      fVal = (*function)(fVxExternal.data()); // value of function at given points
   }

   fDfmin = 8. * fPrecision.Eps2() * (std::abs(fVal) + fUp);
   fVrysml = 8. * fPrecision.Eps() * fPrecision.Eps();
}

DerivatorElement NumericalDerivator::PartialDerivative(const ROOT::Math::IBaseFunctionMultiDim *function,
                                                       const double *x,
                                                       const std::vector<ROOT::Fit::ParameterSettings> &parameters,
                                                       unsigned int i_component, DerivatorElement previous)
{
   SetupDifferentiate(function, x, parameters);
   return FastPartialDerivative(function, parameters, i_component, previous);
}

// leaves the parameter setup to the caller
DerivatorElement NumericalDerivator::FastPartialDerivative(const ROOT::Math::IBaseFunctionMultiDim *function,
                                                           const std::vector<ROOT::Fit::ParameterSettings> &parameters,
                                                           unsigned int i_component, const DerivatorElement &previous)
{
   DerivatorElement deriv{previous};

   double xtf = fVx[i_component];
   double epspri = fPrecision.Eps2() + std::abs(deriv.derivative * fPrecision.Eps2());
   double step_old = 0.;

   for (unsigned int j = 0; j < fNCycles; ++j) {
      double optstp = std::sqrt(fDfmin / (std::abs(deriv.second_derivative) + epspri));
      double step = std::max(optstp, std::abs(0.1 * deriv.step_size));

      if (parameters[i_component].IsBound()) {
         if (step > 0.5)
            step = 0.5;
      }

      double stpmax = 10. * std::abs(deriv.step_size);
      if (step > stpmax)
         step = stpmax;

      double stpmin = std::max(fVrysml, 8. * std::abs(fPrecision.Eps2() * fVx[i_component]));
      if (step < stpmin)
         step = stpmin;
      if (std::abs((step - step_old) / step) < fStepTolerance) {
         break;
      }
      deriv.step_size = step;
      step_old = step;
      fVx[i_component] = xtf + step;
      fVxExternal[i_component] = Int2ext(parameters[i_component], fVx[i_component]);
      double fs1 = (*function)(fVxExternal.data());

      fVx[i_component] = xtf - step;
      fVxExternal[i_component] = Int2ext(parameters[i_component], fVx[i_component]);
      double fs2 = (*function)(fVxExternal.data());

      fVx[i_component] = xtf;
      fVxExternal[i_component] = Int2ext(parameters[i_component], fVx[i_component]);

      double fGrd_old = deriv.derivative;
      deriv.derivative = 0.5 * (fs1 - fs2) / step;

      deriv.second_derivative = (fs1 + fs2 - 2. * fVal) / step / step;

      if (std::abs(fGrd_old - deriv.derivative) / (std::abs(deriv.derivative) + fDfmin / step) < fGradTolerance) {
         break;
      }
   }
   return deriv;
}

DerivatorElement NumericalDerivator::operator()(const ROOT::Math::IBaseFunctionMultiDim *function, const double *x,
                                                const std::vector<ROOT::Fit::ParameterSettings> &parameters,
                                                unsigned int i_component, const DerivatorElement &previous)
{
   return PartialDerivative(function, x, parameters, i_component, previous);
}

std::vector<DerivatorElement>
NumericalDerivator::Differentiate(const ROOT::Math::IBaseFunctionMultiDim *function, const double *cx,
                                  const std::vector<ROOT::Fit::ParameterSettings> &parameters,
                                  const std::vector<DerivatorElement> &previous_gradient)
{
   SetupDifferentiate(function, cx, parameters);

   std::vector<DerivatorElement> gradient;
   gradient.reserve(function->NDim());

   for (unsigned int ix = 0; ix < function->NDim(); ++ix) {
      gradient.emplace_back(FastPartialDerivative(function, parameters, ix, previous_gradient[ix]));
   }

   return gradient;
}

double NumericalDerivator::Int2ext(const ROOT::Fit::ParameterSettings &parameter, double val) const
{
   // return external value from internal value for parameter i
   if (parameter.IsBound()) {
      if (parameter.IsDoubleBound()) {
         return fDoubleLimTrafo.Int2ext(val, parameter.UpperLimit(), parameter.LowerLimit());
      } else if (parameter.HasUpperLimit() && !parameter.HasLowerLimit()) {
         return fUpperLimTrafo.Int2ext(val, parameter.UpperLimit());
      } else {
         return fLowerLimTrafo.Int2ext(val, parameter.LowerLimit());
      }
   }

   return val;
}

double NumericalDerivator::Ext2int(const ROOT::Fit::ParameterSettings &parameter, double val) const
{
   // return the internal value for parameter i with external value val

   if (parameter.IsBound()) {
      if (parameter.IsDoubleBound()) {
         return fDoubleLimTrafo.Ext2int(val, parameter.UpperLimit(), parameter.LowerLimit(), fPrecision);
      } else if (parameter.HasUpperLimit() && !parameter.HasLowerLimit()) {
         return fUpperLimTrafo.Ext2int(val, parameter.UpperLimit(), fPrecision);
      } else {
         return fLowerLimTrafo.Ext2int(val, parameter.LowerLimit(), fPrecision);
      }
   }

   return val;
}

double NumericalDerivator::DInt2Ext(const ROOT::Fit::ParameterSettings &parameter, double val) const
{
   // return the derivative of the int->ext transformation: dPext(i) / dPint(i)
   // for the parameter i with value val

   double dd = 1.;
   if (parameter.IsBound()) {
      if (parameter.IsDoubleBound()) {
         dd = fDoubleLimTrafo.DInt2Ext(val, parameter.UpperLimit(), parameter.LowerLimit());
      } else if (parameter.HasUpperLimit() && !parameter.HasLowerLimit()) {
         dd = fUpperLimTrafo.DInt2Ext(val, parameter.UpperLimit());
      } else {
         dd = fLowerLimTrafo.DInt2Ext(val, parameter.LowerLimit());
      }
   }

   return dd;
}

// MODIFIED:
/// This function was not implemented as in Minuit2. Now it copies the behavior
/// of InitialGradientCalculator. See https://github.com/roofit-dev/root/issues/10
void NumericalDerivator::SetInitialGradient(const std::vector<ROOT::Fit::ParameterSettings> &parameters,
                                            std::vector<DerivatorElement> &gradient)
{
   // set an initial gradient using some given steps
   // (used in the first iteration)

   double eps2 = fPrecision.Eps2();

   unsigned ix = 0;
   for (auto parameter = parameters.begin(); parameter != parameters.end(); ++parameter, ++ix) {
      // What Minuit2 calls "Error" is stepsize on the ROOT side.
      double werr = parameter->StepSize();

      // Actually, sav in Minuit2 is the external parameter value, so that is
      // what we called var before and var is unnecessary here.
      double sav = parameter->Value();

      // However, we do need var below, so let's calculate it using Ext2int:
      double var = Ext2int(*parameter, sav);

      if (fAlwaysExactlyMimicMinuit2) {
         // this transformation can lose a few bits, but Minuit2 does it too
         sav = Int2ext(*parameter, var);
      }

      double sav2 = sav + werr;

      if (parameter->HasUpperLimit() && sav2 > parameter->UpperLimit()) {
         sav2 = parameter->UpperLimit();
      }

      double var2 = Ext2int(*parameter, sav2);
      double vplu = var2 - var;

      sav2 = sav - werr;

      if (parameter->HasLowerLimit() && sav2 < parameter->LowerLimit()) {
         sav2 = parameter->LowerLimit();
      }

      var2 = Ext2int(*parameter, sav2);
      double vmin = var2 - var;
      double gsmin = 8. * eps2 * (std::abs(var) + eps2);
      // protect against very small step sizes which can cause dirin to zero and then nan values in grd
      double dirin = std::max(0.5 * (std::abs(vplu) + std::abs(vmin)), gsmin);
      double g2 = 2.0 * fUp / (dirin * dirin);
      double gstep = std::max(gsmin, 0.1 * dirin);
      double grd = g2 * dirin;

      if (parameter->IsBound()) {
         gstep = std::min(gstep, 0.5);
      }

      gradient[ix].derivative = grd;
      gradient[ix].second_derivative = g2;
      gradient[ix].step_size = gstep;
   }
}

std::ostream &operator<<(std::ostream &out, const DerivatorElement &value)
{
   return out << "(derivative: " << value.derivative << ", second_derivative: " << value.second_derivative
              << ", step_size: " << value.step_size << ")";
}

} // namespace Minuit2
} // namespace ROOT
