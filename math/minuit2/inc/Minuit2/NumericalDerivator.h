// @(#)root/mathcore:$Id$
// Authors: L. Moneta, J.T. Offermann, E.G.P. Bos    2013-2018
//
/**********************************************************************
 *                                                                    *
 * Copyright (c) 2013 , LCG ROOT MathLib Team                         *
 *                                                                    *
 **********************************************************************/
/*
 * NumericalDerivator.h
 *
 *  Original version created on: Aug 14, 2013
 *      Authors: L. Moneta, J. T. Offermann
 *  Modified version created on: Sep 27, 2017
 *      Author: E. G. P. Bos
 */

#ifndef ROOT_Minuit2_NumericalDerivator
#define ROOT_Minuit2_NumericalDerivator

#include <Math/IFunctionfwd.h>

#include <vector>
#include "Fit/ParameterSettings.h"
#include "Minuit2/SinParameterTransformation.h"
#include "Minuit2/SqrtUpParameterTransformation.h"
#include "Minuit2/SqrtLowParameterTransformation.h"
#include "Minuit2/MnMachinePrecision.h"

namespace ROOT {
namespace Minuit2 {

// Holds all necessary derivatives and associated numbers (per parameter) used in the NumericalDerivator class.
struct DerivatorElement {
   double derivative;
   double second_derivative;
   double step_size;
};

class NumericalDerivator {
public:
   explicit NumericalDerivator(bool always_exactly_mimic_minuit2 = true);
   NumericalDerivator(const NumericalDerivator &other);
   NumericalDerivator(double step_tolerance, double grad_tolerance, unsigned int ncycles, double error_level,
                      bool always_exactly_mimic_minuit2 = true);

   void SetupDifferentiate(const ROOT::Math::IBaseFunctionMultiDim *function, const double *cx,
                           const std::vector<ROOT::Fit::ParameterSettings> &parameters);
   std::vector<DerivatorElement> Differentiate(const ROOT::Math::IBaseFunctionMultiDim *function, const double *x,
                                               const std::vector<ROOT::Fit::ParameterSettings> &parameters,
                                               const std::vector<DerivatorElement> &previous_gradient);

   DerivatorElement PartialDerivative(const ROOT::Math::IBaseFunctionMultiDim *function, const double *x,
                                      const std::vector<ROOT::Fit::ParameterSettings> &parameters,
                                      unsigned int i_component, DerivatorElement previous);
   DerivatorElement FastPartialDerivative(const ROOT::Math::IBaseFunctionMultiDim *function,
                                          const std::vector<ROOT::Fit::ParameterSettings> &parameters,
                                          unsigned int i_component, const DerivatorElement &previous);
   DerivatorElement operator()(const ROOT::Math::IBaseFunctionMultiDim *function, const double *x,
                               const std::vector<ROOT::Fit::ParameterSettings> &parameters, unsigned int i_component,
                               const DerivatorElement &previous);

   double GetValue() const { return fVal; }
   inline void SetStepTolerance(double value) { fStepTolerance = value; }
   inline void SetGradTolerance(double value) { fGradTolerance = value; }
   inline void SetNCycles(unsigned int value) { fNCycles = value; }
   inline void SetErrorLevel(double value) { fUp = value; }

   double Int2ext(const ROOT::Fit::ParameterSettings &parameter, double val) const;
   double Ext2int(const ROOT::Fit::ParameterSettings &parameter, double val) const;
   double DInt2Ext(const ROOT::Fit::ParameterSettings &parameter, double val) const;

   void SetInitialGradient(const std::vector<ROOT::Fit::ParameterSettings> &parameters,
                           std::vector<DerivatorElement> &gradient);

   inline bool AlwaysExactlyMimicMinuit2() const { return fAlwaysExactlyMimicMinuit2; }
   inline void SetAlwaysExactlyMimicMinuit2(bool flag) { fAlwaysExactlyMimicMinuit2 = flag; }

private:
   double fStepTolerance = 0.5;
   double fGradTolerance = 0.1;
   double fUp = 1;
   double fVal = 0;

   std::vector<double> fVx;
   std::vector<double> fVxExternal;
   std::vector<double> fVxFValCache;
   double fDfmin;
   double fVrysml;

   // MODIFIED: Minuit2 determines machine precision in a slightly different way than
   // std::numeric_limits<double>::epsilon()). We go with the Minuit2 one.
   ROOT::Minuit2::MnMachinePrecision fPrecision;

   ROOT::Minuit2::SinParameterTransformation fDoubleLimTrafo;
   ROOT::Minuit2::SqrtUpParameterTransformation fUpperLimTrafo;
   ROOT::Minuit2::SqrtLowParameterTransformation fLowerLimTrafo;

   unsigned int fNCycles = 2;
   bool fAlwaysExactlyMimicMinuit2;

};

std::ostream &operator<<(std::ostream &out, const DerivatorElement &value);

} // namespace Minuit2
} // namespace ROOT

#endif // ROOT_Minuit2_NumericalDerivator
