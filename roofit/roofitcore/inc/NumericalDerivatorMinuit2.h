// @(#)root/mathcore:$Id$
// Authors: L. Moneta, J.T. Offermann, E.G.P. Bos    2013-2018
//
/**********************************************************************
 *                                                                    *
 * Copyright (c) 2013 , LCG ROOT MathLib Team                         *
 * Copyright (c) 2017 Patrick Bos, Netherlands eScience Center        *
 *                                                                    *
 **********************************************************************/
/*
 * NumericalDerivatorMinuit2.h
 *
 *  Original version (NumericalDerivator) created on: Aug 14, 2013
 *      Authors: L. Moneta, J. T. Offermann
 *  Modified version (NumericalDerivatorMinuit2) created on: Sep 27, 2017
 *      Author: E. G. P. Bos
 */

#ifndef RooFit_NumericalDerivatorMinuit2
#define RooFit_NumericalDerivatorMinuit2

#ifndef ROOT_Math_IFunctionfwd
#include <Math/IFunctionfwd.h>
#endif

#include <vector>
#include "Fit/ParameterSettings.h"
#include "Minuit2/SinParameterTransformation.h"
#include "Minuit2/SqrtUpParameterTransformation.h"
#include "Minuit2/SqrtLowParameterTransformation.h"
#include "Minuit2/MnMachinePrecision.h"

namespace RooFit {

// Holds all necessary derivatives and associated numbers (per parameter) used in the NumericalDerivatorMinuit2 class.
struct MinuitDerivatorElement {
   double derivative;
   double second_derivative;
   double step_size;
};

class NumericalDerivatorMinuit2 {
public:
   explicit NumericalDerivatorMinuit2(bool always_exactly_mimic_minuit2 = true);
   NumericalDerivatorMinuit2(const NumericalDerivatorMinuit2 &other);
   NumericalDerivatorMinuit2(double step_tolerance, double grad_tolerance, unsigned int ncycles, double error_level,
                             bool always_exactly_mimic_minuit2 = true);
   virtual ~NumericalDerivatorMinuit2();

   void setup_differentiate(const ROOT::Math::IBaseFunctionMultiDim *function, const double *cx,
                            const std::vector<ROOT::Fit::ParameterSettings> &parameters);
   std::vector<MinuitDerivatorElement>
   Differentiate(const ROOT::Math::IBaseFunctionMultiDim *function, const double *x,
                 const std::vector<ROOT::Fit::ParameterSettings> &parameters,
                 const std::vector<MinuitDerivatorElement>& previous_gradient);

   MinuitDerivatorElement
   partial_derivative(const ROOT::Math::IBaseFunctionMultiDim *function, const double *x,
                      const std::vector<ROOT::Fit::ParameterSettings> &parameters, unsigned int i_component, MinuitDerivatorElement previous);
   MinuitDerivatorElement fast_partial_derivative(const ROOT::Math::IBaseFunctionMultiDim *function,
                                const std::vector<ROOT::Fit::ParameterSettings> &parameters, unsigned int i_component, const MinuitDerivatorElement& previous);
   MinuitDerivatorElement operator()(const ROOT::Math::IBaseFunctionMultiDim *function, const double *x,
                                                 const std::vector<ROOT::Fit::ParameterSettings> &parameters,
                                                 unsigned int i_component, const MinuitDerivatorElement& previous);

   double GetFValue() const { return fVal; }
   //   const double *GetG2() const { return fG.G2().Data(); }
   void SetStepTolerance(double value);
   void SetGradTolerance(double value);
   void SetNCycles(int value);

   double Int2ext(const ROOT::Fit::ParameterSettings &parameter, double val) const;
   double Ext2int(const ROOT::Fit::ParameterSettings &parameter, double val) const;

   void SetInitialGradient(const ROOT::Math::IBaseFunctionMultiDim *function,
                           const std::vector<ROOT::Fit::ParameterSettings> &parameters,
                           std::vector<MinuitDerivatorElement>& gradient);

   void set_step_tolerance(double step_tolerance);
   void set_grad_tolerance(double grad_tolerance);
   void set_ncycles(unsigned int ncycles);
   void set_error_level(double error_level);

private:
   double fStepTolerance = 0.5;
   double fGradTolerance = 0.1;
   unsigned int fNCycles = 2;
   double Up = 1;
   double fVal = 0;

   std::vector<double> vx, vx_external;
   double dfmin;
   double vrysml;

   // MODIFIED: Minuit2 determines machine precision in a slightly different way than
   // std::numeric_limits<double>::epsilon()). We go with the Minuit2 one.
   ROOT::Minuit2::MnMachinePrecision precision;

   ROOT::Minuit2::SinParameterTransformation fDoubleLimTrafo;
   ROOT::Minuit2::SqrtUpParameterTransformation fUpperLimTrafo;
   ROOT::Minuit2::SqrtLowParameterTransformation fLowerLimTrafo;

private:
   bool _always_exactly_mimic_minuit2;

public:
   bool always_exactly_mimic_minuit2() const;
   void set_always_exactly_mimic_minuit2(bool flag);

private:
   std::vector<double> vx_fVal_cache;
//#ifndef NDEBUG
//   std::size_t fVal_eval_counter = 0; //!
//#endif
};

std::ostream& operator<<(std::ostream& out, const MinuitDerivatorElement value);

} // namespace RooFit

#endif /* NumericalDerivatorMinuit2_H_ */