/*
 * Project: RooFit
 * Authors:
 *   PB, Patrick Bos, Netherlands eScience Center, p.bos@esciencecenter.nl
 *
 * Copyright (c) 2016-2020, Netherlands eScience Center
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 */
#ifndef ROOT_ROOFIT_TESTSTATISTICS_LikelihoodGradientWrapper
#define ROOT_ROOFIT_TESTSTATISTICS_LikelihoodGradientWrapper

#include <vector>
#include <memory> // shared_ptr
#include <Fit/ParameterSettings.h>
#include <Math/IFunctionfwd.h>
#include "Math/MinimizerOptions.h"

// forward declaration
class RooMinimizer;

namespace RooFit {
namespace TestStatistics {

// forward declaration
class RooAbsL;
struct WrapperCalculationCleanFlags;

class LikelihoodGradientWrapper {
public:
   LikelihoodGradientWrapper(std::shared_ptr<RooAbsL> likelihood, std::shared_ptr<WrapperCalculationCleanFlags> calculation_is_clean, std::size_t N_dim, RooMinimizer* minimizer);
   virtual ~LikelihoodGradientWrapper() = default;
   virtual LikelihoodGradientWrapper* clone() const = 0;

   virtual void fill_gradient(double *grad) = 0;
   virtual void fill_second_derivative(double *g2) = 0;
   virtual void fill_step_size(double *gstep) = 0;

   // synchronize minimizer settings with calculators in child classes
   virtual void synchronize_with_minimizer(const ROOT::Math::MinimizerOptions &options);
   virtual void synchronize_parameter_settings(const std::vector<ROOT::Fit::ParameterSettings> &parameter_settings);
   virtual void synchronize_parameter_settings(ROOT::Math::IMultiGenFunction* function, const std::vector<ROOT::Fit::ParameterSettings> &parameter_settings) = 0;
   // Minuit passes in parameter values that may not conform to RooFit internal standards (like applying range clipping),
   // but that the specific calculator does need. This function can be implemented to receive these Minuit-internal values:
   virtual void update_minuit_internal_parameter_values(const std::vector<double>& minuit_internal_x);
   virtual void update_minuit_external_parameter_values(const std::vector<double>& minuit_external_x);

   // completely depends on the implementation, so pure virtual
   virtual bool uses_minuit_internal_values() = 0;

protected:
   std::shared_ptr<RooAbsL> likelihood;
   RooMinimizer* _minimizer;
   std::shared_ptr<WrapperCalculationCleanFlags> calculation_is_clean;
//   MinuitFcnGrad* _minimizer_fcn;
};

} // namespace TestStatistics
} // namespace RooFit

#endif // ROOT_ROOFIT_TESTSTATISTICS_LikelihoodGradientWrapper
