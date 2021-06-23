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
#include <TestStatistics/LikelihoodGradientWrapper.h>
#include "RooMinimizer.h"

namespace RooFit {
namespace TestStatistics {

LikelihoodGradientWrapper::LikelihoodGradientWrapper(std::shared_ptr<RooAbsL> likelihood,
                                                     std::shared_ptr<WrapperCalculationCleanFlags> calculation_is_clean,
                                                     std::size_t /*N_dim*/, RooMinimizer *minimizer)
   : likelihood(std::move(likelihood)), _minimizer(minimizer),
     calculation_is_clean(std::move(calculation_is_clean)) /*, _minimizer_fcn(minimizer_fcn)*/
{
   // Note to future maintainers: take care when storing the minimizer_fcn pointer. The
   // RooAbsMinimizerFcn subclasses may get cloned inside MINUIT, which means the pointer
   // should also somehow be updated in this class.
}

void LikelihoodGradientWrapper::synchronize_with_minimizer(const ROOT::Math::MinimizerOptions & /*options*/) {}

void LikelihoodGradientWrapper::synchronize_parameter_settings(
   const std::vector<ROOT::Fit::ParameterSettings> &parameter_settings)
{
   synchronize_parameter_settings(_minimizer->getMultiGenFcn(), parameter_settings);
}

void LikelihoodGradientWrapper::update_minuit_internal_parameter_values(const std::vector<double>& /*minuit_internal_x*/) {}
void LikelihoodGradientWrapper::update_minuit_external_parameter_values(const std::vector<double>& /*minuit_external_x*/) {}

} // namespace TestStatistics
} // namespace RooFit
