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
#ifndef ROOT_ROOFIT_TESTSTATISTICS_MinuitFcnGrad
#define ROOT_ROOFIT_TESTSTATISTICS_MinuitFcnGrad

#include <Fit/ParameterSettings.h>
#include "ROOT/RMakeUnique.hxx"
#include "Math/IFunction.h" // ROOT::Math::IMultiGradFunction
#include "RooArgList.h"
#include "RooRealVar.h"
#include "TestStatistics/RooAbsL.h"
#include "TestStatistics/LikelihoodWrapper.h"
#include "TestStatistics/LikelihoodGradientWrapper.h"
#include "TestStatistics/LikelihoodJob.h"
#include "TestStatistics/LikelihoodGradientJob.h"
#include "RooAbsMinimizerFcn.h"

// forward declaration
class RooAbsReal;
class RooMinimizer;

namespace RooFit {
namespace TestStatistics {

// -- for communication with wrappers: --
struct WrapperCalculationCleanFlags {
   // indicate whether that part has been calculated since the last parameter update
   bool likelihood = false;
   bool gradient = false;

   void set_all(bool value) {
      likelihood = value;
      gradient = value;
   }
};

class MinuitFcnGrad : public ROOT::Math::IMultiGradFunction, public RooAbsMinimizerFcn {
public:
   // factory
   template <typename LikelihoodWrapperT = RooFit::TestStatistics::LikelihoodJob,
             typename LikelihoodGradientWrapperT = RooFit::TestStatistics::LikelihoodGradientJob>
   static MinuitFcnGrad *create(const std::shared_ptr<RooFit::TestStatistics::RooAbsL> &likelihood,
                                RooMinimizer *context, bool verbose = false);

//   MinuitFcnGrad(const MinuitFcnGrad &other);
   ROOT::Math::IMultiGradFunction *Clone() const override;

   // override to include gradient strategy synchronization:
   Bool_t Synchronize(std::vector<ROOT::Fit::ParameterSettings> &parameter_settings, Bool_t optConst,
                      Bool_t verbose = kFALSE) override;

   // used inside Minuit:
   bool returnsInMinuit2ParameterSpace() const override;

   void setOptimizeConst(Int_t flag) override;

private:
   // IMultiGradFunction overrides necessary for Minuit: DoEval, Gradient, (has)G2ndDerivative and (has)GStepSize
   double DoEval(const double *x) const override;

public:
   void Gradient(const double *x, double *grad) const override;

   // part of IMultiGradFunction interface, used widely both in Minuit and in RooFit:
   unsigned int NDim() const override;

   std::string getFunctionName() const override;
   std::string getFunctionTitle() const override;

   void enable_likelihood_offsetting(bool flag);

private:
   template <typename LikelihoodWrapperT = RooFit::TestStatistics::LikelihoodJob,
             typename LikelihoodGradientWrapperT = RooFit::TestStatistics::LikelihoodGradientJob>
   MinuitFcnGrad(const std::shared_ptr<RooFit::TestStatistics::RooAbsL> &_likelihood, RooMinimizer *context,
                 bool verbose,
                 LikelihoodWrapperT * /* used only for template deduction */ =
                    static_cast<RooFit::TestStatistics::LikelihoodJob *>(nullptr),
                 LikelihoodGradientWrapperT * /* used only for template deduction */ =
                    static_cast<RooFit::TestStatistics::LikelihoodGradientJob *>(nullptr));

   // The following override will not actually be used in this class, so it will throw:
   double DoDerivative(const double *x, unsigned int icoord) const override;

   void optimizeConstantTerms(bool constStatChange, bool constValChange) override;

   bool sync_parameter_values_from_minuit_calls(const double *x, bool minuit_internal) const;

   // members
   std::shared_ptr<LikelihoodWrapper> likelihood;
   std::shared_ptr<LikelihoodGradientWrapper> gradient;

public:
   mutable std::shared_ptr<WrapperCalculationCleanFlags> calculation_is_clean;
private:
   mutable std::vector<double> minuit_internal_x_;
   mutable std::vector<double> minuit_external_x_;
public:
   mutable bool minuit_internal_roofit_x_mismatch_ = false;
};

} // namespace TestStatistics
} // namespace RooFit


// include here to avoid circular dependency issues in class definitions
#include "RooMinimizer.h"


namespace RooFit {
namespace TestStatistics {

template <typename LikelihoodWrapperT, typename LikelihoodGradientWrapperT>
MinuitFcnGrad::MinuitFcnGrad(const std::shared_ptr<RooFit::TestStatistics::RooAbsL> &_likelihood, RooMinimizer *context,
                             bool verbose, LikelihoodWrapperT * /* value unused */,
                             LikelihoodGradientWrapperT * /* value unused */)
   : RooAbsMinimizerFcn(RooArgList(*_likelihood->getParameters()), context, verbose), minuit_internal_x_(NDim(), 0),
     minuit_external_x_(NDim(), 0)
{
   auto parameters = _context->fitter()->Config().ParamsSettings();
   synchronize_parameter_settings(parameters, kTRUE, verbose);

   calculation_is_clean = std::make_shared<WrapperCalculationCleanFlags>();
   likelihood = std::make_shared<LikelihoodWrapperT>(_likelihood, calculation_is_clean/*, _context*/);
   gradient = std::make_shared<LikelihoodGradientWrapperT>(_likelihood, calculation_is_clean, get_nDim(), _context);

   likelihood->synchronize_parameter_settings(parameters);
   gradient->synchronize_parameter_settings(this, parameters);

   // Note: can be different than RooGradMinimizerFcn, where default options are passed (ROOT::Math::MinimizerOptions::DefaultStrategy() and ROOT::Math::MinimizerOptions::DefaultErrorDef())
   likelihood->synchronize_with_minimizer(ROOT::Math::MinimizerOptions());
   gradient->synchronize_with_minimizer(ROOT::Math::MinimizerOptions());
}

// static function
template <typename LikelihoodWrapperT, typename LikelihoodGradientWrapperT>
MinuitFcnGrad *MinuitFcnGrad::create(const std::shared_ptr<RooFit::TestStatistics::RooAbsL>& likelihood,
                                     RooMinimizer *context, bool verbose)
{
   return new MinuitFcnGrad(likelihood, context, verbose, static_cast<LikelihoodWrapperT *>(nullptr),
                            static_cast<LikelihoodGradientWrapperT *>(nullptr));
}

} // namespace TestStatistics
} // namespace RooFit

#endif // ROOT_ROOFIT_TESTSTATISTICS_MinuitFcnGrad
