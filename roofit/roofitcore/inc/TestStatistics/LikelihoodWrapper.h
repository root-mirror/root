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
#ifndef ROOT_ROOFIT_TESTSTATISTICS_LikelihoodWrapper
#define ROOT_ROOFIT_TESTSTATISTICS_LikelihoodWrapper

#include <memory>  // shared_ptr
#include <string>
#include <Fit/ParameterSettings.h>
#include "Math/MinimizerOptions.h"
#include "RooArgSet.h"
#include "RooAbsArg.h"  // enum ConstOpCode

// forward declaration
class RooMinimizer;

namespace RooFit {
namespace TestStatistics {

// forward declaration
class RooAbsL;
struct WrapperCalculationCleanFlags;

enum class LikelihoodType {
   unbinned,
   binned,
   subsidiary,
   sum
};

// Previously, offsetting was only implemented for RooNLLVar components of a likelihood,
// not for RooConstraintSum terms. To emulate this behavior, use OffsettingMode::legacy. To
// also offset the RooSubsidiaryL component (equivalent of RooConstraintSum) of RooSumL
// likelihoods, use OffsettingMode::full.
enum class OffsettingMode {
   legacy,
   full
};

class LikelihoodWrapper {
public:
   LikelihoodWrapper(std::shared_ptr<RooAbsL> likelihood, std::shared_ptr<WrapperCalculationCleanFlags> calculation_is_clean/*, RooMinimizer* minimizer*/);
   virtual ~LikelihoodWrapper() = default;
   virtual LikelihoodWrapper* clone() const = 0;

   virtual void evaluate() = 0;
   virtual double return_result() const = 0;

   // synchronize minimizer settings with calculators in child classes
   virtual void synchronize_with_minimizer(const ROOT::Math::MinimizerOptions & options);
   virtual void synchronize_parameter_settings(const std::vector<ROOT::Fit::ParameterSettings> &parameter_settings);
   // Minuit passes in parameter values that may not conform to RooFit internal standards (like applying range clipping),
   // but that the specific calculator does need. This function can be implemented to receive these Minuit-internal values:
   virtual void update_minuit_internal_parameter_values(const std::vector<double>& minuit_internal_x);
   virtual void update_minuit_external_parameter_values(const std::vector<double>& minuit_external_x);

   // necessary from MinuitFcnGrad to reach likelihood properties:
   void constOptimizeTestStatistic(RooAbsArg::ConstOpCode opcode, bool doAlsoTrackingOpt);

   double defaultErrorLevel() const;

   virtual std::string GetName() const;
   virtual std::string GetTitle() const;

   virtual bool is_offsetting() const;
   virtual void enable_offsetting(bool flag);
   void set_offsetting_mode(OffsettingMode mode);
   double offset();
   double offset_carry();
   void set_apply_weight_squared(bool flag);

protected:
   std::shared_ptr<RooAbsL> likelihood_;
//   RooMinimizer* minimizer_;
//   RooAbsMinimizerFcn* minimizer_fcn_;
   std::shared_ptr<WrapperCalculationCleanFlags> calculation_is_clean_;

   bool do_offset_ = false;
   double offset_ = 0;
   double offset_carry_ = 0;
   double offset_save_ = 0;      //!
   double offset_carry_save_ = 0; //!
   OffsettingMode offsetting_mode_ = OffsettingMode::legacy;
   void apply_offsetting(double &current_value, double &carry);
   void swap_offsets();
};

}
}

#endif // ROOT_ROOFIT_TESTSTATISTICS_LikelihoodWrapper
