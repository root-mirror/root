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
#include <TestStatistics/LikelihoodWrapper.h>
#include <TestStatistics/RooAbsL.h> // need complete type for likelihood->...
#include <TestStatistics/MinuitFcnGrad.h>
#include <TestStatistics/RooUnbinnedL.h>
#include <TestStatistics/RooSumL.h> // need complete type for dynamic cast

namespace RooFit {
namespace TestStatistics {

LikelihoodWrapper::LikelihoodWrapper(std::shared_ptr<RooAbsL> likelihood,
                                     std::shared_ptr<WrapperCalculationCleanFlags> calculation_is_clean/*,
                                     RooMinimizer *minimizer*/)
   : likelihood_(std::move(likelihood)),/* minimizer_(minimizer),*/
     calculation_is_clean_(std::move(calculation_is_clean)) /*, minimizer_fcn_(minimizer_fcn)*/
{
   // Note to future maintainers: take care when storing the minimizer_fcn pointer. The
   // RooAbsMinimizerFcn subclasses may get cloned inside MINUIT, which means the pointer
   // should also somehow be updated in this class.
}

void LikelihoodWrapper::synchronize_with_minimizer(const ROOT::Math::MinimizerOptions & /*options*/) {}

void LikelihoodWrapper::constOptimizeTestStatistic(RooAbsArg::ConstOpCode opcode, bool doAlsoTrackingOpt)
{
   likelihood_->constOptimizeTestStatistic(opcode, doAlsoTrackingOpt);
}

void LikelihoodWrapper::synchronize_parameter_settings(
   const std::vector<ROOT::Fit::ParameterSettings> & /*parameter_settings*/)
{
}

std::string LikelihoodWrapper::GetName() const
{
   return likelihood_->GetName();
}
std::string LikelihoodWrapper::GetTitle() const
{
   return likelihood_->GetTitle();
}
double LikelihoodWrapper::defaultErrorLevel() const
{
   return likelihood_->defaultErrorLevel();
}
bool LikelihoodWrapper::is_offsetting() const
{
   return do_offset_;
}

double LikelihoodWrapper::offset()
{
   return offset_;
}

double LikelihoodWrapper::offset_carry()
{
   return offset_carry_;
}

void LikelihoodWrapper::enable_offsetting(bool flag)
{
   do_offset_ = flag;
   // Clear offset if feature is disabled so that it is recalculated next time it is enabled
   if (!do_offset_) {
      offset_ = 0;
      offset_carry_ = 0;
   }
}

void LikelihoodWrapper::set_offsetting_mode(OffsettingMode mode)
{
   offsetting_mode_ = mode;
   if (is_offsetting()) {
      oocoutI(static_cast<RooAbsArg *>(nullptr), Minimization) << "LikelihoodWrapper::set_offsetting_mode(" << GetName() << "): changed offsetting mode while offsetting was enabled; resetting offset values" << std::endl;
      offset_ = 0;
      offset_carry_ = 0;
   }
}

void LikelihoodWrapper::apply_offsetting(double &current_value, double &carry)
{
   if (do_offset_) {

      // If no offset is stored enable this feature now
      if (offset_ == 0 && current_value != 0) {
         offset_ = current_value;
         offset_carry_ = carry;
         if (offsetting_mode_ == OffsettingMode::legacy) {
            auto sum_likelihood = dynamic_cast<RooSumL*>(likelihood_.get());
            if (sum_likelihood != nullptr) {
               double subsidiary_value, subsidiary_carry;
               std::tie(subsidiary_value, subsidiary_carry) = sum_likelihood->get_subsidiary_value();
               // "undo" the addition of the subsidiary value to emulate legacy behavior
               offset_ -= subsidiary_value;
               offset_carry_ -= subsidiary_carry;
               // then add 0 in Kahan summation way to make sure the carry gets taken up into the value if it should be
               double y = 0 - offset_carry_;
               double t = offset_ + y;
               offset_carry_ = (t - offset_) - y;
               offset_ = t;
               // also set carry to this value, again to emulate legacy behavior
               carry = offset_carry_;
            }
         }
         oocoutI(static_cast<RooAbsArg *>(nullptr), Minimization)
            << "LikelihoodWrapper::apply_offsetting(" << GetName() << "): Likelihood offset now set to " << offset_
            << std::endl;
      }

      // Subtract offset
      // old method:
//      {
//         double y = -offset_ - (carry + offset_carry_);
//         double t = current_value + y;
//         carry = (t - current_value) - y;
//         current_value = t;
//      }
      // Jonas method:
      {
         double new_value = current_value - offset_;
         double new_carry = carry - offset_carry_;
         // then add 0 in Kahan summation way to make sure the carry gets taken up into the value if it should be
         double y = 0 - new_carry;
         double t = new_value + y;
         carry = (t - new_value) - y;
         current_value = t;
      }
   }
}

/// When calculating an unbinned likelihood with square weights applied, a different offset
/// is necessary. Similar situations may ask for a separate offset as well. This function
/// switches between the two sets of offset values.
void LikelihoodWrapper::swap_offsets()
{
   std::swap(offset_, offset_save_);
   std::swap(offset_carry_, offset_carry_save_);
}

void LikelihoodWrapper::set_apply_weight_squared(bool flag)
{
   RooUnbinnedL *unbinned_likelihood = dynamic_cast<RooUnbinnedL*>(likelihood_.get());
   if (unbinned_likelihood == nullptr) {
      throw std::logic_error("LikelihoodWrapper::set_apply_weight_squared can only be used on unbinned likelihoods, but the wrapped likelihood_ member is not a RooUnbinnedL!");
   }
   bool flag_was_changed = unbinned_likelihood->set_apply_weight_squared(flag);

   if (flag_was_changed) {
      swap_offsets();
   }
}

void LikelihoodWrapper::update_minuit_internal_parameter_values(const std::vector<double>& /*minuit_internal_x*/) {}
void LikelihoodWrapper::update_minuit_external_parameter_values(const std::vector<double>& /*minuit_external_x*/) {}

} // namespace TestStatistics
} // namespace RooFit
