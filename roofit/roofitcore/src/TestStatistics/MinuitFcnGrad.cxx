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

#include "RooMsgService.h"
#include "RooAbsPdf.h"
#include "TestStatistics/MinuitFcnGrad.h"
#include "RooMinimizer.h"

#define DEBUG_STREAM(var) << " " #var "=" << var
#include <sys/types.h>
#include <unistd.h>

namespace RooFit {
namespace TestStatistics {

//MinuitFcnGrad::MinuitFcnGrad(const MinuitFcnGrad &other)
//   : RooAbsMinimizerFcn(other), likelihood(other.likelihood->clone()), gradient(other.gradient->clone()){};

// IMultiGradFunction overrides necessary for Minuit: DoEval, Gradient, G2ndDerivative and GStepSize
// The likelihood and gradient wrappers do the actual calculations.

double MinuitFcnGrad::DoEval(const double *x) const
{
//   auto get_time = []() {
//      return std::chrono::duration_cast<std::chrono::nanoseconds>(
//         std::chrono::high_resolution_clock::now().time_since_epoch())
//         .count();
//   };
//   decltype(get_time()) t1 = get_time(), t2 = 0, t3 = 0, t4 = 0;

   Bool_t parameters_changed = sync_parameter_values_from_minuit_calls(x, false);
//   t2 = get_time();

//   std::cout << "MinuitFcnGrad::DoEval @ PID" << getpid() << ": " DEBUG_STREAM(parameters_changed);

   // Calculate the function for these parameters
//   RooAbsReal::setHideOffset(kFALSE);
   likelihood->evaluate();
//   t3 = get_time();
   double fvalue = likelihood->return_result();
//   t4 = get_time();
   calculation_is_clean->likelihood = true;
//   RooAbsReal::setHideOffset(kTRUE);

//   std::cout DEBUG_STREAM(fvalue) << std::endl;

//   printf("wallclock [worker] DoEval parts: 1. %f 2. %f 3. %f\n", (t2 - t1) / 1.e9, (t3 - t2) / 1.e9, (t4 - t3) / 1.e9);

   if (!parameters_changed) {
      return fvalue;
   }

   if (!std::isfinite(fvalue) || RooAbsReal::numEvalErrors() > 0 || fvalue > 1e30) {

      if (_printEvalErrors >= 0) {

         if (_doEvalErrorWall) {
            oocoutW(static_cast<RooAbsArg *>(nullptr), Eval)
               << "RooGradMinimizerFcn: Minimized function has error status." << std::endl
               << "Returning maximum FCN so far (" << _maxFCN
               << ") to force MIGRAD to back out of this region. Error log follows" << std::endl;
         } else {
            oocoutW(static_cast<RooAbsArg *>(nullptr), Eval)
               << "RooGradMinimizerFcn: Minimized function has error status but is ignored" << std::endl;
         }

         TIterator *iter = _floatParamList->createIterator();
         RooRealVar *var;
         Bool_t first(kTRUE);
         ooccoutW(static_cast<RooAbsArg *>(nullptr), Eval) << "Parameter values: ";
         while ((var = (RooRealVar *)iter->Next())) {
            if (first) {
               first = kFALSE;
            } else
               ooccoutW(static_cast<RooAbsArg *>(nullptr), Eval) << ", ";
            ooccoutW(static_cast<RooAbsArg *>(nullptr), Eval) << var->GetName() << "=" << var->getVal();
         }
         delete iter;
         ooccoutW(static_cast<RooAbsArg *>(nullptr), Eval) << std::endl;

         RooAbsReal::printEvalErrors(ooccoutW(static_cast<RooAbsArg *>(nullptr), Eval), _printEvalErrors);
         ooccoutW(static_cast<RooAbsArg *>(nullptr), Eval) << std::endl;
      }

      if (_doEvalErrorWall) {
         fvalue = _maxFCN + 1;
      }

      RooAbsReal::clearEvalErrorLog();
      _numBadNLL++;
   } else if (fvalue > _maxFCN) {
      _maxFCN = fvalue;
   }

   // Optional logging
   if (_verbose) {
      std::cout << "\nprevFCN" << (likelihood->is_offsetting() ? "-offset" : "") << " = " << std::setprecision(10)
                << fvalue << std::setprecision(4) << "  ";
      std::cout.flush();
   }

   _evalCounter++;
   //#ifndef NDEBUG
   //  std::cout << "RooGradMinimizerFcn " << this << " evaluations (in DoEval): " << _evalCounter <<
   //  std::endl;
   //#endif
   return fvalue;
}

/// Minuit calls (via FcnAdapters etc) DoEval or Gradient/G2ndDerivative/GStepSize with a set of parameters x.
/// This function syncs these values to the proper places in RooFit.
///
/// The first twist, and reason this function is more complicated than one may imagine, is that Minuit internally uses a
/// transformed parameter space to account for parameter boundaries. Whether we receive these Minuit "internal"
/// parameter values or "regular"/untransformed RooFit parameter space values depends on the situation.
/// - The values that arrive here via DoEval are always "normal" parameter values, since Minuit transforms these
///   back into regular space before passing to DoEval (see MnUserFcn::operator() which wraps the Fcn(Gradient)Base
///   in ModularFunctionMinimizer::Minimize and is used for direct function calls from that point on in the minimizer).
///   These can thus always be safely synced with this function's RooFit parameters using SetPdfParamVal.
/// - The values that arrive here via Gradient/G2ndDerivative/GStepSize will be in internal coordinates if that is
///   what this class expects, and indeed this is the case for MinuitFcnGrad's current implementation. This is
///   communicated to Minuit via MinuitFcnGrad::returnsInMinuit2ParameterSpace. Inside Minuit, that function determines
///   whether this class's gradient calculator is wrapped inside a AnalyticalGradientCalculator, to which Minuit passes
///   "external" parameter values, or as an ExternalInternalGradientCalculator, which gets "internal" parameter values.
///   Long story short: when MinuitFcnGrad::returnsInMinuit2ParameterSpace() returns true, Minuit will pass "internal"
///   values to Gradient/G2ndDerivative/GStepSize. These cannot be synced with this function's RooFit parameters using
///   SetPdfParamVal, unless a manual transformation step is performed in advance. However, they do need to be passed
///   on to the gradient calculator, since indeed we expect values there to be in "internal" space. However, this is
///   calculator dependent. Note that in the current MinuitFcnGrad implementation we do not actually allow for
///   calculators in "external" (i.e. regular RooFit parameter space) values, since
///   MinuitFcnGrad::returnsInMinuit2ParameterSpace is hardcoded to true. This should in a future version be changed so
///   that the calculator (the wrapper) is queried for this information.
/// Because some gradient calculators may also use the regular RooFit parameters (e.g. for calculating the likelihood's
/// value itself), this information is also passed on to the gradient wrapper. Vice versa, when updated "internal"
/// parameters are passed to Gradient/G2ndDerivative/GStepSize, the likelihood may be affected as well. Even though a
/// transformation from internal to "external" may be necessary before the values can be used, the likelihood can at
/// least log that its parameter values are possibly no longer in sync with those of the gradient.
///
/// The second twist is that the Minuit external parameters may still be different from the ones used in RooFit. This
/// happens when Minuit tries out values that lay outside the RooFit parameter's range(s). RooFit's setVal (called
/// inside SetPdfParamVal) then clips the RooAbsArg's value to one of the range limits, instead of setting it to the
/// value Minuit intended. When this happens, i.e. sync_parameter_values_from_minuit_calls is called with
/// minuit_internal = false and the values do not match the previous values stored in minuit_internal_x_ *but* the
/// values after SetPdfParamVal did not get set to the intended value, the minuit_internal_roofit_x_mismatch_ flag is
/// set. This information can be used by calculators, if desired, for instance when a calculator does not want to make
/// use of the range information in the RooAbsArg parameters.
bool MinuitFcnGrad::sync_parameter_values_from_minuit_calls(const double *x, bool minuit_internal) const
{
   bool a_parameter_has_been_updated = false;
   if (minuit_internal) {
      if (!returnsInMinuit2ParameterSpace()) {
         throw std::logic_error("Updating Minuit-internal parameters only makes sense for (gradient) calculators that are defined in Minuit-internal parameter space.");
      }

      for (std::size_t ix = 0; ix < NDim(); ++ix) {
         bool parameter_changed = (x[ix] != minuit_internal_x_[ix]);
         if (parameter_changed) {
            minuit_internal_x_[ix] = x[ix];
         }
         a_parameter_has_been_updated |= parameter_changed;
      }

      if(a_parameter_has_been_updated) {
         calculation_is_clean->set_all(false);
         likelihood->update_minuit_internal_parameter_values(minuit_internal_x_);
         gradient->update_minuit_internal_parameter_values(minuit_internal_x_);
      }
   } else {
      bool a_parameter_is_mismatched = false;

      for (std::size_t ix = 0; ix < NDim(); ++ix) {
         // Note: the return value of SetPdfParamVal does not always mean that the parameter's value in the RooAbsReal changed since last
         // time! If the value was out of range bin, setVal was still called, but the value was not updated.
         SetPdfParamVal(ix, x[ix]);
         minuit_external_x_[ix] = x[ix];
         // The above is why we need minuit_external_x_. The minuit_external_x_ vector can also be passed to
         // LikelihoodWrappers, if needed, but typically they will make use of the RooFit parameters directly. However,
         // we log in the flag below whether they are different so that calculators can use this information.
         bool parameter_changed = (x[ix] != minuit_external_x_[ix]);
         a_parameter_has_been_updated |= parameter_changed;
         a_parameter_is_mismatched |= (((RooRealVar *)_floatParamList->at(ix))->getVal() != minuit_external_x_[ix]);
      }

      minuit_internal_roofit_x_mismatch_ = a_parameter_is_mismatched;

      if(a_parameter_has_been_updated) {
         calculation_is_clean->set_all(false);
         likelihood->update_minuit_external_parameter_values(minuit_external_x_);
         gradient->update_minuit_external_parameter_values(minuit_external_x_);
      }
   }
   return a_parameter_has_been_updated;
}


void MinuitFcnGrad::Gradient(const double *x, double *grad) const
{
   printf("Calling MinuitFcnGrad::Gradient\n");
   sync_parameter_values_from_minuit_calls(x, returnsInMinuit2ParameterSpace());
   gradient->fill_gradient(grad);
}

ROOT::Math::IMultiGradFunction *MinuitFcnGrad::Clone() const
{
   return new MinuitFcnGrad(*this);
}

double MinuitFcnGrad::DoDerivative(const double * /*x*/, unsigned int /*icoord*/) const
{
   throw std::runtime_error("MinuitFcnGrad::DoDerivative is not implemented, please use Gradient instead.");
}

unsigned int MinuitFcnGrad::NDim() const
{
   return _nDim;
}

bool MinuitFcnGrad::returnsInMinuit2ParameterSpace() const
{
   return gradient->uses_minuit_internal_values();
}

void MinuitFcnGrad::optimizeConstantTerms(bool constStatChange, bool constValChange)
{
   if (constStatChange) {

      RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::CollectErrors);

      oocoutI(static_cast<RooAbsArg *>(nullptr), Eval)
         << "MinuitFcnGrad::synchronize: set of constant parameters changed, rerunning const optimizer" << std::endl;
      likelihood->constOptimizeTestStatistic(RooAbsArg::ConfigChange, true);
   } else if (constValChange) {
      oocoutI(static_cast<RooAbsArg *>(nullptr), Eval)
         << "MinuitFcnGrad::synchronize: constant parameter values changed, rerunning const optimizer" << std::endl;
      likelihood->constOptimizeTestStatistic(RooAbsArg::ValueChange, true);
   }

   RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::PrintErrors);
}

Bool_t
MinuitFcnGrad::Synchronize(std::vector<ROOT::Fit::ParameterSettings> &parameters, Bool_t optConst, Bool_t verbose)
{
   Bool_t returnee = synchronize_parameter_settings(parameters, optConst, verbose);
   likelihood->synchronize_parameter_settings(parameters);
   gradient->synchronize_parameter_settings(parameters);

   likelihood->synchronize_with_minimizer(_context->fitter()->Config().MinimizerOptions());
   gradient->synchronize_with_minimizer(_context->fitter()->Config().MinimizerOptions());
   return returnee;
}

std::string MinuitFcnGrad::getFunctionName() const
{
   return likelihood->GetName();
}

std::string MinuitFcnGrad::getFunctionTitle() const
{
   return likelihood->GetTitle();
}

void MinuitFcnGrad::setOptimizeConst(Int_t flag)
{
   RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::CollectErrors);

   if (_optConst && !flag) {
      if (_context->getPrintLevel() > -1)
         oocoutI(_context, Minimization) << "MinuitFcnGrad::setOptimizeConst: deactivating const optimization"
                                         << std::endl;
      likelihood->constOptimizeTestStatistic(RooAbsArg::DeActivate, true);
      _optConst = flag;
   } else if (!_optConst && flag) {
      if (_context->getPrintLevel() > -1)
         oocoutI(_context, Minimization) << "MinuitFcnGrad::setOptimizeConst: activating const optimization"
                                         << std::endl;
      likelihood->constOptimizeTestStatistic(RooAbsArg::Activate, flag > 1);
      _optConst = flag;
   } else if (_optConst && flag) {
      if (_context->getPrintLevel() > -1)
         oocoutI(_context, Minimization) << "MinuitFcnGrad::setOptimizeConst: const optimization already active"
                                         << std::endl;
   } else {
      if (_context->getPrintLevel() > -1)
         oocoutI(_context, Minimization) << "MinuitFcnGrad::setOptimizeConst: const optimization wasn't active"
                                         << std::endl;
   }

   RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::PrintErrors);
}

void MinuitFcnGrad::enable_likelihood_offsetting(bool flag) {
   likelihood->enable_offsetting(flag);
}

} // namespace TestStatistics
} // namespace RooFit
