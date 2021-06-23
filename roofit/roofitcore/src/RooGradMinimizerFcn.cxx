/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 * @(#)root/roofitcore:$Id$
 * Authors:                                                                  *
 *   AL, Alfio Lazzaro,   INFN Milan,        alfio.lazzaro@mi.infn.it        *
 *   PB, Patrick Bos,     NL eScience Center, p.bos@esciencecenter.nl        *
 *   VC, Vince Croft,     DIANA / NYU,        vincent.croft@cern.ch          *
 *                                                                           *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

#ifndef __ROOFIT_NOROOGRADMINIMIZER

//////////////////////////////////////////////////////////////////////////////
//
// RooGradMinimizerFcn is am interface class to the ROOT::Math function
// for minization. See RooGradMinimizer.cxx for more information.
//

#include <iostream>

#include "RooFit.h"

#include "Riostream.h"

#include "TIterator.h"
#include "TClass.h"

#include "RooAbsArg.h"
#include "RooAbsPdf.h"
#include "RooArgSet.h"
#include "RooRealVar.h"
#include "RooAbsRealLValue.h"
#include "RooMsgService.h"

#include "RooMinimizer.h"
#include "RooGradMinimizerFcn.h"

#include "Fit/Fitter.h"
#include "Math/Minimizer.h"
//#include "Minuit2/MnMachinePrecision.h"

#include <algorithm> // std::equal

RooGradMinimizerFcn::RooGradMinimizerFcn(RooAbsReal *funct, RooMinimizer *context, bool verbose)
   : RooAbsMinimizerFcn(RooArgList(*funct->getParameters(RooArgSet())), context, verbose),
     _grad(get_nDim()), _grad_params(get_nDim()), _funct(funct),
     has_been_calculated(get_nDim())
{
   // TODO: added "parameters" after rewrite in april 2020, check if correct
   auto parameters = _context->fitter()->Config().ParamsSettings();
   synchronize_parameter_settings(parameters, kTRUE, verbose);
   synchronize_gradient_parameter_settings(parameters);
   set_strategy(ROOT::Math::MinimizerOptions::DefaultStrategy());
   set_error_level(ROOT::Math::MinimizerOptions::DefaultErrorDef());
}

RooGradMinimizerFcn::RooGradMinimizerFcn(const RooGradMinimizerFcn &other)
   : RooAbsMinimizerFcn(other), _grad(other._grad), _grad_params(other._grad_params), _gradf(other._gradf), _funct(other._funct),
     has_been_calculated(other.has_been_calculated), none_have_been_calculated(other.none_have_been_calculated)
{
}

ROOT::Math::IMultiGradFunction *RooGradMinimizerFcn::Clone() const
{
   return new RooGradMinimizerFcn(*this);
}

void RooGradMinimizerFcn::synchronize_gradient_parameter_settings(
   std::vector<ROOT::Fit::ParameterSettings> &parameter_settings) const
{
   _gradf.SetInitialGradient(_context->getMultiGenFcn(), parameter_settings, _grad);
}

////////////////////////////////////////////////////////////////////////////////

#define DEBUG_STREAM(var) << " " #var "=" << var
#include <sys/types.h>
#include <unistd.h>

double RooGradMinimizerFcn::DoEval(const double *x) const
{
      Bool_t parameters_changed = kFALSE;

   // Set the parameter values for this iteration
   for (unsigned index = 0; index < NDim(); index++) {
      // also check whether the function was already evaluated for this set of parameters
      parameters_changed |= SetPdfParamVal(index, x[index]);
   }

//   std::cout << "RooGradMinimizerFcn::DoEval @ PID" << getpid() << ": " DEBUG_STREAM(parameters_changed);

   // Calculate the function for these parameters
   RooAbsReal::setHideOffset(kFALSE);
   double fvalue = _funct->getVal();
   RooAbsReal::setHideOffset(kTRUE);

//   std::cout DEBUG_STREAM(fvalue) << std::endl;

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
      std::cout << "\nprevFCN" << (_funct->isOffsetting() ? "-offset" : "") << " = " << std::setprecision(10) << fvalue
                << std::setprecision(4) << "  ";
      std::cout.flush();
   }

   _evalCounter++;
   //#ifndef NDEBUG
   //  std::cout << "RooGradMinimizerFcn " << this << " evaluations (in DoEval): " << _evalCounter <<
   //  std::endl;
   //#endif
   return fvalue;
}

void RooGradMinimizerFcn::reset_has_been_calculated_flags() const
{
   for (auto it = has_been_calculated.begin(); it != has_been_calculated.end(); ++it) {
      *it = false;
   }
   none_have_been_calculated = true;
}

bool RooGradMinimizerFcn::sync_parameter(double x, std::size_t ix) const
{
   bool parameter_has_changed = (_grad_params[ix] != x);

   if (parameter_has_changed) {
      _grad_params[ix] = x;
      // Set the parameter values for this iteration
      // TODO: this is already done in DoEval as well; find efficient way to do only once
      SetPdfParamVal(ix, x);

      if (!none_have_been_calculated) {
         reset_has_been_calculated_flags();
      }
   }

   return parameter_has_changed;
}

bool RooGradMinimizerFcn::sync_parameters(const double *x) const
{
   bool has_been_synced = false;

   for (std::size_t ix = 0; ix < NDim(); ++ix) {
      bool parameter_has_changed = (_grad_params[ix] != x[ix]);

      if (parameter_has_changed) {
         _grad_params[ix] = x[ix];
         // Set the parameter values for this iteration
         // TODO: this is already done in DoEval as well; find efficient way to do only once
         SetPdfParamVal(ix, x[ix]);
      }

      has_been_synced |= parameter_has_changed;
   }

   if (has_been_synced) {
      reset_has_been_calculated_flags();
   }

   return has_been_synced;
}

void RooGradMinimizerFcn::run_derivator(unsigned int i_component) const
{
   // check whether the derivative was already calculated for this set of parameters
   if (!has_been_calculated[i_component]) {
      // Calculate the derivative etc for these parameters
      _grad[i_component] =
         _gradf.partial_derivative(_context->getMultiGenFcn(), _grad_params.data(),
                                   _context->fitter()->Config().ParamsSettings(), i_component, _grad[i_component]);
      has_been_calculated[i_component] = true;
      none_have_been_calculated = false;
   }
}

double RooGradMinimizerFcn::DoDerivative(const double *x, unsigned int i_component) const
{
   sync_parameters(x);
   run_derivator(i_component);
   return _grad[i_component].derivative;
}

bool RooGradMinimizerFcn::returnsInMinuit2ParameterSpace() const
{
   return true;
}

unsigned int RooGradMinimizerFcn::NDim() const
{
   return get_nDim();
}

////////////////////////////////////////////////////////////////////////////////

void RooGradMinimizerFcn::set_strategy(int istrat)
{
   assert(istrat >= 0);
   ROOT::Minuit2::MnStrategy strategy(static_cast<unsigned int>(istrat));

   set_step_tolerance(strategy.GradientStepTolerance());
   set_grad_tolerance(strategy.GradientTolerance());
   set_ncycles(strategy.GradientNCycles());
}

Bool_t
RooGradMinimizerFcn::Synchronize(std::vector<ROOT::Fit::ParameterSettings> &parameters, Bool_t optConst, Bool_t verbose)
{
   Bool_t returnee = synchronize_parameter_settings(parameters, optConst, verbose);
   synchronize_gradient_parameter_settings(parameters);
   set_strategy(_context->fitter()->Config().MinimizerOptions().Strategy());
   set_error_level(_context->fitter()->Config().MinimizerOptions().ErrorDef());
   return returnee;
}

void RooGradMinimizerFcn::optimizeConstantTerms(bool constStatChange, bool constValChange)
{
   if (constStatChange) {

      RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::CollectErrors);

      oocoutI(static_cast<RooAbsArg *>(nullptr), Eval)
         << "RooGradMinimizerFcn::::synchronize: set of constant parameters changed, rerunning const optimizer"
         << std::endl;
      _funct->constOptimizeTestStatistic(RooAbsArg::ConfigChange, true);
   } else if (constValChange) {
      oocoutI(static_cast<RooAbsArg *>(nullptr), Eval)
         << "RooGradMinimizerFcn::::synchronize: constant parameter values changed, rerunning const optimizer"
         << std::endl;
      _funct->constOptimizeTestStatistic(RooAbsArg::ValueChange, true);
   }

   RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::PrintErrors);
}

void RooGradMinimizerFcn::set_step_tolerance(double step_tolerance) const
{
   _gradf.set_step_tolerance(step_tolerance);
}
void RooGradMinimizerFcn::set_grad_tolerance(double grad_tolerance) const
{
   _gradf.set_grad_tolerance(grad_tolerance);
}
void RooGradMinimizerFcn::set_ncycles(unsigned int ncycles) const
{
   _gradf.set_ncycles(ncycles);
}
void RooGradMinimizerFcn::set_error_level(double error_level) const
{
   _gradf.set_error_level(error_level);
}

std::string RooGradMinimizerFcn::getFunctionName() const
{
   return _funct->GetName();
}

std::string RooGradMinimizerFcn::getFunctionTitle() const
{
   return _funct->GetTitle();
}

void RooGradMinimizerFcn::setOptimizeConst(Int_t flag)
{
   RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::CollectErrors);

   if (_optConst && !flag) {
      if (_context->getPrintLevel() > -1)
         oocoutI(_context, Minimization) << "RooGradMinimizerFcn::setOptimizeConst: deactivating const optimization" << std::endl;
      _funct->constOptimizeTestStatistic(RooAbsArg::DeActivate, true);
      _optConst = flag;
   } else if (!_optConst && flag) {
      if (_context->getPrintLevel() > -1)
         oocoutI(_context, Minimization) << "RooGradMinimizerFcn::setOptimizeConst: activating const optimization" << std::endl;
      _funct->constOptimizeTestStatistic(RooAbsArg::Activate, flag > 1);
      _optConst = flag;
   } else if (_optConst && flag) {
      if (_context->getPrintLevel() > -1)
         oocoutI(_context, Minimization) << "RooGradMinimizerFcn::setOptimizeConst: const optimization already active" << std::endl;
   } else {
      if (_context->getPrintLevel() > -1)
         oocoutI(_context, Minimization) << "RooGradMinimizerFcn::setOptimizeConst: const optimization wasn't active" << std::endl;
   }

   RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::PrintErrors);
}

void RooGradMinimizerFcn::Gradient(const double *x, double *grad) const
{
   ROOT::Math::IMultiGradFunction::Gradient(x, grad);
}

#endif
