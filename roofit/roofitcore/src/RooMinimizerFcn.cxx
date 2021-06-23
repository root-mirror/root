/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 * @(#)root/roofitcore:$Id$
 * Authors:                                                                  *
 *   AL, Alfio Lazzaro,   INFN Milan,        alfio.lazzaro@mi.infn.it        *
 *                                                                           *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

#ifndef __ROOFIT_NOROOMINIMIZER

//////////////////////////////////////////////////////////////////////////////
/// \class RooMinimizerFcn
/// RooMinimizerFcn is an interface to the ROOT::Math::IBaseFunctionMultiDim,
/// a function that ROOT's minimisers use to carry out minimisations.
///

#include "RooMinimizerFcn.h"

#include "RooAbsArg.h"
#include "RooAbsPdf.h"
#include "RooArgSet.h"
#include "RooRealVar.h"
#include "RooAbsRealLValue.h"
#include "RooMsgService.h"
#include "RooMinimizer.h"
#include "RooGaussMinimizer.h"
#include "RooNaNPacker.h"

#include "TClass.h"
#include "TMatrixDSym.h"

#include <fstream>
#include <iomanip>

using namespace std;

RooMinimizerFcn::RooMinimizerFcn(RooAbsReal *funct, RooMinimizer* context,
			   bool verbose) :
  RooAbsMinimizerFcn(*funct->getParameters(RooArgSet()), context, verbose), _funct(funct)
{}



RooMinimizerFcn::RooMinimizerFcn(const RooMinimizerFcn& other) : RooAbsMinimizerFcn(other), ROOT::Math::IBaseFunctionMultiDim(other),
  _funct(other._funct)
{}


RooMinimizerFcn::~RooMinimizerFcn()
{}


ROOT::Math::IBaseFunctionMultiDim* RooMinimizerFcn::Clone() const
{
  return new RooMinimizerFcn(*this) ;
}

void RooMinimizerFcn::setOptimizeConst(Int_t flag)
{
   RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::CollectErrors);

   if (_optConst && !flag) {
      if (_context->getPrintLevel() > -1)
         oocoutI(_context, Minimization) << "RooMinimizerFcn::setOptimizeConst: deactivating const optimization" << endl;
      _funct->constOptimizeTestStatistic(RooAbsArg::DeActivate, true);
      _optConst = flag;
   } else if (!_optConst && flag) {
      if (_context->getPrintLevel() > -1)
         oocoutI(_context, Minimization) << "RooMinimizerFcn::setOptimizeConst: activating const optimization" << endl;
      _funct->constOptimizeTestStatistic(RooAbsArg::Activate, flag > 1);
      _optConst = flag;
   } else if (_optConst && flag) {
      if (_context->getPrintLevel() > -1)
         oocoutI(_context, Minimization) << "RooMinimizerFcn::setOptimizeConst: const optimization already active" << endl;
   } else {
      if (_context->getPrintLevel() > -1)
         oocoutI(_context, Minimization) << "RooMinimizerFcn::setOptimizeConst: const optimization wasn't active" << endl;
   }

   RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::PrintErrors);
}

void RooMinimizerFcn::optimizeConstantTerms(bool constStatChange, bool constValChange) {
   if (constStatChange) {

      RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::CollectErrors) ;

      oocoutI(_context,Minimization) << "RooMinimizerFcn::optimizeConstantTerms: set of constant parameters changed, rerunning const optimizer" << endl ;
      _funct->constOptimizeTestStatistic(RooAbsArg::ConfigChange, true) ;
   } else if (constValChange) {
      oocoutI(_context,Minimization) << "RooMinimizerFcn::optimizeConstantTerms: constant parameter values changed, rerunning const optimizer" << endl ;
      _funct->constOptimizeTestStatistic(RooAbsArg::ValueChange, true) ;
   }

   RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::PrintErrors) ;
}


/// Evaluate function given the parameters in `x`.
double RooMinimizerFcn::DoEval(const double *x) const {

  // Set the parameter values for this iteration
  for (unsigned index = 0; index < _nDim; index++) {
    if (_logfile) (*_logfile) << x[index] << " " ;
    SetPdfParamVal(index,x[index]);
  }

  // Calculate the function for these parameters
  RooAbsReal::setHideOffset(kFALSE) ;
  double fvalue = _funct->getVal();
  RooAbsReal::setHideOffset(kTRUE) ;

  if (!std::isfinite(fvalue) || RooAbsReal::numEvalErrors() > 0 || fvalue > 1e30) {
    printEvalErrors();
    RooAbsReal::clearEvalErrorLog() ;
    _numBadNLL++ ;

    if (_doEvalErrorWall) {
      const double badness = RooNaNPacker::unpackNaN(fvalue);
      fvalue = (std::isfinite(_maxFCN) ? _maxFCN : 0.) + _recoverFromNaNStrength * badness;
    }
  } else {
    if (_evalCounter > 0 && _evalCounter == _numBadNLL) {
      // This is the first time we get a valid function value; while before, the
      // function was always invalid. For invalid  cases, we returned values > 0.
      // Now, we offset valid values such that they are < 0.
      _funcOffset = -fvalue;
    }
    fvalue += _funcOffset;
    _maxFCN = std::max(fvalue, _maxFCN);
  }

  // Optional logging
  if (_logfile)
    (*_logfile) << setprecision(15) << fvalue << setprecision(4) << endl;
  if (_verbose) {
    cout << "\nprevFCN" << (_funct->isOffsetting()?"-offset":"") << " = " << setprecision(10)
         << fvalue << setprecision(4) << "  " ;
    cout.flush() ;
  }

  _evalCounter++ ;

  return fvalue;
}

std::string RooMinimizerFcn::getFunctionName() const
{
   return _funct->GetName();
}

std::string RooMinimizerFcn::getFunctionTitle() const
{
   return _funct->GetTitle();
}

#endif
