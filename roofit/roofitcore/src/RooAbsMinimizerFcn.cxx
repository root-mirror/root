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
//
// RooAbsMinimizerFcn is an interface class to the ROOT::Math function
// for minimization. It contains only the "logistics" of synchronizing
// between Minuit and RooFit. Its subclasses implement actual interfacing
// to Minuit by subclassing IMultiGenFunction or IMultiGradFunction.
//

#include "RooAbsMinimizerFcn.h"

#include "RooAbsArg.h"
#include "RooAbsPdf.h"
#include "RooArgSet.h"
#include "RooRealVar.h"
#include "RooAbsRealLValue.h"
#include "RooMsgService.h"
#include "RooMinimizer.h"
#include "RooNaNPacker.h"

#include "TMatrixDSym.h"

#include <fstream>
#include <iomanip>

using namespace std;

RooAbsMinimizerFcn::RooAbsMinimizerFcn(RooArgList paramList, RooMinimizer *context, bool verbose)
   : _context(context), _verbose(verbose)
{
   // Examine parameter list
   _floatParamList = (RooArgList *)paramList.selectByAttrib("Constant", kFALSE);
   if (_floatParamList->getSize() > 1) {
      _floatParamList->sort();
   }
   _floatParamList->setName("floatParamList");

   _constParamList = (RooArgList *)paramList.selectByAttrib("Constant", kTRUE);
   if (_constParamList->getSize() > 1) {
      _constParamList->sort();
   }
   _constParamList->setName("constParamList");

  // Remove all non-RooRealVar parameters from list (MINUIT cannot handle them)
  for (unsigned int i = 0; i < _floatParamList->size(); ) { // Note: Counting loop, since removing from collection!
    const RooAbsArg* arg = (*_floatParamList).at(i);
    if (!arg->IsA()->InheritsFrom(RooAbsRealLValue::Class())) {
      oocoutW(_context,Minimization) << "RooMinimizerFcn::RooMinimizerFcn: removing parameter "
				     << arg->GetName() << " from list because it is not of type RooRealVar" << endl;
      _floatParamList->remove(*arg);
    } else {
      ++i;
    }
  }

   _nDim = _floatParamList->getSize();

   // Save snapshot of initial lists
   _initFloatParamList = (RooArgList *)_floatParamList->snapshot(kFALSE);
   _initConstParamList = (RooArgList *)_constParamList->snapshot(kFALSE);
}

RooAbsMinimizerFcn::RooAbsMinimizerFcn(const RooAbsMinimizerFcn &other)
   : _context(other._context), _maxFCN(other._maxFCN),
     _funcOffset(other._funcOffset),
     _recoverFromNaNStrength(other._recoverFromNaNStrength),
     _numBadNLL(other._numBadNLL),
     _printEvalErrors(other._printEvalErrors), _evalCounter(other._evalCounter),
     _nDim(other._nDim), _optConst(other._optConst), _logfile(other._logfile), _doEvalErrorWall(other._doEvalErrorWall), _verbose(other._verbose)
{
   _floatParamList = new RooArgList(*other._floatParamList);
   _constParamList = new RooArgList(*other._constParamList);
   _initFloatParamList = (RooArgList *)other._initFloatParamList->snapshot(kFALSE);
   _initConstParamList = (RooArgList *)other._initConstParamList->snapshot(kFALSE);
}

RooAbsMinimizerFcn::~RooAbsMinimizerFcn()
{
   delete _floatParamList;
   delete _initFloatParamList;
   delete _constParamList;
   delete _initConstParamList;
}

/// Internal function to synchronize TMinimizer with current
/// information in RooAbsReal function parameters
Bool_t RooAbsMinimizerFcn::synchronize_parameter_settings(std::vector<ROOT::Fit::ParameterSettings> &parameters, Bool_t optConst, Bool_t verbose)
{
   Bool_t constValChange(kFALSE);
   Bool_t constStatChange(kFALSE);

   Int_t index(0);

   // Handle eventual migrations from constParamList -> floatParamList
   for (index = 0; index < _constParamList->getSize(); index++) {

      RooRealVar *par = dynamic_cast<RooRealVar *>(_constParamList->at(index));
      if (!par)
         continue;

      RooRealVar *oldpar = dynamic_cast<RooRealVar *>(_initConstParamList->at(index));
      if (!oldpar)
         continue;

      // Test if constness changed
      if (!par->isConstant()) {

         // Remove from constList, add to floatList
         _constParamList->remove(*par);
         _floatParamList->add(*par);
         _initFloatParamList->addClone(*oldpar);
         _initConstParamList->remove(*oldpar);
         constStatChange = kTRUE;
         _nDim++;

         if (verbose) {
            oocoutI(_context, Minimization)
               << "RooAbsMinimizerFcn::synchronize: parameter " << par->GetName() << " is now floating." << endl;
         }
      }

      // Test if value changed
      if (par->getVal() != oldpar->getVal()) {
         constValChange = kTRUE;
         if (verbose) {
            oocoutI(_context, Minimization)
               << "RooAbsMinimizerFcn::synchronize: value of constant parameter " << par->GetName() << " changed from "
               << oldpar->getVal() << " to " << par->getVal() << endl;
         }
      }
   }

   // Update reference list
   *_initConstParamList = *_constParamList;

   // Synchronize MINUIT with function state
   // Handle floatParamList
   for (index = 0; index < _floatParamList->getSize(); index++) {
      RooRealVar *par = dynamic_cast<RooRealVar *>(_floatParamList->at(index));

      if (!par)
         continue;

      Double_t pstep(0);
      Double_t pmin(0);
      Double_t pmax(0);

      if (!par->isConstant()) {

         // Verify that floating parameter is indeed of type RooRealVar
         if (!par->IsA()->InheritsFrom(RooRealVar::Class())) {
            oocoutW(_context, Minimization) << "RooAbsMinimizerFcn::fit: Error, non-constant parameter "
                                            << par->GetName() << " is not of type RooRealVar, skipping" << endl;
            _floatParamList->remove(*par);
            index--;
            _nDim--;
            continue;
         }
         // make sure the parameter are in dirty state to enable
         // a real NLL computation when the minimizer calls the function the first time
         // (see issue #7659)
         par->setValueDirty();

         // Set the limits, if not infinite
         if (par->hasMin())
            pmin = par->getMin();
         if (par->hasMax())
            pmax = par->getMax();

         // Calculate step size
         pstep = par->getError();
         if (pstep <= 0) {
            // Floating parameter without error estitimate
            if (par->hasMin() && par->hasMax()) {
               pstep = 0.1 * (pmax - pmin);

               // Trim default choice of error if within 2 sigma of limit
               if (pmax - par->getVal() < 2 * pstep) {
                  pstep = (pmax - par->getVal()) / 2;
               } else if (par->getVal() - pmin < 2 * pstep) {
                  pstep = (par->getVal() - pmin) / 2;
               }

               // If trimming results in zero error, restore default
               if (pstep == 0) {
                  pstep = 0.1 * (pmax - pmin);
               }

            } else {
               pstep = 1;
            }
            if (verbose) {
               oocoutW(_context, Minimization)
                  << "RooAbsMinimizerFcn::synchronize: WARNING: no initial error estimate available for "
                  << par->GetName() << ": using " << pstep << endl;
            }
         }
      } else {
         pmin = par->getVal();
         pmax = par->getVal();
      }

      // new parameter
      if (index >= Int_t(parameters.size())) {

         if (par->hasMin() && par->hasMax()) {
            parameters.emplace_back(par->GetName(), par->getVal(), pstep, pmin, pmax);
         } else {
            parameters.emplace_back(par->GetName(), par->getVal(), pstep);
            if (par->hasMin())
               parameters.back().SetLowerLimit(pmin);
            else if (par->hasMax())
               parameters.back().SetUpperLimit(pmax);
         }

         continue;
      }

      Bool_t oldFixed = parameters[index].IsFixed();
      Double_t oldVar = parameters[index].Value();
      Double_t oldVerr = parameters[index].StepSize();
      Double_t oldVlo = parameters[index].LowerLimit();
      Double_t oldVhi = parameters[index].UpperLimit();

      if (par->isConstant() && !oldFixed) {

         // Parameter changes floating -> constant : update only value if necessary
         if (oldVar != par->getVal()) {
            parameters[index].SetValue(par->getVal());
            if (verbose) {
               oocoutI(_context, Minimization)
                  << "RooAbsMinimizerFcn::synchronize: value of parameter " << par->GetName() << " changed from "
                  << oldVar << " to " << par->getVal() << endl;
            }
         }
         parameters[index].Fix();
         constStatChange = kTRUE;
         if (verbose) {
            oocoutI(_context, Minimization)
               << "RooAbsMinimizerFcn::synchronize: parameter " << par->GetName() << " is now fixed." << endl;
         }

      } else if (par->isConstant() && oldFixed) {

         // Parameter changes constant -> constant : update only value if necessary
         if (oldVar != par->getVal()) {
            parameters[index].SetValue(par->getVal());
            constValChange = kTRUE;

            if (verbose) {
               oocoutI(_context, Minimization)
                  << "RooAbsMinimizerFcn::synchronize: value of fixed parameter " << par->GetName() << " changed from "
                  << oldVar << " to " << par->getVal() << endl;
            }
         }

      } else {
         // Parameter changes constant -> floating
         if (!par->isConstant() && oldFixed) {
            parameters[index].Release();
            constStatChange = kTRUE;

            if (verbose) {
               oocoutI(_context, Minimization)
                  << "RooAbsMinimizerFcn::synchronize: parameter " << par->GetName() << " is now floating." << endl;
            }
         }

         // Parameter changes constant -> floating : update all if necessary
         if (oldVar != par->getVal() || oldVlo != pmin || oldVhi != pmax || oldVerr != pstep) {
            parameters[index].SetValue(par->getVal());
            parameters[index].SetStepSize(pstep);
            if (par->hasMin() && par->hasMax())
               parameters[index].SetLimits(pmin, pmax);
            else if (par->hasMin())
               parameters[index].SetLowerLimit(pmin);
            else if (par->hasMax())
               parameters[index].SetUpperLimit(pmax);
         }

         // Inform user about changes in verbose mode
         if (verbose) {
            // if ierr<0, par was moved from the const list and a message was already printed

            if (oldVar != par->getVal()) {
               oocoutI(_context, Minimization)
                  << "RooAbsMinimizerFcn::synchronize: value of parameter " << par->GetName() << " changed from "
                  << oldVar << " to " << par->getVal() << endl;
            }
            if (oldVlo != pmin || oldVhi != pmax) {
               oocoutI(_context, Minimization)
                  << "RooAbsMinimizerFcn::synchronize: limits of parameter " << par->GetName() << " changed from ["
                  << oldVlo << "," << oldVhi << "] to [" << pmin << "," << pmax << "]" << endl;
            }

            // If oldVerr=0, then parameter was previously fixed
            if (oldVerr != pstep && oldVerr != 0) {
               oocoutI(_context, Minimization)
                  << "RooAbsMinimizerFcn::synchronize: error/step size of parameter " << par->GetName()
                  << " changed from " << oldVerr << " to " << pstep << endl;
            }
         }
      }
   }

   if (optConst) {
      optimizeConstantTerms(constStatChange, constValChange);
   }

   return 0;
}

Bool_t
RooAbsMinimizerFcn::Synchronize(std::vector<ROOT::Fit::ParameterSettings> &parameters, Bool_t optConst, Bool_t verbose) {
   return synchronize_parameter_settings(parameters, optConst, verbose);
}

/// Modify PDF parameter error by ordinal index (needed by MINUIT)
void RooAbsMinimizerFcn::SetPdfParamErr(Int_t index, Double_t value)
{
   static_cast<RooRealVar*>(_floatParamList->at(index))->setError(value);
}

/// Modify PDF parameter error by ordinal index (needed by MINUIT)
void RooAbsMinimizerFcn::ClearPdfParamAsymErr(Int_t index)
{
   static_cast<RooRealVar*>(_floatParamList->at(index))->removeAsymError();
}

/// Modify PDF parameter error by ordinal index (needed by MINUIT)
void RooAbsMinimizerFcn::SetPdfParamErr(Int_t index, Double_t loVal, Double_t hiVal)
{
   static_cast<RooRealVar*>(_floatParamList->at(index))->setAsymError(loVal, hiVal);
}

/// Transfer MINUIT fit results back into RooFit objects.
void RooAbsMinimizerFcn::BackProp(const ROOT::Fit::FitResult &results)
{
   for (unsigned int index = 0; index < _nDim; index++) {
      Double_t value = results.Value(index);
      SetPdfParamVal(index, value);

      // Set the parabolic error
      Double_t err = results.Error(index);
      SetPdfParamErr(index, err);

      Double_t eminus = results.LowerError(index);
      Double_t eplus = results.UpperError(index);

      if (eplus > 0 || eminus < 0) {
         // Store the asymmetric error, if it is available
         SetPdfParamErr(index, eminus, eplus);
      } else {
         // Clear the asymmetric error
         ClearPdfParamAsymErr(index);
      }
   }
}

/// Change the file name for logging of a RooMinimizer of all MINUIT steppings
/// through the parameter space. If inLogfile is null, the current log file
/// is closed and logging is stopped.
Bool_t RooAbsMinimizerFcn::SetLogFile(const char *inLogfile)
{
   if (_logfile) {
      oocoutI(_context, Minimization) << "RooAbsMinimizerFcn::setLogFile: closing previous log file" << endl;
      _logfile->close();
      delete _logfile;
      _logfile = 0;
   }
   _logfile = new ofstream(inLogfile);
   if (!_logfile->good()) {
      oocoutI(_context, Minimization) << "RooAbsMinimizerFcn::setLogFile: cannot open file " << inLogfile << endl;
      _logfile->close();
      delete _logfile;
      _logfile = 0;
   }

   return kFALSE;
}

/// Apply results of given external covariance matrix. i.e. propagate its errors
/// to all RRV parameter representations and give this matrix instead of the
/// HESSE matrix at the next save() call
void RooAbsMinimizerFcn::ApplyCovarianceMatrix(TMatrixDSym &V)
{
   for (unsigned int i = 0; i < _nDim; i++) {
      // Skip fixed parameters
      if (_floatParamList->at(i)->isConstant()) {
         continue;
      }
      SetPdfParamErr(i, sqrt(V(i, i)));
   }
}

/// Set value of parameter i.
Bool_t RooAbsMinimizerFcn::SetPdfParamVal(int index, double value) const
{
  auto par = static_cast<RooRealVar*>(&(*_floatParamList)[index]);

  if (par->getVal()!=value) {
    if (_verbose) cout << par->GetName() << "=" << value << ", " ;
    
    par->setVal(value);
    return kTRUE;
  }

  return kFALSE;
}


/// Print information about why evaluation failed.
/// Using _printEvalErrors, the number of errors printed can be steered.
/// Negative values disable printing.
void RooAbsMinimizerFcn::printEvalErrors() const {
  if (_printEvalErrors < 0)
    return;

  std::ostringstream msg;
  if (_doEvalErrorWall) {
    msg << "RooMinimizerFcn: Minimized function has error status." << endl
        << "Returning maximum FCN so far (" << _maxFCN
        << ") to force MIGRAD to back out of this region. Error log follows.\n";
  } else {
    msg << "RooMinimizerFcn: Minimized function has error status but is ignored.\n";
  }

  msg << "Parameter values: " ;
  for (const auto par : *_floatParamList) {
    auto var = static_cast<const RooRealVar*>(par);
    msg << "\t" << var->GetName() << "=" << var->getVal() ;
  }
  msg << std::endl;

  RooAbsReal::printEvalErrors(msg, _printEvalErrors);
  ooccoutW(_context,Minimization) << msg.str() << endl;
}


////////////////////////////////////////////////////////////////////////////////
/// Logistics

RooArgList *RooAbsMinimizerFcn::GetFloatParamList()
{
   return _floatParamList;
}
RooArgList *RooAbsMinimizerFcn::GetConstParamList()
{
   return _constParamList;
}
RooArgList *RooAbsMinimizerFcn::GetInitFloatParamList()
{
   return _initFloatParamList;
}
RooArgList *RooAbsMinimizerFcn::GetInitConstParamList()
{
   return _initConstParamList;
}

void RooAbsMinimizerFcn::SetEvalErrorWall(Bool_t flag)
{
   _doEvalErrorWall = flag;
}
void RooAbsMinimizerFcn::SetPrintEvalErrors(Int_t numEvalErrors)
{
   _printEvalErrors = numEvalErrors;
}

Double_t &RooAbsMinimizerFcn::GetMaxFCN()
{
   return _maxFCN;
}
Int_t RooAbsMinimizerFcn::GetNumInvalidNLL() const
{
   return _numBadNLL;
}

Int_t RooAbsMinimizerFcn::evalCounter() const
{
   return _evalCounter;
}
void RooAbsMinimizerFcn::zeroEvalCount()
{
   _evalCounter = 0;
}

void RooAbsMinimizerFcn::SetVerbose(Bool_t flag)
{
   _verbose = flag;
}

bool RooAbsMinimizerFcn::getOptConst()
{
   return _optConst;
}


//Double_t RooAbsMinimizerFcn::GetPdfParamVal(Int_t index)
//{
//   // Access PDF parameter value by ordinal index (needed by MINUIT)
//
//   return ((RooRealVar *)_floatParamList->at(index))->getVal();
//}
//
//Double_t RooAbsMinimizerFcn::GetPdfParamErr(Int_t index)
//{
//   // Access PDF parameter error by ordinal index (needed by MINUIT)
//   return ((RooRealVar *)_floatParamList->at(index))->getError();
//}

std::vector<double> RooAbsMinimizerFcn::get_parameter_values() const
{
   // TODO: make a cache for this somewhere so it doesn't have to be recreated on each call
   std::vector<double> values;
   values.reserve(_nDim);

   for (std::size_t index = 0; index < _nDim; ++index) {
      RooRealVar *par = (RooRealVar *)_floatParamList->at(index);
      values.push_back(par->getVal());
   }

   return values;
}

#endif
