/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 * @(#)root/roofitcore:$Id$
 * Authors:                                                                  *
 *   WV, Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu       *
 *   DK, David Kirkby,    UC Irvine,         dkirkby@uci.edu                 *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

/**
\file RooMinuit.cxx
\class RooMinuit
\ingroup Roofitcore

RooMinuit is a wrapper class around TFitter/TMinuit that
provides a seamless interface between the MINUIT functionality
and the native RooFit interface.
RooMinuit can minimize any RooAbsReal function with respect to
its parameters. Usual choices for minimization are RooNLLVar
and RooChi2Var
RooMinuit has methods corresponding to MINUIT functions like
hesse(), migrad(), minos() etc. In each of these function calls
the state of the MINUIT engine is synchronized with the state
of the RooFit variables: any change in variables, change
in the constant status etc is forwarded to MINUIT prior to
execution of the MINUIT call. Afterwards the RooFit objects
are resynchronized with the output state of MINUIT: changes
parameter values, errors are propagated.
Various methods are available to control verbosity, profiling,
automatic PDF optimization.
**/

#include "RooFit.h"

#include "TClass.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include "TH2.h"
#include "TMarker.h"
#include "TGraph.h"
#include "TStopwatch.h"
#include "TFitter.h"
#include "TMinuit.h"
#include "TMatrixDSym.h"
#include "RooMinuit.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooAbsReal.h"
#include "RooAbsRealLValue.h"
#include "RooRealVar.h"
#include "RooFitResult.h"
#include "RooAbsPdf.h"
#include "RooSentinel.h"
#include "RooMsgService.h"
#include "RooPlot.h"


#if (__GNUC__==3&&__GNUC_MINOR__==2&&__GNUC_PATCHLEVEL__==3)
char* operator+( streampos&, char* );
#endif

using namespace std;

ClassImp(RooMinuit);

TVirtualFitter *RooMinuit::_theFitter = nullptr;



////////////////////////////////////////////////////////////////////////////////
/// Cleanup method called by atexit handler installed by RooSentinel
/// to delete all global heap objects when the program is terminated

void RooMinuit::cleanup()
{
  if (_theFitter) {
    delete _theFitter ;
    _theFitter =0 ;
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Construct MINUIT interface to given function. Function can be anything,
/// but is typically a -log(likelihood) implemented by RooNLLVar or a chi^2
/// (implemented by RooChi2Var). Other frequent use cases are a RooAddition
/// of a RooNLLVar plus a penalty or constraint term. This class propagates
/// all RooFit information (floating parameters, their values and errors)
/// to MINUIT before each MINUIT call and propagates all MINUIT information
/// back to the RooFit object at the end of each call (updated parameter
/// values, their (asymmetric errors) etc. The default MINUIT error level
/// for HESSE and MINOS error analysis is taken from the defaultErrorLevel()
/// value of the input function.

RooMinuit::RooMinuit(RooAbsReal& function)
{
  RooSentinel::activate() ;

  // Store function reference
  _evalCounter = 0 ;
  _extV = 0 ;
  _func = &function ;
  _logfile = 0 ;
  _optConst = kFALSE ;
  _verbose = kFALSE ;
  _profile = kFALSE ;
  _handleLocalErrors = kTRUE ;
  _printLevel = 1 ;
  _printEvalErrors = 10 ;
  _warnLevel = -999 ;
  _maxEvalMult = 500 ;
  _doEvalErrorWall = kTRUE ;

  // Examine parameter list
  RooArgSet* paramSet = function.getParameters(RooArgSet()) ;
  RooArgList paramList(*paramSet) ;
  delete paramSet ;

  _floatParamList = (RooArgList*) paramList.selectByAttrib("Constant",kFALSE) ;
  if (_floatParamList->getSize()>1) {
    _floatParamList->sort() ;
  }
  _floatParamList->setName("floatParamList") ;

  _constParamList = (RooArgList*) paramList.selectByAttrib("Constant",kTRUE) ;
  if (_constParamList->getSize()>1) {
    _constParamList->sort() ;
  }
  _constParamList->setName("constParamList") ;

  // Remove all non-RooRealVar parameters from list (MINUIT cannot handle them)
  TIterator* pIter = _floatParamList->createIterator() ;
  RooAbsArg* arg ;
  while((arg=(RooAbsArg*)pIter->Next())) {
    if (!arg->IsA()->InheritsFrom(RooAbsRealLValue::Class())) {
      coutW(Minimization) << "RooMinuit::RooMinuit: removing parameter " << arg->GetName()
			  << " from list because it is not of type RooRealVar" << endl ;
      _floatParamList->remove(*arg) ;
    }
  }
  _nPar      = _floatParamList->getSize() ;
  delete pIter ;

  updateFloatVec() ;

  // Save snapshot of initial lists
  _initFloatParamList = (RooArgList*) _floatParamList->snapshot(kFALSE) ;
  _initConstParamList = (RooArgList*) _constParamList->snapshot(kFALSE) ;

  // Initialize MINUIT
  Int_t nPar= _floatParamList->getSize() + _constParamList->getSize() ;
  if (_theFitter) delete _theFitter ;
  _theFitter = new TFitter(nPar*2+1) ; //WVE Kludge, nPar*2 works around TMinuit memory allocation bug
  _theFitter->SetObjectFit(this) ;

  // Shut up for now
  setPrintLevel(-1) ;
  _theFitter->Clear();

  // Tell MINUIT to use our global glue function
  _theFitter->SetFCN(RooMinuitGlue);

  // Use +0.5 for 1-sigma errors
  setErrorLevel(function.defaultErrorLevel()) ;

  // Declare our parameters to MINUIT
  synchronize(kFALSE) ;

  // Reset the *largest* negative log-likelihood value we have seen so far
  _maxFCN= -1e30 ;
  _numBadNLL = 0 ;

  // Now set default verbosity
  if (RooMsgService::instance().silentMode()) {
    setWarnLevel(-1) ;
    setPrintLevel(-1) ;
  } else {
    setWarnLevel(1) ;
    setPrintLevel(1) ;
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Destructor

RooMinuit::~RooMinuit()
{
  delete _floatParamList ;
  delete _initFloatParamList ;
  delete _constParamList ;
  delete _initConstParamList ;
  if (_extV) {
    delete _extV ;
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Change MINUIT strategy to istrat. Accepted codes
/// are 0,1,2 and represent MINUIT strategies for dealing
/// most efficiently with fast FCNs (0), expensive FCNs (2)
/// and 'intermediate' FCNs (1)

void RooMinuit::setStrategy(Int_t istrat)
{
  Double_t stratArg(istrat) ;
  _theFitter->ExecuteCommand("SET STR",&stratArg,1) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Set the level for MINUIT error analysis to the given
/// value. This function overrides the default value
/// that is taken in the RooMinuit constructor from
/// the defaultErrorLevel() method of the input function

void RooMinuit::setErrorLevel(Double_t level)
{
  _theFitter->ExecuteCommand("SET ERR",&level,1);
}



////////////////////////////////////////////////////////////////////////////////
/// Change MINUIT epsilon

void RooMinuit::setEps(Double_t eps)
{
  _theFitter->ExecuteCommand("SET EPS",&eps,1) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Enable internal likelihood offsetting for enhanced numeric precision

void RooMinuit::setOffsetting(Bool_t flag) 
{ 
  _func->enableOffsetting(flag) ; 
}


////////////////////////////////////////////////////////////////////////////////
/// Parse traditional RooAbsPdf::fitTo driver options
///
///  s - Run Hesse first to estimate initial step size
///  m - Run Migrad only
///  h - Run Hesse to estimate errors
///  v - Verbose mode
///  l - Log parameters after each Minuit steps to file
///  t - Activate profile timer
///  r - Save fit result
///  0 - Run Migrad with strategy 0

RooFitResult* RooMinuit::fit(const char* options)
{
  if (_floatParamList->getSize()==0) {
    return 0 ;
  }

  _theFitter->SetObjectFit(this) ;

  TString opts(options) ;
  opts.ToLower() ;

  // Initial configuration
  if (opts.Contains("v")) setVerbose(1) ;
  if (opts.Contains("t")) setProfile(1) ;
  if (opts.Contains("l")) setLogFile(Form("%s.log",_func->GetName())) ;
  if (opts.Contains("c")) optimizeConst(1) ;

  // Fitting steps
  if (opts.Contains("s")) hesse() ;
  if (opts.Contains("0")) setStrategy(0) ;
  migrad() ;
  if (opts.Contains("0")) setStrategy(1) ;
  if (opts.Contains("h")||!opts.Contains("m")) hesse() ;
  if (!opts.Contains("m")) minos() ;

  return (opts.Contains("r")) ? save() : 0 ;
}



////////////////////////////////////////////////////////////////////////////////
/// Execute MIGRAD. Changes in parameter values
/// and calculated errors are automatically
/// propagated back the RooRealVars representing
/// the floating parameters in the MINUIT operation

Int_t RooMinuit::migrad()
{
  if (_floatParamList->getSize()==0) {
    return -1 ;
  }

  _theFitter->SetObjectFit(this) ;

  Double_t arglist[2];
  arglist[0]= _maxEvalMult*_nPar; // maximum iterations
  arglist[1]= 1.0;       // tolerance

  synchronize(_verbose) ;
  profileStart() ;
  RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::CollectErrors) ;
  RooAbsReal::clearEvalErrorLog() ;
  _status= _theFitter->ExecuteCommand("MIGRAD",arglist,2);
  RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::PrintErrors) ;
  profileStop() ;
  backProp() ;

  saveStatus("MIGRAD",_status) ;

  return _status ;
}



////////////////////////////////////////////////////////////////////////////////
/// Execute HESSE. Changes in parameter values
/// and calculated errors are automatically
/// propagated back the RooRealVars representing
/// the floating parameters in the MINUIT operation

Int_t RooMinuit::hesse()
{
  if (_floatParamList->getSize()==0) {
    return -1 ;
  }

  _theFitter->SetObjectFit(this) ;

  Double_t arglist[2];
  arglist[0]= _maxEvalMult*_nPar; // maximum iterations

  synchronize(_verbose) ;
  profileStart() ;
  RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::CollectErrors) ;
  RooAbsReal::clearEvalErrorLog() ;
  _status= _theFitter->ExecuteCommand("HESSE",arglist,1);
  RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::PrintErrors) ;
  profileStop() ;
  backProp() ;

  saveStatus("HESSE",_status) ;

  return _status ;
}



////////////////////////////////////////////////////////////////////////////////
/// Execute MINOS. Changes in parameter values
/// and calculated errors are automatically
/// propagated back the RooRealVars representing
/// the floating parameters in the MINUIT operation

Int_t RooMinuit::minos()
{
  if (_floatParamList->getSize()==0) {
    return -1 ;
  }

  _theFitter->SetObjectFit(this) ;

  Double_t arglist[2];
  arglist[0]= _maxEvalMult*_nPar; // maximum iterations

  synchronize(_verbose) ;
  profileStart() ;
  RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::CollectErrors) ;
  RooAbsReal::clearEvalErrorLog() ;
  _status= _theFitter->ExecuteCommand("MINOS",arglist,1);
  // check also the status of Minos looking at fCstatu
  if (_status == 0 && gMinuit->fCstatu != "SUCCESSFUL") {
    if (gMinuit->fCstatu == "FAILURE" || 
	gMinuit->fCstatu == "PROBLEMS") _status = 5; 
    _status = 6; 
  } 
  
  RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::PrintErrors) ;
  profileStop() ;
  backProp() ;

  saveStatus("MINOS",_status) ;
  return _status ;
}


// added FMV, 08/18/03

////////////////////////////////////////////////////////////////////////////////
/// Execute MINOS for given list of parameters. Changes in parameter values
/// and calculated errors are automatically
/// propagated back the RooRealVars representing
/// the floating parameters in the MINUIT operation

Int_t RooMinuit::minos(const RooArgSet& minosParamList)
{
  if (_floatParamList->getSize()==0) {
    return -1 ;
  }

  _theFitter->SetObjectFit(this) ;

  Int_t nMinosPar(0) ;
  Double_t* arglist = new Double_t[_nPar+1];

  if (minosParamList.getSize()>0) {
    TIterator* aIter = minosParamList.createIterator() ;
    RooAbsArg* arg ;
    while((arg=(RooAbsArg*)aIter->Next())) {
      RooAbsArg* par = _floatParamList->find(arg->GetName());
      if (par && !par->isConstant()) {
	Int_t index = _floatParamList->index(par);
	nMinosPar++;
        arglist[nMinosPar]=index+1;
      }
    }
    delete aIter ;
  }
  arglist[0]= _maxEvalMult*_nPar; // maximum iterations

  synchronize(_verbose) ;
  profileStart() ;
  RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::CollectErrors) ;
  RooAbsReal::clearEvalErrorLog() ;
  _status= _theFitter->ExecuteCommand("MINOS",arglist,1+nMinosPar);
  // check also the status of Minos looking at fCstatu
  if (_status == 0 && gMinuit->fCstatu != "SUCCESSFUL") {
    if (gMinuit->fCstatu == "FAILURE" || 
	gMinuit->fCstatu == "PROBLEMS") _status = 5; 
    _status = 6; 
  } 
  RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::PrintErrors) ;
  profileStop() ;
  backProp() ;

  delete[] arglist ;
  
  saveStatus("MINOS",_status) ;

  return _status ;
}



////////////////////////////////////////////////////////////////////////////////
/// Execute SEEK. Changes in parameter values
/// and calculated errors are automatically
/// propagated back the RooRealVars representing
/// the floating parameters in the MINUIT operation

Int_t RooMinuit::seek()
{
  if (_floatParamList->getSize()==0) {
    return -1 ;
  }

  _theFitter->SetObjectFit(this) ;

  Double_t arglist[2];
  arglist[0]= _maxEvalMult*_nPar; // maximum iterations

  synchronize(_verbose) ;
  profileStart() ;
  RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::CollectErrors) ;
  RooAbsReal::clearEvalErrorLog() ;
  _status= _theFitter->ExecuteCommand("SEEK",arglist,1);
  RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::PrintErrors) ;
  profileStop() ;
  backProp() ;

  saveStatus("SEEK",_status) ;
  
  return _status ;
}



////////////////////////////////////////////////////////////////////////////////
/// Execute SIMPLEX. Changes in parameter values
/// and calculated errors are automatically
/// propagated back the RooRealVars representing
/// the floating parameters in the MINUIT operation

Int_t RooMinuit::simplex()
{
  if (_floatParamList->getSize()==0) {
    return -1 ;
  }

  _theFitter->SetObjectFit(this) ;

  Double_t arglist[2];
  arglist[0]= _maxEvalMult*_nPar; // maximum iterations
  arglist[1]= 1.0;       // tolerance

  synchronize(_verbose) ;
  profileStart() ;
  RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::CollectErrors) ;
  RooAbsReal::clearEvalErrorLog() ;
  _status= _theFitter->ExecuteCommand("SIMPLEX",arglist,2);
  RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::PrintErrors) ;
  profileStop() ;
  backProp() ;

  saveStatus("SIMPLEX",_status) ;
  
  return _status ;
}



////////////////////////////////////////////////////////////////////////////////
/// Execute IMPROVE. Changes in parameter values
/// and calculated errors are automatically
/// propagated back the RooRealVars representing
/// the floating parameters in the MINUIT operation

Int_t RooMinuit::improve()
{
  if (_floatParamList->getSize()==0) {
    return -1 ;
  }

  _theFitter->SetObjectFit(this) ;

  Double_t arglist[2];
  arglist[0]= _maxEvalMult*_nPar; // maximum iterations

  synchronize(_verbose) ;
  profileStart() ;
  RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::CollectErrors) ;
  RooAbsReal::clearEvalErrorLog() ;
  _status= _theFitter->ExecuteCommand("IMPROVE",arglist,1);
  RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::PrintErrors) ;
  profileStop() ;
  backProp() ;

  saveStatus("IMPROVE",_status) ;
  
  return _status ;
}



////////////////////////////////////////////////////////////////////////////////
/// Change the MINUIT internal printing level

Int_t RooMinuit::setPrintLevel(Int_t newLevel)
{
  Int_t ret = _printLevel ;
  Double_t arg(newLevel) ;
  _theFitter->ExecuteCommand("SET PRINT",&arg,1);
  _printLevel = newLevel ;
  return ret ;
}



////////////////////////////////////////////////////////////////////////////////
/// Instruct MINUIT to suppress warnings

void RooMinuit::setNoWarn()
{
  Double_t arg(0) ;
  _theFitter->ExecuteCommand("SET NOWARNINGS",&arg,1);
  _warnLevel = -1 ;
}



////////////////////////////////////////////////////////////////////////////////
/// Set MINUIT warning level to given level

Int_t RooMinuit::setWarnLevel(Int_t newLevel)
{
  if (newLevel==_warnLevel) {
    return _warnLevel ;
  }

  Int_t ret = _warnLevel ;
  Double_t arg(newLevel) ;

  if (newLevel>=0) {
    _theFitter->ExecuteCommand("SET WARNINGS",&arg,1);
  } else {
    Double_t arg2(0) ;
    _theFitter->ExecuteCommand("SET NOWARNINGS",&arg2,1);
  }
  _warnLevel = newLevel ;

  return ret ;
}



////////////////////////////////////////////////////////////////////////////////
/// Internal function to synchronize TMinuit with current
/// information in RooAbsReal function parameters

Bool_t RooMinuit::synchronize(Bool_t verbose)
{
  Int_t oldPrint = setPrintLevel(-1) ;
  gMinuit->fNwrmes[0] = 0;  // to clear buffer
  Int_t oldWarn = setWarnLevel(-1) ;

  Bool_t constValChange(kFALSE) ;
  Bool_t constStatChange(kFALSE) ;

  Int_t index(0) ;

  // Handle eventual migrations from constParamList -> floatParamList
  for(index= 0; index < _constParamList->getSize() ; index++) {
    RooRealVar *par= dynamic_cast<RooRealVar*>(_constParamList->at(index)) ;
    if (!par) continue ;

    RooRealVar *oldpar= dynamic_cast<RooRealVar*>(_initConstParamList->at(index)) ;
    if (!oldpar) continue ;

    // Test if constness changed
    if (!par->isConstant()) {

      // Remove from constList, add to floatList
      _constParamList->remove(*par) ;
      _floatParamList->add(*par) ;
      _initFloatParamList->addClone(*oldpar) ;
      _initConstParamList->remove(*oldpar) ;
      constStatChange=kTRUE ;
      _nPar++ ;

      if (verbose) {
	coutI(Minimization) << "RooMinuit::synchronize: parameter " << par->GetName() << " is now floating." << endl ;
      }
    }

    // Test if value changed
    if (par->getVal()!= oldpar->getVal()) {
      constValChange=kTRUE ;
      if (verbose) {
	coutI(Minimization) << "RooMinuit::synchronize: value of constant parameter " << par->GetName()
			    << " changed from " << oldpar->getVal() << " to " << par->getVal() << endl ;
      }
    }

  }

  // Update reference list
  *_initConstParamList = *_constParamList ;


  // Synchronize MINUIT with function state
  for(index= 0; index < _nPar; index++) {
    RooRealVar *par= dynamic_cast<RooRealVar*>(_floatParamList->at(index)) ;
    if (!par) continue ;

    Double_t pstep(0) ;
    Double_t pmin(0) ;
    Double_t pmax(0) ;
    
    if(!par->isConstant()) {

      // Verify that floating parameter is indeed of type RooRealVar
      if (!par->IsA()->InheritsFrom(RooRealVar::Class())) {
	coutW(Minimization) << "RooMinuit::fit: Error, non-constant parameter " << par->GetName()
			    << " is not of type RooRealVar, skipping" << endl ;
	continue ;
      }

      // Set the limits, if not infinite
      if (par->hasMin() && par->hasMax()) {
	pmin = par->getMin();
	pmax = par->getMax();
      }

      // Calculate step size
      pstep= par->getError();
      if(pstep <= 0) {
	// Floating parameter without error estitimate
	if (par->hasMin() && par->hasMax()) {
	  pstep= 0.1*(pmax-pmin);

	  // Trim default choice of error if within 2 sigma of limit
	  if (pmax - par->getVal() < 2*pstep) {
	    pstep = (pmax - par->getVal())/2 ;
	  } else if (par->getVal() - pmin < 2*pstep) {
	    pstep = (par->getVal() - pmin )/2 ;
	  }

	  // If trimming results in zero error, restore default
	  if (pstep==0) {
	    pstep= 0.1*(pmax-pmin);
	  }

	} else {
	  pstep=1 ;
	}
	if(_verbose) {
	  coutW(Minimization) << "RooMinuit::synchronize: WARNING: no initial error estimate available for "
			      << par->GetName() << ": using " << pstep << endl;
	}
      }
    } else {
      pmin = par->getVal() ;
      pmax = par->getVal() ;
    }

    // Extract previous information
    Double_t oldVar,oldVerr,oldVlo,oldVhi ;
    char oldParname[100] ;
    Int_t ierr = _theFitter->GetParameter(index,oldParname,oldVar,oldVerr,oldVlo,oldVhi)  ;

    // Determine if parameters is currently fixed in MINUIT

    Int_t ix ;
    Bool_t oldFixed(kFALSE) ;
    if (ierr>=0) {
      for (ix = 1; ix <= gMinuit->fNpfix; ++ix) {
        if (gMinuit->fIpfix[ix-1] == index+1) oldFixed=kTRUE ;
      }
    }

    if (par->isConstant() && !oldFixed) {

      // Parameter changes floating -> constant : update only value if necessary
      if (oldVar!=par->getVal()) {
	Double_t arglist[2] ;
	arglist[0] = index+1 ;
	arglist[1] = par->getVal() ;
	_theFitter->ExecuteCommand("SET PAR",arglist,2) ;
	if (verbose) {
	  coutI(Minimization) << "RooMinuit::synchronize: value of parameter " << par->GetName() << " changed from " << oldVar << " to " << par->getVal() << endl ;
	}
      }

      _theFitter->FixParameter(index) ;
      constStatChange=kTRUE ;
      if (verbose) {
	coutI(Minimization) << "RooMinuit::synchronize: parameter " << par->GetName() << " is now fixed." << endl ;
      }

    } else if (par->isConstant() && oldFixed) {

      // Parameter changes constant -> constant : update only value if necessary
      if (oldVar!=par->getVal()) {
	Double_t arglist[2] ;
	arglist[0] = index+1 ;
	arglist[1] = par->getVal() ;
	_theFitter->ExecuteCommand("SET PAR",arglist,2) ;
	constValChange=kTRUE ;

	if (verbose) {
	  coutI(Minimization) << "RooMinuit::synchronize: value of fixed parameter " << par->GetName() << " changed from " << oldVar << " to " << par->getVal() << endl ;
	}
      }

    } else {

      if (!par->isConstant() && oldFixed) {
	_theFitter->ReleaseParameter(index) ;
	constStatChange=kTRUE ;

	if (verbose) {
	  coutI(Minimization) << "RooMinuit::synchronize: parameter " << par->GetName() << " is now floating." << endl ;
	}
      }

      // Parameter changes constant -> floating : update all if necessary
      if (oldVar!=par->getVal() || oldVlo!=pmin || oldVhi != pmax || oldVerr!=pstep) {
	_theFitter->SetParameter(index, par->GetName(), par->getVal(), pstep, pmin, pmax);
      }

      // Inform user about changes in verbose mode
      if (verbose && ierr>=0) {
	// if ierr<0, par was moved from the const list and a message was already printed

	if (oldVar!=par->getVal()) {
	  coutI(Minimization) << "RooMinuit::synchronize: value of parameter " << par->GetName() << " changed from " << oldVar << " to " << par->getVal() << endl ;
	}
	if (oldVlo!=pmin || oldVhi!=pmax) {
	  coutI(Minimization) << "RooMinuit::synchronize: limits of parameter " << par->GetName() << " changed from [" << oldVlo << "," << oldVhi
	       << "] to [" << pmin << "," << pmax << "]" << endl ;
	}

	// If oldVerr=0, then parameter was previously fixed
	if (oldVerr!=pstep && oldVerr!=0) {
	coutI(Minimization) << "RooMinuit::synchronize: error/step size of parameter " << par->GetName() << " changed from " << oldVerr << " to " << pstep << endl ;
	}
      }
    }
  }


  gMinuit->fNwrmes[0] = 0;  // to clear buffer
  oldWarn = setWarnLevel(oldWarn) ;
  oldPrint = setPrintLevel(oldPrint) ;

  if (_optConst) {
    if (constStatChange) {

      RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::CollectErrors) ;

      coutI(Minimization) << "RooMinuit::synchronize: set of constant parameters changed, rerunning const optimizer" << endl ;
      _func->constOptimizeTestStatistic(RooAbsArg::ConfigChange) ;
    } else if (constValChange) {
      coutI(Minimization) << "RooMinuit::synchronize: constant parameter values changed, rerunning const optimizer" << endl ;
      _func->constOptimizeTestStatistic(RooAbsArg::ValueChange) ;
    }

    RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::PrintErrors) ;

  }

  updateFloatVec() ;

  return 0 ;
}




////////////////////////////////////////////////////////////////////////////////
/// If flag is true, perform constant term optimization on
/// function being minimized.

void RooMinuit::optimizeConst(Int_t flag)
{
  RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::CollectErrors) ;

  if (_optConst && !flag){
    if (_printLevel>-1) coutI(Minimization) << "RooMinuit::optimizeConst: deactivating const optimization" << endl ;
    _func->constOptimizeTestStatistic(RooAbsArg::DeActivate,flag>1) ;
    _optConst = flag ;
  } else if (!_optConst && flag) {
    if (_printLevel>-1) coutI(Minimization) << "RooMinuit::optimizeConst: activating const optimization" << endl ;
    _func->constOptimizeTestStatistic(RooAbsArg::Activate,flag>1) ;
    _optConst = flag ;
  } else if (_optConst && flag) {
    if (_printLevel>-1) coutI(Minimization) << "RooMinuit::optimizeConst: const optimization already active" << endl ;
  } else {
    if (_printLevel>-1) coutI(Minimization) << "RooMinuit::optimizeConst: const optimization wasn't active" << endl ;
  }

  RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::PrintErrors) ;

}



////////////////////////////////////////////////////////////////////////////////
/// Save and return a RooFitResult snaphot of current minimizer status.
/// This snapshot contains the values of all constant parameters,
/// the value of all floating parameters at RooMinuit construction and
/// after the last MINUIT operation, the MINUIT status, variance quality,
/// EDM setting, number of calls with evaluation problems, the minimized
/// function value and the full correlation matrix

RooFitResult* RooMinuit::save(const char* userName, const char* userTitle)
{
  TString name,title ;
  name = userName ? userName : Form("%s", _func->GetName()) ;
  title = userTitle ? userTitle : Form("%s", _func->GetTitle()) ;

  if (_floatParamList->getSize()==0) {
    RooFitResult* fitRes = new RooFitResult(name,title) ;
    fitRes->setConstParList(*_constParamList) ;
    fitRes->setInitParList(RooArgList()) ;
    fitRes->setFinalParList(RooArgList()) ;
    fitRes->setStatus(-999) ;
    fitRes->setCovQual(-999) ;
    fitRes->setMinNLL(_func->getVal()) ;
    fitRes->setNumInvalidNLL(0) ;
    fitRes->setEDM(-999) ;
    return fitRes ;
  }

  RooFitResult* fitRes = new RooFitResult(name,title) ;

  // Move eventual fixed parameters in floatList to constList
  Int_t i ;
  RooArgList saveConstList(*_constParamList) ;
  RooArgList saveFloatInitList(*_initFloatParamList) ;
  RooArgList saveFloatFinalList(*_floatParamList) ;
  for (i=0 ; i<_floatParamList->getSize() ; i++) {
    RooAbsArg* par = _floatParamList->at(i) ;
    if (par->isConstant()) {
      saveFloatInitList.remove(*saveFloatInitList.find(par->GetName()),kTRUE) ;
      saveFloatFinalList.remove(*par) ;
      saveConstList.add(*par) ;
    }
  }
  saveConstList.sort() ;

  fitRes->setConstParList(saveConstList) ;
  fitRes->setInitParList(saveFloatInitList) ;

  Double_t edm, errdef, minVal;
  Int_t nvpar, nparx;
  Int_t icode = _theFitter->GetStats(minVal, edm, errdef, nvpar, nparx);
  fitRes->setStatus(_status) ;
  fitRes->setCovQual(icode) ;
  fitRes->setMinNLL(minVal) ;
  fitRes->setNumInvalidNLL(_numBadNLL) ;
  fitRes->setEDM(edm) ;
  fitRes->setFinalParList(saveFloatFinalList) ;
  if (!_extV) {
    fitRes->fillCorrMatrix() ;
  } else {
    fitRes->setCovarianceMatrix(*_extV) ;
  }

  fitRes->setStatusHistory(_statusHistory) ;

  return fitRes ;
}




////////////////////////////////////////////////////////////////////////////////
/// Create and draw a TH2 with the error contours in parameters var1 and v2 at up to 6 'sigma' settings
/// where 'sigma' is calculated as n*n*errorLevel

RooPlot* RooMinuit::contour(RooRealVar& var1, RooRealVar& var2, Double_t n1, Double_t n2, Double_t n3, Double_t n4, Double_t n5, Double_t n6)
{

  _theFitter->SetObjectFit(this) ;

  RooArgList* paramSave = (RooArgList*) _floatParamList->snapshot() ;

  // Verify that both variables are floating parameters of PDF
  Int_t index1= _floatParamList->index(&var1);
  if(index1 < 0) {
    coutE(Minimization) << "RooMinuit::contour(" << GetName()
			<< ") ERROR: " << var1.GetName() << " is not a floating parameter of " << _func->GetName() << endl ;
    return 0;
  }

  Int_t index2= _floatParamList->index(&var2);
  if(index2 < 0) {
    coutE(Minimization) << "RooMinuit::contour(" << GetName()
			<< ") ERROR: " << var2.GetName() << " is not a floating parameter of PDF " << _func->GetName() << endl ;
    return 0;
  }

  // create and draw a frame
  RooPlot* frame = new RooPlot(var1,var2) ;

  // draw a point at the current parameter values
  TMarker *point= new TMarker(var1.getVal(), var2.getVal(), 8);
  frame->addObject(point) ;

  // remember our original value of ERRDEF
  Double_t errdef= gMinuit->fUp;

  Double_t n[6] ;
  n[0] = n1 ; n[1] = n2 ; n[2] = n3 ; n[3] = n4 ; n[4] = n5 ; n[5] = n6 ;


  for (Int_t ic = 0 ; ic<6 ; ic++) {
    if(n[ic] > 0) {
      // set the value corresponding to an n1-sigma contour
      gMinuit->SetErrorDef(n[ic]*n[ic]*errdef);
      // calculate and draw the contour
      TGraph* graph= (TGraph*)gMinuit->Contour(50, index1, index2);
      if (!graph) {
	coutE(Minimization) << "RooMinuit::contour(" << GetName() << ") ERROR: MINUIT did not return a contour graph for n=" << n[ic] << endl ;
      } else {
	graph->SetName(Form("contour_%s_n%f",_func->GetName(),n[ic])) ;
	graph->SetLineStyle(ic+1) ;
	graph->SetLineWidth(2) ;
	graph->SetLineColor(kBlue) ;
	frame->addObject(graph,"L") ;
      }
    }
  }

  // restore the original ERRDEF
  gMinuit->SetErrorDef(errdef);

  // restore parameter values
  *_floatParamList = *paramSave ;
  delete paramSave ;


  return frame ;
}



////////////////////////////////////////////////////////////////////////////////
/// Change the file name for logging of a RooMinuit of all MINUIT steppings
/// through the parameter space. If inLogfile is null, the current log file
/// is closed and logging is stopped.

Bool_t RooMinuit::setLogFile(const char* inLogfile)
{
  if (_logfile) {
    coutI(Minimization) << "RooMinuit::setLogFile: closing previous log file" << endl ;
    _logfile->close() ;
    delete _logfile ;
    _logfile = 0 ;
  }
  _logfile = new ofstream(inLogfile) ;
  if (!_logfile->good()) {
    coutI(Minimization) << "RooMinuit::setLogFile: cannot open file " << inLogfile << endl ;
    _logfile->close() ;
    delete _logfile ;
    _logfile= 0;
  }
  return kFALSE ;
}



////////////////////////////////////////////////////////////////////////////////
/// Access PDF parameter value by ordinal index (needed by MINUIT)

Double_t RooMinuit::getPdfParamVal(Int_t index)
{
  return ((RooRealVar*)_floatParamList->at(index))->getVal() ;
}



////////////////////////////////////////////////////////////////////////////////
/// Access PDF parameter error by ordinal index (needed by MINUIT)

Double_t RooMinuit::getPdfParamErr(Int_t index)
{
  return ((RooRealVar*)_floatParamList->at(index))->getError() ;
}



////////////////////////////////////////////////////////////////////////////////
/// Modify PDF parameter value by ordinal index (needed by MINUIT)

Bool_t RooMinuit::setPdfParamVal(Int_t index, Double_t value, Bool_t verbose)
{
  //RooRealVar* par = (RooRealVar*)_floatParamList->at(index) ;
  RooRealVar* par = (RooRealVar*)_floatParamVec[index] ;

  if (par->getVal()!=value) {
    if (verbose) cout << par->GetName() << "=" << value << ", " ;
    par->setVal(value) ;
    return kTRUE ;
  }

  return kFALSE ;
}



////////////////////////////////////////////////////////////////////////////////
/// Modify PDF parameter error by ordinal index (needed by MINUIT)

void RooMinuit::setPdfParamErr(Int_t index, Double_t value)
{
  ((RooRealVar*)_floatParamList->at(index))->setError(value) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Modify PDF parameter error by ordinal index (needed by MINUIT)

void RooMinuit::clearPdfParamAsymErr(Int_t index)
{
  ((RooRealVar*)_floatParamList->at(index))->removeAsymError() ;
}


////////////////////////////////////////////////////////////////////////////////
/// Modify PDF parameter error by ordinal index (needed by MINUIT)

void RooMinuit::setPdfParamErr(Int_t index, Double_t loVal, Double_t hiVal)
{
  ((RooRealVar*)_floatParamList->at(index))->setAsymError(loVal,hiVal) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Start profiling timer

void RooMinuit::profileStart()
{
  if (_profile) {
    _timer.Start() ;
    _cumulTimer.Start(kFALSE) ;
  }
}




////////////////////////////////////////////////////////////////////////////////
/// Stop profiling timer and report results of last session

void RooMinuit::profileStop()
{
  if (_profile) {
    _timer.Stop() ;
    _cumulTimer.Stop() ;
    coutI(Minimization) << "Command timer: " ; _timer.Print() ;
    coutI(Minimization) << "Session timer: " ; _cumulTimer.Print() ;
  }
}





////////////////////////////////////////////////////////////////////////////////
/// Transfer MINUIT fit results back into RooFit objects

void RooMinuit::backProp()
{
  Double_t val,err,vlo,vhi, eplus, eminus, eparab, globcc;
  char buffer[64000];
  Int_t index ;
  for(index= 0; index < _nPar; index++) {
    _theFitter->GetParameter(index, buffer, val, err, vlo, vhi);
    setPdfParamVal(index, val);
    _theFitter->GetErrors(index, eplus, eminus, eparab, globcc);

    // Set the parabolic error
    setPdfParamErr(index, err);

    if(eplus > 0 || eminus < 0) {
      // Store the asymmetric error, if it is available
      setPdfParamErr(index, eminus,eplus);
    } else {
      // Clear the asymmetric error
      clearPdfParamAsymErr(index) ;
    }
  }
}


////////////////////////////////////////////////////////////////////////////////

void RooMinuit::updateFloatVec() 
{
  _floatParamVec.clear() ;
  RooFIter iter = _floatParamList->fwdIterator() ;
  RooAbsArg* arg ;
  _floatParamVec.resize(_floatParamList->getSize()) ;
  Int_t i(0) ;
  while((arg=iter.next())) {
    _floatParamVec[i++] = arg ;
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Apply results of given external covariance matrix. i.e. propagate its errors
/// to all RRV parameter representations and give this matrix instead of the
/// HESSE matrix at the next save() call

void RooMinuit::applyCovarianceMatrix(TMatrixDSym& V)
{
  _extV = (TMatrixDSym*) V.Clone() ;

  for (Int_t i=0 ; i<getNPar() ; i++) {
    // Skip fixed parameters
    if (_floatParamList->at(i)->isConstant()) {
      continue ;
    }
    RooMinuit* context = (RooMinuit*) RooMinuit::_theFitter->GetObjectFit() ;
    if (context && context->_verbose)
       cout << "setting parameter " << i << " error to " << sqrt((*_extV)(i,i)) << endl ;
    setPdfParamErr(i, sqrt((*_extV)(i,i))) ;
  }

}




void RooMinuitGlue(Int_t& /*np*/, Double_t* /*gin*/,
		   Double_t &f, Double_t *par, Int_t /*flag*/)
{
  // Static function that interfaces minuit with RooMinuit

  // Retrieve fit context and its components
  RooMinuit* context = (RooMinuit*) RooMinuit::_theFitter->GetObjectFit() ;
  ofstream* logf   = context->logfile() ;
  Double_t& maxFCN = context->maxFCN() ;
  Bool_t verbose   = context->_verbose ;

  // Set the parameter values for this iteration
  Int_t nPar= context->getNPar();
  for(Int_t index= 0; index < nPar; index++) {
    if (logf) (*logf) << par[index] << " " ;
    context->setPdfParamVal(index, par[index],verbose);
  }

  // Calculate the function for these parameters
  RooAbsReal::setHideOffset(kFALSE) ;
  f= context->_func->getVal() ;
  RooAbsReal::setHideOffset(kTRUE) ;
  context->_evalCounter++ ;
  if (RooAbsReal::numEvalErrors()>0 || f>1e30) {

    if (context->_printEvalErrors>=0) {

      if (context->_doEvalErrorWall) {
	oocoutW(context,Minimization) << "RooMinuitGlue: Minimized function has error status." << endl
				      << "Returning maximum FCN so far (" << maxFCN
				      << ") to force MIGRAD to back out of this region. Error log follows" << endl ;
      } else {
	oocoutW(context,Minimization) << "RooMinuitGlue: Minimized function has error status but is ignored" << endl ;
      }

      TIterator* iter = context->_floatParamList->createIterator() ;
      RooRealVar* var ;
      Bool_t first(kTRUE) ;
      ooccoutW(context,Minimization) << "Parameter values: " ;
      while((var=(RooRealVar*)iter->Next())) {
	if (first) { first = kFALSE ; } else ooccoutW(context,Minimization) << ", " ;
	ooccoutW(context,Minimization) << var->GetName() << "=" << var->getVal() ;
      }
      delete iter ;
      ooccoutW(context,Minimization) << endl ;

      RooAbsReal::printEvalErrors(ooccoutW(context,Minimization),context->_printEvalErrors) ;
      ooccoutW(context,Minimization) << endl ;
    }

    if (context->_doEvalErrorWall) {
      f = maxFCN+1 ;
    }

    RooAbsReal::clearEvalErrorLog() ;
    context->_numBadNLL++ ;
  } else if (f>maxFCN) {
    maxFCN = f ;
  }
  
  // Optional logging
  if (logf) (*logf) << setprecision(15) << f << setprecision(4) << endl;
  if (verbose) {
    cout << "\nprevFCN" << (context->_func->isOffsetting()?"-offset":"") << " = " << setprecision(10) << f << setprecision(4) << "  " ;
    cout.flush() ;
  }
}
