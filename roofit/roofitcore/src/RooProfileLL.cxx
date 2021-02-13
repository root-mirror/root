 /***************************************************************************** 
  * Project: RooFit                                                           * 
  *                                                                           * 
  * Copyright (c) 2000-2005, Regents of the University of California          * 
  *                          and Stanford University. All rights reserved.    * 
  *                                                                           * 
  * Redistribution and use in source and binary forms,                        * 
  * with or without modification, are permitted according to the terms        * 
  * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             * 
  *****************************************************************************/ 

/**
\file RooProfileLL.cxx
\class RooProfileLL
\ingroup Roofitcore

Class RooProfileLL implements the profile likelihood estimator for
a given likelihood and set of parameters of interest. The value return by 
RooProfileLL is the input likelihood nll minimized w.r.t all nuisance parameters
(which are all parameters except for those listed in the constructor) minus
the -log(L) of the best fit. Note that this function is slow to evaluate
as a MIGRAD minimization step is executed for each function evaluation
**/

#include "Riostream.h" 

#include "RooFit.h"
#include "RooProfileLL.h" 
#include "RooAbsReal.h" 
#include "RooMinuit.h"
#include "RooMinimizer.h"
#include "RooMsgService.h"
#include "RooRealVar.h"

using namespace std ;

ClassImp(RooProfileLL); 


////////////////////////////////////////////////////////////////////////////////
/// Default constructor 
/// Should only be used by proof. 

 RooProfileLL::RooProfileLL() : 
   RooAbsReal("RooProfileLL","RooProfileLL"), 
   _nll(), 
   _obs("paramOfInterest","Parameters of interest",this), 
   _par("nuisanceParam","Nuisance parameters",this,kFALSE,kFALSE), 
   _startFromMin(kTRUE), 
   _minimizer(0), 
   _absMinValid(kFALSE), 
   _absMin(0),
   _neval(0)
{ 
  _piter = _par.createIterator() ; 
  _oiter = _obs.createIterator() ; 
} 


////////////////////////////////////////////////////////////////////////////////
/// Constructor of profile likelihood given input likelihood nll w.r.t
/// the given set of variables. The input log likelihood is minimized w.r.t
/// to all other variables of the likelihood at each evaluation and the
/// value of the global log likelihood minimum is always subtracted.

RooProfileLL::RooProfileLL(const char *name, const char *title, 
			   RooAbsReal& nllIn, const RooArgSet& observables) :
  RooAbsReal(name,title), 
  _nll("input","-log(L) function",this,nllIn),
  _obs("paramOfInterest","Parameters of interest",this),
  _par("nuisanceParam","Nuisance parameters",this,kFALSE,kFALSE),
  _startFromMin(kTRUE),
  _minimizer(0),
  _absMinValid(kFALSE),
  _absMin(0),
  _neval(0)
{ 
  // Determine actual parameters and observables
  RooArgSet* actualObs = nllIn.getObservables(observables) ;
  RooArgSet* actualPars = nllIn.getParameters(observables) ;

  _obs.add(*actualObs) ;
  _par.add(*actualPars) ;

  delete actualObs ;
  delete actualPars ;

  _piter = _par.createIterator() ;
  _oiter = _obs.createIterator() ;
} 



////////////////////////////////////////////////////////////////////////////////
/// Copy constructor

RooProfileLL::RooProfileLL(const RooProfileLL& other, const char* name) :  
  RooAbsReal(other,name), 
  _nll("nll",this,other._nll),
  _obs("obs",this,other._obs),
  _par("par",this,other._par),
  _startFromMin(other._startFromMin),
  _minimizer(0),
  _absMinValid(kFALSE),
  _absMin(0),
  _paramFixed(other._paramFixed),
  _neval(0)
{ 
  _piter = _par.createIterator() ;
  _oiter = _obs.createIterator() ;

  _paramAbsMin.addClone(other._paramAbsMin) ;
  _obsAbsMin.addClone(other._obsAbsMin) ;
    
} 



////////////////////////////////////////////////////////////////////////////////
/// Destructor

RooProfileLL::~RooProfileLL()
{
  // Delete instance of minuit if it was ever instantiated
  if (_minimizer) {
    delete _minimizer ;
  }

  delete _piter ;
  delete _oiter ;
}




////////////////////////////////////////////////////////////////////////////////

const RooArgSet& RooProfileLL::bestFitParams() const 
{
  validateAbsMin() ;
  return _paramAbsMin ;
}


////////////////////////////////////////////////////////////////////////////////

const RooArgSet& RooProfileLL::bestFitObs() const 
{
  validateAbsMin() ;
  return _obsAbsMin ;
}




////////////////////////////////////////////////////////////////////////////////
/// Optimized implementation of createProfile for profile likelihoods.
/// Return profile of original function in terms of stated parameters 
/// of interest rather than profiling recursively.

RooAbsReal* RooProfileLL::createProfile(const RooArgSet& paramsOfInterest) 
{
  return nll().createProfile(paramsOfInterest) ;
}




////////////////////////////////////////////////////////////////////////////////

void RooProfileLL::initializeMinimizer() const
{
  coutI(Minimization) << "RooProfileLL::evaluate(" << GetName() << ") Creating instance of MINUIT" << endl ;
  
  Bool_t smode = RooMsgService::instance().silentMode() ;
  RooMsgService::instance().setSilentMode(kTRUE) ;
  _minimizer = new MINIMIZER(const_cast<RooAbsReal&>(_nll.arg())) ;
  if (!smode) RooMsgService::instance().setSilentMode(kFALSE) ;
  
  //_minimizer->setPrintLevel(-999) ;
  //_minimizer->setNoWarn() ;
  //_minimizer->setVerbose(1) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Evaluate profile likelihood by minimizing likelihood w.r.t. all
/// parameters that are not considered observables of this profile
/// likelihood object.

Double_t RooProfileLL::evaluate() const 
{ 
  // Instantiate minimizer if we haven't done that already
  if (!_minimizer) {
    initializeMinimizer() ;
  }

  // Save current value of observables
  RooArgSet* obsSetOrig = (RooArgSet*) _obs.snapshot() ;

  validateAbsMin() ;

  
  // Set all observables constant in the minimization
  const_cast<RooSetProxy&>(_obs).setAttribAll("Constant",kTRUE) ;  
  ccoutP(Eval) << "." ; ccoutP(Eval).flush() ;

  // If requested set initial parameters to those corresponding to absolute minimum
  if (_startFromMin) {
    const_cast<RooProfileLL&>(*this)._par = _paramAbsMin ;
  }

  _minimizer->zeroEvalCount() ;

  //TString minim=::ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str();
  //TString algorithm = ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str();
  //if (algorithm == "Migrad") algorithm = "Minimize"; // prefer to use Minimize instead of Migrad
  //_minimizer->minimize(minim.Data(),algorithm.Data());
  
  _minimizer->migrad() ;
  _neval = _minimizer->evalCounter() ;

  // Restore original values and constant status of observables
  TIterator* iter = obsSetOrig->createIterator() ;
  RooRealVar* var ;
  while((var=dynamic_cast<RooRealVar*>(iter->Next()) ) )  {
    RooRealVar* target = (RooRealVar*) _obs.find(var->GetName()) ;
    target->setVal(var->getVal()) ;
    target->setConstant(var->isConstant()) ;
  }
  delete iter ;
  delete obsSetOrig ;

  return _nll - _absMin ; 
} 



////////////////////////////////////////////////////////////////////////////////
/// Check that parameters and likelihood value for 'best fit' are still valid. If not,
/// because the best fit has never been calculated, or because constant parameters have
/// changed value or parameters have changed const/float status, the minimum is recalculated

void RooProfileLL::validateAbsMin() const 
{
  // Check if constant status of any of the parameters have changed
  if (_absMinValid) {
    _piter->Reset() ;
    RooAbsArg* par ;
    while((par=(RooAbsArg*)_piter->Next())) {
      if (_paramFixed[par->GetName()] != par->isConstant()) {
	cxcoutI(Minimization) << "RooProfileLL::evaluate(" << GetName() << ") constant status of parameter " << par->GetName() << " has changed from " 
				<< (_paramFixed[par->GetName()]?"fixed":"floating") << " to " << (par->isConstant()?"fixed":"floating") 
				<< ", recalculating absolute minimum" << endl ;
	_absMinValid = kFALSE ;
	break ;
      }
    }
  }


  // If we don't have the absolute minimum w.r.t all observables, calculate that first
  if (!_absMinValid) {

    cxcoutI(Minimization) << "RooProfileLL::evaluate(" << GetName() << ") determining minimum likelihood for current configurations w.r.t all observable" << endl ;

    
    if (!_minimizer) {
      initializeMinimizer() ;
    }
    
    // Save current values of non-marginalized parameters
    RooArgSet* obsStart = (RooArgSet*) _obs.snapshot(kFALSE) ;

    // Start from previous global minimum 
    if (_paramAbsMin.getSize()>0) {
      const_cast<RooSetProxy&>(_par).assignValueOnly(_paramAbsMin) ;
    }
    if (_obsAbsMin.getSize()>0) {
      const_cast<RooSetProxy&>(_obs).assignValueOnly(_obsAbsMin) ;
    }

    // Find minimum with all observables floating
    const_cast<RooSetProxy&>(_obs).setAttribAll("Constant",kFALSE) ;  
    
    //TString minim=::ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str();
    //TString algorithm = ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str();
    //if (algorithm == "Migrad") algorithm = "Minimize"; // prefer to use Minimize instead of Migrad
    //_minimizer->minimize(minim.Data(),algorithm.Data());
    _minimizer->migrad() ;

    // Save value and remember
    _absMin = _nll ;
    _absMinValid = kTRUE ;

    // Save parameter values at abs minimum as well
    _paramAbsMin.removeAll() ;

    // Only store non-constant parameters here!
    RooArgSet* tmp = (RooArgSet*) _par.selectByAttrib("Constant",kFALSE) ;
    _paramAbsMin.addClone(*tmp) ;
    delete tmp ;

    _obsAbsMin.addClone(_obs) ;

    // Save constant status of all parameters
    _piter->Reset() ;
    RooAbsArg* par ;
    while((par=(RooAbsArg*)_piter->Next())) {
      _paramFixed[par->GetName()] = par->isConstant() ;
    }
    
    if (dologI(Minimization)) {
      cxcoutI(Minimization) << "RooProfileLL::evaluate(" << GetName() << ") minimum found at (" ;

      RooAbsReal* arg ;
      Bool_t first=kTRUE ;
      _oiter->Reset() ;
      while ((arg=(RooAbsReal*)_oiter->Next())) {
	ccxcoutI(Minimization) << (first?"":", ") << arg->GetName() << "=" << arg->getVal() ;	
	first=kFALSE ;
      }      
      ccxcoutI(Minimization) << ")" << endl ;            
    }

    // Restore original parameter values
    const_cast<RooSetProxy&>(_obs) = *obsStart ;
    delete obsStart ;

  }
}



////////////////////////////////////////////////////////////////////////////////

Bool_t RooProfileLL::redirectServersHook(const RooAbsCollection& /*newServerList*/, Bool_t /*mustReplaceAll*/, 
					 Bool_t /*nameChange*/, Bool_t /*isRecursive*/) 
{ 
  if (_minimizer) {
    delete _minimizer ;
    _minimizer = 0 ;
  }
  return kFALSE ;
} 


