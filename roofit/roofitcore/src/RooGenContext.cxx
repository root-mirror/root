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
\file RooGenContext.cxx
\class RooGenContext
\ingroup Roofitcore

Class RooGenContext implement a universal generator context for all
RooAbsPdf classes that do not have or need a specialized generator
context. This generator context queries the input p.d.f which observables
it can generate internally and delegates generation of those observables
to the p.d.f if it deems that safe. The other observables are generated
use a RooAcceptReject sampling technique.
**/


#include "RooFit.h"
#include "RooMsgService.h"
#include "Riostream.h"

#include "RooGenContext.h"
#include "RooAbsPdf.h"
#include "RooDataSet.h"
#include "RooRealIntegral.h"
#include "RooAcceptReject.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooErrorHandler.h"
#include "RooNumGenConfig.h"
#include "RooNumGenFactory.h"

#include "TString.h"
#include "TIterator.h"
#include "TClass.h"



using namespace std;

ClassImp(RooGenContext);
  ;



////////////////////////////////////////////////////////////////////////////////
/// Initialize a new context for generating events with the specified
/// variables, using the specified PDF model. A prototype dataset (if provided)
/// is not cloned and still belongs to the caller. The contents and shape
/// of this dataset can be changed between calls to generate() as long as the
/// expected columns to be copied to the generated dataset are present.
/// Any argument supplied in the forceDirect RooArgSet are always offered
/// for internal generation to the p.d.f., even if this is deemed unsafe by
/// the logic of RooGenContext.

RooGenContext::RooGenContext(const RooAbsPdf &model, const RooArgSet &vars,
			     const RooDataSet *prototype, const RooArgSet* auxProto,
			     Bool_t verbose, const RooArgSet* forceDirect) :  
  RooAbsGenContext(model,vars,prototype,auxProto,verbose),
  _cloneSet(0), _pdfClone(0), _acceptRejectFunc(0), _generator(0),
  _maxVar(0), _uniIter(0), _updateFMaxPerEvent(0) 
{
  cxcoutI(Generation) << "RooGenContext::ctor() setting up event generator context for p.d.f. " << model.GetName() 
			<< " for generation of observable(s) " << vars ;
  if (prototype) ccxcoutI(Generation) << " with prototype data for " << *prototype->get() ;
  if (auxProto && auxProto->getSize()>0)  ccxcoutI(Generation) << " with auxiliary prototypes " << *auxProto ;
  if (forceDirect && forceDirect->getSize()>0)  ccxcoutI(Generation) << " with internal generation forced for observables " << *forceDirect ;
  ccxcoutI(Generation) << endl ;


  // Clone the model and all nodes that it depends on so that this context
  // is independent of any existing objects.
  RooArgSet nodes(model,model.GetName());
  _cloneSet= (RooArgSet*) nodes.snapshot(kTRUE);
  if (!_cloneSet) {
    coutE(Generation) << "RooGenContext::RooGenContext(" << GetName() << ") Couldn't deep-clone PDF, abort," << endl ;
    RooErrorHandler::softAbort() ;
  }

  // Find the clone in the snapshot list
  _pdfClone = (RooAbsPdf*)_cloneSet->find(model.GetName());
  _pdfClone->setOperMode(RooAbsArg::ADirty,kTRUE) ;

  // Optionally fix RooAddPdf normalizations
  if (prototype&&_pdfClone->dependsOn(*prototype->get())) {
    RooArgSet fullNormSet(vars) ;
    fullNormSet.add(*prototype->get()) ;
    _pdfClone->fixAddCoefNormalization(fullNormSet) ;
  }
      
  // Analyze the list of variables to generate...
  _isValid= kTRUE;
  TIterator *iterator= vars.createIterator();
  TIterator *servers= _pdfClone->serverIterator();
  const RooAbsArg *tmp = 0;
  const RooAbsArg *arg = 0;
  while((_isValid && (tmp= (const RooAbsArg*)iterator->Next()))) {
    // is this argument derived?
    if(!tmp->isFundamental()) {
      coutE(Generation) << "RooGenContext::ctor(): cannot generate values for derived \""  << tmp->GetName() << "\"" << endl;
      _isValid= kFALSE;
      continue;
    }
    // lookup this argument in the cloned set of PDF dependents
    arg= (const RooAbsArg*)_cloneSet->find(tmp->GetName());
    if(0 == arg) {
      coutI(Generation) << "RooGenContext::ctor() WARNING model does not depend on \"" << tmp->GetName() 
			  << "\" which will have uniform distribution" << endl;
      _uniformVars.add(*tmp);
    }
    else {

      // does the model depend on this variable directly, ie, like "x" in
      // f(x) or f(x,g(x,y)) or even f(x,x) ?      
      const RooAbsArg *direct= arg ;
      if (forceDirect==0 || !forceDirect->find(direct->GetName())) {
	if (!_pdfClone->isDirectGenSafe(*arg)) {
	  cxcoutD(Generation) << "RooGenContext::ctor() observable " << arg->GetName() << " has been determined to be unsafe for internal generation" << endl;
	  direct=0 ;
	}
      }
      
      // does the model depend indirectly on this variable through an lvalue chain?	
      // otherwise, this variable will have to be generated with accept/reject
      if(direct) { 
	_directVars.add(*arg);
      } else {
	_otherVars.add(*arg);
      }
    }
  }
  delete servers;
  delete iterator;
  if(!isValid()) {
    coutE(Generation) << "RooGenContext::ctor() constructor failed with errors" << endl;
    return;
  }

  // At this point directVars are all variables that are safe to be generated directly
  //               otherVars are all variables that are _not_ safe to be generated directly
  if (_directVars.getSize()>0) {
    cxcoutD(Generation) << "RooGenContext::ctor() observables " << _directVars << " are safe for internal generator (if supported by p.d.f)" << endl ;
  }
  if (_otherVars.getSize()>0) {
    cxcoutD(Generation) << "RooGenContext::ctor() observables " << _otherVars << " are NOT safe for internal generator (if supported by p.d.f)" << endl ;
  }

  // If PDF depends on prototype data, direct generator cannot use static initialization
  // in initGenerator()
  Bool_t staticInitOK = !_pdfClone->dependsOn(_protoVars) ;
  if (!staticInitOK) {
    cxcoutD(Generation) << "RooGenContext::ctor() Model depends on supplied protodata observables, static initialization of internal generator will not be allowed" << endl ;
  }
  
  // Can the model generate any of the direct variables itself?
  RooArgSet generatedVars;
  if (_directVars.getSize()>0) {
    _code= _pdfClone->getGenerator(_directVars,generatedVars,staticInitOK);
    
    cxcoutD(Generation) << "RooGenContext::ctor() Model indicates that it can internally generate observables " 
			  << generatedVars << " with configuration identifier " << _code << endl ;
  } else {
    _code = 0 ;
  }

  // Move variables which cannot be generated into the list to be generated with accept/reject
  _directVars.remove(generatedVars);
  _otherVars.add(_directVars);

  // Update _directVars to only include variables that will actually be directly generated
  _directVars.removeAll();
  _directVars.add(generatedVars);

  cxcoutI(Generation) << "RooGenContext::ctor() Context will" ;
  if (_directVars.getSize()>0) ccxcoutI(Generation) << " let model internally generate observables " << _directVars ;
  if (_directVars.getSize()>0 && _otherVars.getSize()>0)  ccxcoutI(Generation) << " and" ;
  if (_otherVars.getSize()>0)  ccxcoutI(Generation) << " generate variables " << _otherVars << " with accept/reject sampling" ;
  if (_uniformVars.getSize()>0) ccxcoutI(Generation) << ", non-observable variables " << _uniformVars << " will be generated with flat distribution" ;
  ccxcoutI(Generation) << endl ;

  // initialize the accept-reject generator
  RooArgSet *depList= _pdfClone->getObservables(_theEvent);
  depList->remove(_otherVars);

  TString nname(_pdfClone->GetName()) ;
  nname.Append("_AccRej") ;
  TString ntitle(_pdfClone->GetTitle()) ;
  ntitle.Append(" (Accept/Reject)") ;
  

  RooArgSet* protoDeps = model.getObservables(_protoVars) ;
  

  if (_protoVars.getSize()==0) {

    // No prototype variables
    
    if(depList->getSize()==0) {
      // All variable are generated with accept-reject
      
      // Check if PDF supports maximum finding
      Int_t maxFindCode = _pdfClone->getMaxVal(_otherVars) ;
      if (maxFindCode != 0) {
	coutI(Generation) << "RooGenContext::ctor() no prototype data provided, all observables are generated with numerically and "
			      << "model supports analytical maximum findin:, can provide analytical pdf maximum to numeric generator" << endl ;
	Double_t maxVal = _pdfClone->maxVal(maxFindCode) / _pdfClone->getNorm(_theEvent) ;
	_maxVar = new RooRealVar("funcMax","function maximum",maxVal) ;
	cxcoutD(Generation) << "RooGenContext::ctor() maximum value returned by RooAbsPdf::maxVal() is " << maxVal << endl ;
      }
    }

    if (_otherVars.getSize()>0) {
      _pdfClone->getVal(&vars) ; // WVE debug
      _acceptRejectFunc= new RooRealIntegral(nname,ntitle,*_pdfClone,*depList,&vars);
      cxcoutI(Generation) << "RooGenContext::ctor() accept/reject sampling function is " << _acceptRejectFunc->GetName() << endl ;
    } else {
      _acceptRejectFunc = 0 ;
    }

  } else {

    // Generation _with_ prototype variable
    depList->remove(_protoVars,kTRUE,kTRUE) ;
    _acceptRejectFunc= (RooRealIntegral*) _pdfClone->createIntegral(*depList,vars) ;
    cxcoutI(Generation) << "RooGenContext::ctor() accept/reject sampling function is " << _acceptRejectFunc->GetName() << endl ;
    
    // Check if PDF supports maximum finding for the entire phase space
    RooArgSet allVars(_otherVars);
    allVars.add(_directVars);
    Int_t maxFindCode = _pdfClone->getMaxVal(allVars) ;
    if (maxFindCode != 0) {
      // Special case: PDF supports max-finding in otherVars, no need to scan other+proto space for maximum
      coutI(Generation) << "RooGenContext::ctor() prototype data provided, and "
                        << "model supports analytical maximum finding in the full phase space: " 
                        << "can provide analytical pdf maximum to numeric generator" << endl ;
      _maxVar = new RooRealVar("funcMax","function maximum",_pdfClone->maxVal(maxFindCode)) ;
    } else {
      maxFindCode = _pdfClone->getMaxVal(_otherVars) ;
      if (maxFindCode != 0) {
         _updateFMaxPerEvent = maxFindCode ;
         coutI(Generation) << "RooGenContext::ctor() prototype data provided, and "
                           << "model supports analytical maximum finding in the variables which are not"
                           << " internally generated. Can provide analytical pdf maximum to numeric "
                           << "generator" << endl;
         cxcoutD(Generation) << "RooGenContext::ctor() maximum value must be reevaluated for each "
                             << "event with configuration code " << maxFindCode << endl ;
         _maxVar = new RooRealVar("funcMax","function maximum",1) ;    
      }
    }

    if (!_maxVar) {
      
      // Regular case: First find maximum in other+proto space
      RooArgSet otherAndProto(_otherVars) ;

      otherAndProto.add(*protoDeps) ;
      
      if (_otherVars.getSize()>0) {      

	cxcoutD(Generation) << "RooGenContext::ctor() prototype data provided, observables are generated numericaly no " 
			      << "analytical estimate of maximum function value provided by model, must determine maximum value through initial sampling space "
			      << "of accept/reject observables plus prototype observables: " << otherAndProto << endl ;

	// Calculate maximum in other+proto space if there are any accept/reject generated observables
	RooAbsNumGenerator* maxFinder = RooNumGenFactory::instance().createSampler(*_acceptRejectFunc,otherAndProto,RooArgSet(_protoVars),
										   *model.getGeneratorConfig(),_verbose) ;
// 	RooAcceptReject maxFinder(*_acceptRejectFunc,otherAndProto,RooNumGenConfig::defaultConfig(),_verbose) ;
	Double_t max = maxFinder->getFuncMax() ;
	_maxVar = new RooRealVar("funcMax","function maximum",max) ;

	if (max==0) {	  
	  coutE(Generation) << "RooGenContext::ctor(" << model.GetName() 
			    << ") ERROR: generating conditional p.d.f. which requires prior knowledge of function maximum, " 
			    << "but chosen numeric generator (" << maxFinder->IsA()->GetName() << ") does not support maximum finding" << endl ;	
	  delete maxFinder ;	
	  throw string("RooGenContext::ctor()") ;	  
	}	
	delete maxFinder ;	

	cxcoutD(Generation) << "RooGenContext::ctor() maximum function value found through initial sampling is " << max << endl ;
      }
    }
      
  }

  if (_acceptRejectFunc && _otherVars.getSize()>0) {
    _generator = RooNumGenFactory::instance().createSampler(*_acceptRejectFunc,_otherVars,RooArgSet(_protoVars),*model.getGeneratorConfig(),_verbose,_maxVar) ;    
    cxcoutD(Generation) << "RooGenContext::ctor() creating MC sampling generator " << _generator->IsA()->GetName() << "  from function for observables " << _otherVars << endl ;
    //_generator= new RooAcceptReject(*_acceptRejectFunc,_otherVars,RooNumGenConfig::defaultConfig(),_verbose,_maxVar);
  } else {
    _generator = 0 ;
  }

  delete protoDeps ;
  delete depList;
  _otherVars.add(_uniformVars);
}


////////////////////////////////////////////////////////////////////////////////
/// Destructor.

RooGenContext::~RooGenContext() 
{
  // Clean up the cloned objects used in this context.
  delete _cloneSet;

  // Clean up our accept/reject generator
  if (_generator) delete _generator;
  if (_acceptRejectFunc) delete _acceptRejectFunc;
  if (_maxVar) delete _maxVar ;
  if (_uniIter) delete _uniIter ;
}



////////////////////////////////////////////////////////////////////////////////
/// Attach the cloned model to the event buffer we will be filling.

void RooGenContext::attach(const RooArgSet& args) 
{
  _pdfClone->recursiveRedirectServers(args,kFALSE);
  if (_acceptRejectFunc) {
    _acceptRejectFunc->recursiveRedirectServers(args,kFALSE) ; // WVE DEBUG
  }

  // Attach the RooAcceptReject generator the event buffer
  if (_generator) {
    _generator->attachParameters(args) ;
  }

}



////////////////////////////////////////////////////////////////////////////////
/// Perform one-time initialization of the generator context

void RooGenContext::initGenerator(const RooArgSet &theEvent) 
{
  RooFIter iter = theEvent.fwdIterator() ;
  RooAbsArg* arg ;
  while((arg=iter.next())) {
    arg->setOperMode(RooAbsArg::ADirty) ;
  }

  attach(theEvent) ;

  // Reset the cloned model's error counters.
  _pdfClone->resetErrorCounters();

  // Initialize the PDFs internal generator
  if (_directVars.getSize()>0) {
    cxcoutD(Generation) << "RooGenContext::initGenerator() initializing internal generator of model with code " << _code << endl ;
    _pdfClone->initGenerator(_code) ;
  }

  // Create iterator for uniform vars (if any)
  if (_uniformVars.getSize()>0) {
    _uniIter = _uniformVars.createIterator() ;
  }
}


////////////////////////////////////////////////////////////////////////////////
/// Generate one event. The 'remaining' integer is not used other than
/// for printing messages 

void RooGenContext::generateEvent(RooArgSet &theEvent, Int_t remaining) 
{
  if(_otherVars.getSize() > 0) {
    // call the accept-reject generator to generate its variables

    if (_updateFMaxPerEvent!=0) {
      Double_t max = _pdfClone->maxVal(_updateFMaxPerEvent)/_pdfClone->getNorm(_otherVars) ;
      cxcoutD(Generation) << "RooGenContext::initGenerator() reevaluation of maximum function value is required for each event, new value is  " << max << endl ;
      _maxVar->setVal(max) ;
    }

    if (_generator) {
      Double_t resampleRatio(1) ;
      const RooArgSet *subEvent= _generator->generateEvent(remaining,resampleRatio);
      if (resampleRatio<1) {
	coutI(Generation) << "RooGenContext::generateEvent INFO: accept/reject generator requests resampling of previously produced events by factor " 
			  << resampleRatio << " due to increased maximum weight" << endl ; 
	resampleData(resampleRatio) ;
      }
      if(0 == subEvent) {
	coutE(Generation) << "RooGenContext::generateEvent ERROR accept/reject generator failed" << endl ;
	return;
      }
      theEvent.assignValueOnly(*subEvent) ;
      //theEvent= *subEvent;
      
    }
  }

  // Use the model's optimized generator, if one is available.
  // The generator writes directly into our local 'event' since we attached it above.
  if(_directVars.getSize() > 0) {
    _pdfClone->generateEvent(_code);
  }

  // Generate uniform variables (non-dependents)  
  if (_uniIter) {
    _uniIter->Reset() ;
    RooAbsArg* uniVar ;
    while((uniVar=(RooAbsArg*)_uniIter->Next())) {
      RooAbsLValue* arglv = dynamic_cast<RooAbsLValue*>(uniVar) ;
      if (!arglv) {
	coutE(Generation) << "RooGenContext::generateEvent(" << GetName() << ") ERROR: uniform variable " << uniVar->GetName() << " is not an lvalue" << endl ;
	RooErrorHandler::softAbort() ;
      }
      arglv->randomize() ;
    }
    theEvent = _uniformVars ;
  }

}


////////////////////////////////////////////////////////////////////////////////
/// Printing interface

void RooGenContext::printMultiline(ostream &os, Int_t content, Bool_t verbose, TString indent) const
{
  RooAbsGenContext::printMultiline(os,content,verbose,indent);
  os << indent << " --- RooGenContext --- " << endl ;
  os << indent << "Using PDF ";
  _pdfClone->printStream(os,kName|kArgs|kClassName,kSingleLine,indent);
  
  if(verbose) {
    os << indent << "Use PDF generator for " << _directVars << endl ;
    os << indent << "Use MC sampling generator " << (_generator ? _generator->IsA()->GetName() : "<none>") << " for " << _otherVars << endl ;
    if (_protoVars.getSize()>0) {
      os << indent << "Prototype observables are " << _protoVars << endl ;
    }
  }
}
