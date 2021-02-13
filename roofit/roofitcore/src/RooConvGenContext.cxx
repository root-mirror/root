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
\file RooConvGenContext.cxx
\class RooConvGenContext
\ingroup Roofitcore

RooConvGenContext is an efficient implementation of the generator context
specific for RooAbsAnaConvPdf objects. The physics model is generated
with a truth resolution model and the requested resolution model is generated
separately as a PDF. The convolution variable of the physics model is 
subsequently explicitly smeared with the resolution model distribution.
**/

#include "RooMsgService.h"
#include "RooErrorHandler.h"
#include "RooConvGenContext.h"
#include "RooAbsAnaConvPdf.h"
#include "RooNumConvPdf.h"
#include "RooFFTConvPdf.h"
#include "RooProdPdf.h"
#include "RooDataSet.h"
#include "RooArgSet.h"
#include "RooTruthModel.h"
#include "Riostream.h"


using namespace std;

ClassImp(RooConvGenContext);
;
  

////////////////////////////////////////////////////////////////////////////////
/// Constructor for specialized generator context for analytical convolutions. 
/// 
/// Builds a generator for the physics PDF convoluted with the truth model
/// and a generator for the resolution model as PDF. Events are generated
/// by sampling events from the p.d.f and smearings from the resolution model
/// and adding these to obtain a distribution of events consistent with the
/// convolution of these two. The advantage of this procedure is so that
/// both p.d.f and resolution model can take advantage of any internal
/// generators that may be defined.

RooConvGenContext::RooConvGenContext(const RooAbsAnaConvPdf &model, const RooArgSet &vars, 
	 			     const RooDataSet *prototype, const RooArgSet* auxProto, Bool_t verbose) :
  RooAbsGenContext(model,vars,prototype,auxProto,verbose), _pdfVarsOwned(0), _modelVarsOwned(0)
{
  cxcoutI(Generation) << "RooConvGenContext::ctor() setting up special generator context for analytical convolution p.d.f. " << model.GetName() 
		      << " for generation of observable(s) " << vars << endl ;

  // Clone PDF and change model to internal truth model
  _pdfCloneSet = (RooArgSet*) RooArgSet(model).snapshot(kTRUE) ;
  if (!_pdfCloneSet) {
    coutE(Generation) << "RooConvGenContext::RooConvGenContext(" << GetName() << ") Couldn't deep-clone PDF, abort," << endl ;
    RooErrorHandler::softAbort() ;
  }

  RooAbsAnaConvPdf* pdfClone = (RooAbsAnaConvPdf*) _pdfCloneSet->find(model.GetName()) ;
  RooTruthModel truthModel("truthModel","Truth resolution model",*pdfClone->convVar()) ;
  pdfClone->changeModel(truthModel) ;
  auto convV = dynamic_cast<RooRealVar*>(pdfClone->convVar());
  if (!convV) {
    throw std::runtime_error("RooConvGenContext only works with convolution variables of type RooRealVar.");
  }
  convV->removeRange();

  // Create generator for physics X truth model
  _pdfVars = (RooArgSet*) pdfClone->getObservables(&vars) ; ;
  _pdfGen = pdfClone->genContext(*_pdfVars,prototype,auxProto,verbose) ;  

  // Clone resolution model and use as normal PDF
  _modelCloneSet = (RooArgSet*) RooArgSet(*model._convSet.at(0)).snapshot(kTRUE) ;
  if (!_modelCloneSet) {
    coutE(Generation) << "RooConvGenContext::RooConvGenContext(" << GetName() << ") Couldn't deep-clone resolution model, abort," << endl ;
    RooErrorHandler::softAbort() ;
  }
  RooResolutionModel* modelClone = (RooResolutionModel*) 
    _modelCloneSet->find(model._convSet.at(0)->GetName())->Clone("smearing") ;
  _modelCloneSet->addOwned(*modelClone) ;
  modelClone->changeBasis(0) ;
  convV = dynamic_cast<RooRealVar*>(&modelClone->convVar());
  if (!convV) {
    throw std::runtime_error("RooConvGenContext only works with convolution variables of type RooRealVar.");
  }
  convV->removeRange();

  // Create generator for resolution model as PDF
  _modelVars = (RooArgSet*) modelClone->getObservables(&vars) ;

  _modelVars->add(modelClone->convVar()) ;
  _convVarName = modelClone->convVar().GetName() ;
  _modelGen = modelClone->genContext(*_modelVars,prototype,auxProto,verbose) ;

  if (prototype) {
    _pdfVars->add(*prototype->get()) ;
    _modelVars->add(*prototype->get()) ;  
  }

  // WVE ADD FOR DEBUGGING
  if (auxProto) {
    _pdfVars->add(*auxProto) ;
    _modelVars->add(*auxProto) ;
  }

//   cout << "RooConvGenContext::ctor(" << this << "," << GetName() << ") _pdfVars = " << _pdfVars << " "  ; _pdfVars->Print("1") ;
//   cout << "RooConvGenContext::ctor(" << this << "," << GetName() << ") _modelVars = " << _modelVars << " " ; _modelVars->Print("1") ;
}



////////////////////////////////////////////////////////////////////////////////
/// Constructor for specialized generator context for numerical convolutions. 
/// 
/// Builds a generator for the physics PDF convoluted with the truth model
/// and a generator for the resolution model as PDF. Events are generated
/// by sampling events from the p.d.f and smearings from the resolution model
/// and adding these to obtain a distribution of events consistent with the
/// convolution of these two. The advantage of this procedure is so that
/// both p.d.f and resolution model can take advantage of any internal
/// generators that may be defined.

RooConvGenContext::RooConvGenContext(const RooNumConvPdf &model, const RooArgSet &vars, 
	 			     const RooDataSet *prototype, const RooArgSet* auxProto, Bool_t verbose) :
  RooAbsGenContext(model,vars,prototype,auxProto,verbose)
{
  cxcoutI(Generation) << "RooConvGenContext::ctor() setting up special generator context for numeric convolution p.d.f. " << model.GetName() 
			<< " for generation of observable(s) " << vars << endl ;

  // Create generator for physics X truth model
  _pdfVarsOwned = (RooArgSet*) model.conv().clonePdf().getObservables(&vars)->snapshot(kTRUE) ;
  _pdfVars = new RooArgSet(*_pdfVarsOwned) ;
  _pdfGen = ((RooAbsPdf&)model.conv().clonePdf()).genContext(*_pdfVars,prototype,auxProto,verbose) ;  
  _pdfCloneSet = 0 ;

  // Create generator for resolution model as PDF
  _modelVarsOwned = (RooArgSet*) model.conv().cloneModel().getObservables(&vars)->snapshot(kTRUE) ;
  _modelVars = new RooArgSet(*_modelVarsOwned) ;
  _convVarName = model.conv().cloneVar().GetName() ;
  _modelGen = ((RooAbsPdf&)model.conv().cloneModel()).genContext(*_modelVars,prototype,auxProto,verbose) ;
  _modelCloneSet = new RooArgSet ;
  _modelCloneSet->add(model.conv().cloneModel()) ;

  if (prototype) {
    _pdfVars->add(*prototype->get()) ;
    _modelVars->add(*prototype->get()) ;  
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Constructor for specialized generator context for FFT numerical convolutions.
/// 
/// Builds a generator for the physics PDF convoluted with the truth model
/// and a generator for the resolution model as PDF. Events are generated
/// by sampling events from the p.d.f and smearings from the resolution model
/// and adding these to obtain a distribution of events consistent with the
/// convolution of these two. The advantage of this procedure is so that
/// both p.d.f and resolution model can take advantage of any internal
/// generators that may be defined.

RooConvGenContext::RooConvGenContext(const RooFFTConvPdf &model, const RooArgSet &vars, 
	 			     const RooDataSet *prototype, const RooArgSet* auxProto, Bool_t verbose) :
  RooAbsGenContext(model,vars,prototype,auxProto,verbose)
{
  cxcoutI(Generation) << "RooConvGenContext::ctor() setting up special generator context for fft convolution p.d.f. " << model.GetName() 
			<< " for generation of observable(s) " << vars << endl ;

  _convVarName = model._x.arg().GetName() ;

  // Create generator for physics model
  _pdfCloneSet = (RooArgSet*) RooArgSet(model._pdf1.arg()).snapshot(kTRUE) ;
  RooAbsPdf* pdfClone = (RooAbsPdf*) _pdfCloneSet->find(model._pdf1.arg().GetName()) ;
  RooRealVar* cvPdf = (RooRealVar*) _pdfCloneSet->find(model._x.arg().GetName()) ;
  cvPdf->removeRange() ;
  RooArgSet* tmp1 = pdfClone->getObservables(&vars) ;
  _pdfVarsOwned = (RooArgSet*) tmp1->snapshot(kTRUE) ;
  _pdfVars = new RooArgSet(*_pdfVarsOwned) ;
  _pdfGen = pdfClone->genContext(*_pdfVars,prototype,auxProto,verbose) ;  

  // Create generator for resolution model
  _modelCloneSet = (RooArgSet*) RooArgSet(model._pdf2.arg()).snapshot(kTRUE) ;
  RooAbsPdf* modelClone = (RooAbsPdf*) _modelCloneSet->find(model._pdf2.arg().GetName()) ;
  RooRealVar* cvModel = (RooRealVar*) _modelCloneSet->find(model._x.arg().GetName()) ;
  cvModel->removeRange() ;
  RooArgSet* tmp2 = modelClone->getObservables(&vars) ;
  _modelVarsOwned = (RooArgSet*) tmp2->snapshot(kTRUE) ;
  _modelVars = new RooArgSet(*_modelVarsOwned) ;
  _modelGen = modelClone->genContext(*_pdfVars,prototype,auxProto,verbose) ;  

  delete tmp1 ;
  delete tmp2 ;

  if (prototype) {
    _pdfVars->add(*prototype->get()) ;
    _modelVars->add(*prototype->get()) ;  
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Destructor

RooConvGenContext::~RooConvGenContext()
{
  // Destructor. Delete all owned subgenerator contexts
  delete _pdfGen ;
  delete _modelGen ;
  delete _pdfCloneSet ;
  delete _modelCloneSet ;
  delete _modelVars ;
  delete _pdfVars ;
  delete _pdfVarsOwned ;
  delete _modelVarsOwned ;
}



////////////////////////////////////////////////////////////////////////////////
/// Attach given set of arguments to internal clones of
/// pdf and resolution model

void RooConvGenContext::attach(const RooArgSet& args) 
{
  // Find convolution variable in input and output sets
  RooRealVar* cvModel = (RooRealVar*) _modelVars->find(_convVarName) ;
  RooRealVar* cvPdf   = (RooRealVar*) _pdfVars->find(_convVarName) ;

  // Replace all servers in _pdfVars and _modelVars with those in theEvent, except for the convolution variable  
  RooArgSet* pdfCommon = (RooArgSet*) args.selectCommon(*_pdfVars) ;
  pdfCommon->remove(*cvPdf,kTRUE,kTRUE) ;

  RooArgSet* modelCommon = (RooArgSet*) args.selectCommon(*_modelVars) ;
  modelCommon->remove(*cvModel,kTRUE,kTRUE) ;

  _pdfGen->attach(*pdfCommon) ;
  _modelGen->attach(*modelCommon) ;  

  delete pdfCommon ;
  delete modelCommon ;
}


////////////////////////////////////////////////////////////////////////////////
/// One-time initialization of generator context, attaches
/// the context to the supplied event container

void RooConvGenContext::initGenerator(const RooArgSet &theEvent)
{
  // Find convolution variable in input and output sets
  _cvModel = (RooRealVar*) _modelVars->find(_convVarName) ;
  _cvPdf   = (RooRealVar*) _pdfVars->find(_convVarName) ;
  _cvOut   = (RooRealVar*) theEvent.find(_convVarName) ;

  // Replace all servers in _pdfVars and _modelVars with those in theEvent, except for the convolution variable  
  RooArgSet* pdfCommon = (RooArgSet*) theEvent.selectCommon(*_pdfVars) ;
  pdfCommon->remove(*_cvPdf,kTRUE,kTRUE) ;
  _pdfVars->replace(*pdfCommon) ;
  delete pdfCommon ;

  RooArgSet* modelCommon = (RooArgSet*) theEvent.selectCommon(*_modelVars) ;
  modelCommon->remove(*_cvModel,kTRUE,kTRUE) ;
  _modelVars->replace(*modelCommon) ;
  delete modelCommon ;

  // Initialize component generators
  _pdfGen->initGenerator(*_pdfVars) ;
  _modelGen->initGenerator(*_modelVars) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Generate a single event 

void RooConvGenContext::generateEvent(RooArgSet &theEvent, Int_t remaining)
{
  while(1) {

    // Generate pdf and model data
    _modelGen->generateEvent(*_modelVars,remaining) ;
    _pdfGen->generateEvent(*_pdfVars,remaining) ;    
    
    // Construct smeared convolution variable
    Double_t convValSmeared = _cvPdf->getVal() + _cvModel->getVal() ;
    if (_cvOut->isValidReal(convValSmeared)) {
      // Smeared value in acceptance range, transfer values to output set
      theEvent = *_modelVars ;
      theEvent = *_pdfVars ;
      _cvOut->setVal(convValSmeared) ;
      return ;
    }
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Set the traversal order for events in the prototype dataset
/// The argument is a array of integers with a size identical
/// to the number of events in the prototype dataset. Each element
/// should contain an integer in the range 1-N.

void RooConvGenContext::setProtoDataOrder(Int_t* lut)
{
  RooAbsGenContext::setProtoDataOrder(lut) ;
  _modelGen->setProtoDataOrder(lut) ;
  _pdfGen->setProtoDataOrder(lut) ;
}


////////////////////////////////////////////////////////////////////////////////
/// Print the details of this generator context

void RooConvGenContext::printMultiline(ostream &os, Int_t content, Bool_t verbose, TString indent) const 
{
  RooAbsGenContext::printMultiline(os,content,verbose,indent) ;
  os << indent << "--- RooConvGenContext ---" << endl ;
  os << indent << "List of component generators" << endl ;

  TString indent2(indent) ;
  indent2.Append("    ") ;
  
  _modelGen->printMultiline(os,content,verbose,indent2);
  _pdfGen->printMultiline(os,content,verbose,indent2);
}
