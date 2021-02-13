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
\file RooSimGenContext.cxx
\class RooSimGenContext
\ingroup Roofitcore

RooSimGenContext is an efficient implementation of the generator context
specific for RooSimultaneous PDFs when generating more than one of the
component pdfs.
It runs in two modes:
- Proto data with category entries are given: An event from the same category as
in the proto data is created.
- No proto data: A category is chosen randomly.
\note This requires that the PDFs are extended, to determine the relative probabilities
that an event originates from a certain category.
**/

#include "RooFit.h"

#include "RooSimGenContext.h"
#include "RooSimultaneous.h"
#include "RooRealProxy.h"
#include "RooDataSet.h"
#include "Roo1DTable.h"
#include "RooCategory.h"
#include "RooMsgService.h"
#include "RooRandom.h"
#include "RooGlobalFunc.h"

using namespace RooFit;

#include <iostream>
#include <string>

using namespace std;

ClassImp(RooSimGenContext);
;
  

////////////////////////////////////////////////////////////////////////////////
/// Constructor of specialized generator context for RooSimultaneous p.d.f.s. This
/// context creates a dedicated context for each component p.d.f.s and delegates
/// generation of events to the appropriate component generator context

RooSimGenContext::RooSimGenContext(const RooSimultaneous &model, const RooArgSet &vars, 
				   const RooDataSet *prototype, const RooArgSet* auxProto, Bool_t verbose) :
  RooAbsGenContext(model,vars,prototype,auxProto,verbose), _pdf(&model), _protoData(0)
{
  // Determine if we are requested to generate the index category
  RooAbsCategory *idxCat = (RooAbsCategory*) model._indexCat.absArg() ;
  RooArgSet pdfVars(vars) ;

  RooArgSet allPdfVars(pdfVars) ;
  if (prototype) allPdfVars.add(*prototype->get(),kTRUE) ;

  if (!idxCat->isDerived()) {
    pdfVars.remove(*idxCat,kTRUE,kTRUE) ;
    Bool_t doGenIdx = allPdfVars.find(idxCat->GetName())?kTRUE:kFALSE ;

    if (!doGenIdx) {
      oocoutE(_pdf,Generation) << "RooSimGenContext::ctor(" << GetName() << ") ERROR: This context must"
			       << " generate the index category" << endl ;
      _isValid = kFALSE ;
      _numPdf = 0 ;
      _haveIdxProto = kFALSE ;
      return ;
    }
  } else {
    TIterator* sIter = idxCat->serverIterator() ;
    RooAbsArg* server ;
    Bool_t anyServer(kFALSE), allServers(kTRUE) ;
    while((server=(RooAbsArg*)sIter->Next())) {
      if (vars.find(server->GetName())) {
	anyServer=kTRUE ;
	pdfVars.remove(*server,kTRUE,kTRUE) ;
      } else {
	allServers=kFALSE ;
      }
    }
    delete sIter ;    

    if (anyServer && !allServers) {
      oocoutE(_pdf,Generation) << "RooSimGenContext::ctor(" << GetName() << ") ERROR: This context must"
			       << " generate all components of a derived index category" << endl ;
      _isValid = kFALSE ;
      _numPdf = 0 ;
      _haveIdxProto = kFALSE ;
      return ;
    }
  }

  // We must either have the prototype or extended likelihood to determined
  // the relative fractions of the components
  _haveIdxProto = prototype ? kTRUE : kFALSE ;
  _idxCatName = idxCat->GetName() ;
  if (!_haveIdxProto && !model.canBeExtended()) {
    oocoutE(_pdf,Generation) << "RooSimGenContext::ctor(" << GetName() << ") ERROR: Need either extended mode"
			     << " or prototype data to calculate number of events per category" << endl ;
    _isValid = kFALSE ;
    _numPdf = 0 ;
    return ;
  }

  // Initialize fraction threshold array (used only in extended mode)
  _numPdf = model._pdfProxyList.GetSize() ;
  _fracThresh = new Double_t[_numPdf+1] ;
  _fracThresh[0] = 0 ;
  
  // Generate index category and all registered PDFS
  _proxyIter = model._pdfProxyList.MakeIterator() ;
  _allVarsPdf.add(allPdfVars) ;
  RooRealProxy* proxy ;
  RooAbsPdf* pdf ;
  Int_t i(1) ;
  while((proxy=(RooRealProxy*)_proxyIter->Next())) {
    pdf=(RooAbsPdf*)proxy->absArg() ;

    // Create generator context for this PDF
    RooAbsGenContext* cx = pdf->genContext(pdfVars,prototype,auxProto,verbose) ;

    // Name the context after the associated state and add to list
    cx->SetName(proxy->name()) ;
    _gcList.push_back(cx) ;
    _gcIndex.push_back(idxCat->lookupIndex(proxy->name()));

    // Fill fraction threshold array
    _fracThresh[i] = _fracThresh[i-1] + (_haveIdxProto?0:pdf->expectedEvents(&allPdfVars)) ;
    i++ ;
  }   
    
  // Normalize fraction threshold array
  if (!_haveIdxProto) {
    for(i=0 ; i<_numPdf ; i++) 
      _fracThresh[i] /= _fracThresh[_numPdf] ;
  }
  

  // Clone the index category
  _idxCatSet = (RooArgSet*) RooArgSet(model._indexCat.arg()).snapshot(kTRUE) ;
  if (!_idxCatSet) {
    oocoutE(_pdf,Generation) << "RooSimGenContext::RooSimGenContext(" << GetName() << ") Couldn't deep-clone index category, abort," << endl ;
    throw std::string("RooSimGenContext::RooSimGenContext() Couldn't deep-clone index category, abort") ;
  }
  
  _idxCat = (RooAbsCategoryLValue*) _idxCatSet->find(model._indexCat.arg().GetName()) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Destructor. Delete all owned subgenerator contexts

RooSimGenContext::~RooSimGenContext()
{
  delete[] _fracThresh ;
  delete _idxCatSet ;
  for (vector<RooAbsGenContext*>::iterator iter = _gcList.begin() ; iter!=_gcList.end() ; ++iter) {
    delete (*iter) ;
  }
  delete _proxyIter ;
  if (_protoData) delete _protoData ;
}



////////////////////////////////////////////////////////////////////////////////
/// Attach the index category clone to the given event buffer

void RooSimGenContext::attach(const RooArgSet& args) 
{
  if (_idxCat->isDerived()) {
    _idxCat->recursiveRedirectServers(args,kTRUE) ;
  }

  // Forward initGenerator call to all components
  for (vector<RooAbsGenContext*>::iterator iter = _gcList.begin() ; iter!=_gcList.end() ; ++iter) {
    (*iter)->attach(args) ;
  }
  
}


////////////////////////////////////////////////////////////////////////////////
/// Perform one-time initialization of generator context

void RooSimGenContext::initGenerator(const RooArgSet &theEvent)
{
  // Attach the index category clone to the event
  if (_idxCat->isDerived()) {
    _idxCat->recursiveRedirectServers(theEvent,kTRUE) ;
  } else {
    _idxCat = (RooAbsCategoryLValue*) theEvent.find(_idxCat->GetName()) ;
  }
  
  // Update fractions reflecting possible new parameter values
  updateFractions() ;

  // Forward initGenerator call to all components
  for (vector<RooAbsGenContext*>::iterator iter = _gcList.begin() ; iter!=_gcList.end() ; ++iter) {
    (*iter)->initGenerator(theEvent) ;
  }

}


////////////////////////////////////////////////////////////////////////////////
/// Create an empty dataset to hold the events that will be generated

RooDataSet* RooSimGenContext::createDataSet(const char* name, const char* title, const RooArgSet& obs)
{
  
  // If the observables do not contain the index, make a plain dataset
  if (!obs.contains(*_idxCat)) {
    return new RooDataSet(name,title,obs) ;
  }

  if (!_protoData) {
    map<string,RooAbsData*> dmap ;
    for (const auto& nameIdx : *_idxCat) {
      RooAbsPdf* slicePdf = _pdf->getPdf(nameIdx.first.c_str());
      RooArgSet* sliceObs = slicePdf->getObservables(obs) ;
      std::string sliceName = Form("%s_slice_%s", name, nameIdx.first.c_str());
      std::string sliceTitle = Form("%s (index slice %s)", title, nameIdx.first.c_str());
      RooDataSet* dset = new RooDataSet(sliceName.c_str(),sliceTitle.c_str(),*sliceObs) ;
      dmap[nameIdx.first] = dset ;
      delete sliceObs ;
    }
    _protoData = new RooDataSet(name, title, obs, Index((RooCategory&)*_idxCat), Link(dmap), OwnLinked()) ;    
  }

  RooDataSet* emptyClone =  new RooDataSet(*_protoData,name) ;

  return emptyClone ;
}





////////////////////////////////////////////////////////////////////////////////
/// Generate event appropriate for current index state. 
/// The index state is taken either from the prototype
/// or is generated from the fraction threshold table.

void RooSimGenContext::generateEvent(RooArgSet &theEvent, Int_t remaining)
{
  if (_haveIdxProto) {

    // Lookup pdf from selected prototype index state
    Int_t gidx(0), cidx =_idxCat->getCurrentIndex() ;
    for (Int_t i=0 ; i<(Int_t)_gcIndex.size() ; i++) {
      if (_gcIndex[i]==cidx) { gidx = i ; break ; }
    }
    RooAbsGenContext* cx = _gcList[gidx] ;
    if (cx) {      
      cx->generateEvent(theEvent,remaining) ;
    } else {
      oocoutW(_pdf,Generation) << "RooSimGenContext::generateEvent: WARNING, no PDF to generate event of type " << cidx << endl ;
    }    

  
  } else {

    // Throw a random number and select PDF from fraction threshold table
    Double_t rand = RooRandom::uniform() ;
    Int_t i=0 ;
    for (i=0 ; i<_numPdf ; i++) {
      if (rand>_fracThresh[i] && rand<_fracThresh[i+1]) {
        RooAbsGenContext* gen=_gcList[i] ;
        gen->generateEvent(theEvent,remaining) ;
        //Write through to sub-categories because they might be written to dataset:
        _idxCat->setIndex(_gcIndex[i]);
        return ;
      }
    }

  }
}



////////////////////////////////////////////////////////////////////////////////
/// No action needed if we have a proto index

void RooSimGenContext::updateFractions()
{
  if (_haveIdxProto) return ;

  // Generate index category and all registered PDFS
  RooRealProxy* proxy ;
  RooAbsPdf* pdf ;
  Int_t i(1) ;
  _proxyIter->Reset() ;
  while((proxy=(RooRealProxy*)_proxyIter->Next())) {
    pdf=(RooAbsPdf*)proxy->absArg() ;
    
    // Fill fraction threshold array
    _fracThresh[i] = _fracThresh[i-1] + (_haveIdxProto?0:pdf->expectedEvents(&_allVarsPdf)) ;
    i++ ;
  }   
    
  // Normalize fraction threshold array
  if (!_haveIdxProto) {
    for(i=0 ; i<_numPdf ; i++) 
      _fracThresh[i] /= _fracThresh[_numPdf] ;
  }
  
}



////////////////////////////////////////////////////////////////////////////////
/// Set the traversal order of the prototype data to that in the
/// given lookup table. This information is passed to all
/// component generator contexts

void RooSimGenContext::setProtoDataOrder(Int_t* lut)
{
  RooAbsGenContext::setProtoDataOrder(lut) ;

  for (vector<RooAbsGenContext*>::iterator iter = _gcList.begin() ; iter!=_gcList.end() ; ++iter) {
    (*iter)->setProtoDataOrder(lut) ;
  }
}


////////////////////////////////////////////////////////////////////////////////
/// Detailed printing interface

void RooSimGenContext::printMultiline(ostream &os, Int_t content, Bool_t verbose, TString indent) const 
{
  RooAbsGenContext::printMultiline(os,content,verbose,indent) ;
  os << indent << "--- RooSimGenContext ---" << endl ;
  os << indent << "Using PDF ";
  _pdf->printStream(os,kName|kArgs|kClassName,kSingleLine,indent);
  os << indent << "List of component generators" << endl ;

  TString indent2(indent) ;
  indent2.Append("    ") ;

  for (vector<RooAbsGenContext*>::const_iterator iter = _gcList.begin() ; iter!=_gcList.end() ; ++iter) {
    (*iter)->printMultiline(os,content,verbose,indent2);
  }
}
