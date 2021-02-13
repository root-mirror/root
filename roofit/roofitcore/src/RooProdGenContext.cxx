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
\file RooProdGenContext.cxx
\class RooProdGenContext
\ingroup Roofitcore

RooProdGenContext is an efficient implementation of the generator context
specific for RooProdPdf PDFs. The sim-context owns a list of
component generator contexts that are used to generate the dependents
for each component PDF sequentially. 
**/

#include "RooFit.h"
#include "Riostream.h"
#include "RooMsgService.h"

#include "RooProdGenContext.h"
#include "RooProdPdf.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooGlobalFunc.h"



using namespace std;

ClassImp(RooProdGenContext);
;
  

////////////////////////////////////////////////////////////////////////////////

RooProdGenContext::RooProdGenContext(const RooProdPdf &model, const RooArgSet &vars, 
				     const RooDataSet *prototype, const RooArgSet* auxProto, Bool_t verbose) :
  RooAbsGenContext(model,vars,prototype,auxProto,verbose), _uniIter(0), _pdf(&model)
{
  // Constructor of optimization generator context for RooProdPdf objects

  //Build an array of generator contexts for each product component PDF
  cxcoutI(Generation) << "RooProdGenContext::ctor() setting up event special generator context for product p.d.f. " << model.GetName() 
			<< " for generation of observable(s) " << vars ;
  if (prototype) ccxcoutI(Generation) << " with prototype data for " << *prototype->get() ;
  if (auxProto && auxProto->getSize()>0)  ccxcoutI(Generation) << " with auxiliary prototypes " << *auxProto ;
  ccxcoutI(Generation) << endl ;

  // Make full list of dependents (generated & proto)
  RooArgSet deps(vars) ;
  if (prototype) {
    RooArgSet* protoDeps = model.getObservables(*prototype->get()) ;
    deps.remove(*protoDeps,kTRUE,kTRUE) ;
    delete protoDeps ;
  }

  // Factorize product in irreducible terms
  RooLinkedList termList,depsList,impDepList,crossDepList,intList ;
  model.factorizeProduct(deps,RooArgSet(),termList,depsList,impDepList,crossDepList,intList) ;
  TIterator* termIter = termList.MakeIterator() ;
  TIterator* normIter = depsList.MakeIterator() ;
  TIterator* impIter = impDepList.MakeIterator() ;
  
  if (dologD(Generation)) {
    cxcoutD(Generation) << "RooProdGenContext::ctor() factorizing product expression in irriducible terms " ;    
    while(RooArgSet* t=(RooArgSet*)termIter->Next()) {
      ccxcoutD(Generation) << *t ;
    }
    ccxcoutD(Generation) << endl ;
  }

  RooArgSet genDeps ;
  // First add terms that do not import observables
  
  Bool_t anyAction = kTRUE ;
  Bool_t go=kTRUE ; 
  while(go) {

    RooAbsPdf* pdf ;
    RooArgSet* term ;
    RooArgSet* impDeps ;
    RooArgSet* termDeps ;

    termIter->Reset() ;
    impIter->Reset() ;
    normIter->Reset() ;

    Bool_t anyPrevAction=anyAction ;
    anyAction=kFALSE ;

    if (termList.GetSize()==0) {
      break ;
    }

    while((term=(RooArgSet*)termIter->Next())) {

      impDeps = (RooArgSet*)impIter->Next() ;
      termDeps = (RooArgSet*)normIter->Next() ;
      if (impDeps==0 || termDeps==0) {
	break ;
      }

      cxcoutD(Generation) << "RooProdGenContext::ctor() analyzing product term " << *term << " with observable(s) " << *termDeps ;
      if (impDeps->getSize()>0) {
	ccxcoutD(Generation) << " which has dependence of external observable(s) " << *impDeps << " that to be generated first by other terms" ;
      }
      ccxcoutD(Generation) << endl ;

      // Add this term if we have no imported dependents, or imported dependents are already generated
      RooArgSet neededDeps(*impDeps) ;
      neededDeps.remove(genDeps,kTRUE,kTRUE) ;

      if (neededDeps.getSize()>0) {
	if (!anyPrevAction) {
	  cxcoutD(Generation) << "RooProdGenContext::ctor() no convergence in single term analysis loop, terminating loop and process remainder of terms as single unit " << endl ;
	  go=kFALSE ;
	  break ;
	}
	cxcoutD(Generation) << "RooProdGenContext::ctor() skipping this term for now because it needs imported dependents that are not generated yet" << endl ;	
	continue ;
      }

      // Check if this component has any dependents that need to be generated
      // e.g. it can happen that there are none if all dependents of this component are prototyped
      if (termDeps->getSize()==0) {
	cxcoutD(Generation) << "RooProdGenContext::ctor() term has no observables requested to be generated, removing it" << endl ;
	termList.Remove(term) ;
	depsList.Remove(termDeps) ;
	impDepList.Remove(impDeps) ;
	delete term ;
	delete termDeps ;
	delete impDeps ;
	anyAction=kTRUE ;
	continue ;
      }

      TIterator* pdfIter = term->createIterator() ;      
      if (term->getSize()==1) {
	// Simple term
	
	pdf = (RooAbsPdf*) pdfIter->Next() ;
	RooArgSet* pdfDep = pdf->getObservables(termDeps) ;
	if (pdfDep->getSize()>0) {
 	  coutI(Generation) << "RooProdGenContext::ctor() creating subcontext for generation of observables " << *pdfDep << " from model " << pdf->GetName() << endl ;
	  RooArgSet* auxProto2 = pdf->getObservables(impDeps) ;
	  RooAbsGenContext* cx = pdf->genContext(*pdfDep,prototype,auxProto2,verbose) ;
	  delete auxProto2 ;
	  _gcList.push_back(cx) ;
	} 

// 	cout << "adding following dependents to list of generated observables: " ; pdfDep->Print("1") ;
	genDeps.add(*pdfDep) ;

	delete pdfDep ;
	
      } else {
	
	// Composite term
	if (termDeps->getSize()>0) {
	  const std::string name = model.makeRGPPName("PRODGEN_",*term,RooArgSet(),RooArgSet(),0) ;      
	  
	  // Construct auxiliary PDF expressing product of composite terms, 
	  // following Conditional component specification of input model
	  RooLinkedList cmdList ;
	  RooLinkedList pdfSetList ;
	  pdfIter->Reset() ;
	  RooArgSet fullPdfSet ;
	  while((pdf=(RooAbsPdf*)pdfIter->Next())) {

	    RooArgSet* pdfnset = model.findPdfNSet(*pdf) ;
	    RooArgSet* pdfSet = new RooArgSet(*pdf) ;
	    pdfSetList.Add(pdfSet) ;

	    if (pdfnset && pdfnset->getSize()>0) {
	      // This PDF requires a Conditional() construction
	      cmdList.Add(RooFit::Conditional(*pdfSet,*pdfnset).Clone()) ;
//   	      cout << "Conditional " << pdf->GetName() << " " ; pdfnset->Print("1") ;
	    } else {
	      fullPdfSet.add(*pdfSet) ;
	    }
	    
	  }
	  RooProdPdf* multiPdf = new RooProdPdf(name.c_str(),name.c_str(),fullPdfSet,cmdList) ;
	  cmdList.Delete() ;
	  pdfSetList.Delete() ;

	  multiPdf->setOperMode(RooAbsArg::ADirty,kTRUE) ;
	  multiPdf->useDefaultGen(kTRUE) ;
	  _ownedMultiProds.addOwned(*multiPdf) ;
	  
   	  coutI(Generation) << "RooProdGenContext()::ctor creating subcontext for generation of observables " << *termDeps 
			      << "for irriducuble composite term using sub-product object " << multiPdf->GetName() ;
	  RooAbsGenContext* cx = multiPdf->genContext(*termDeps,prototype,auxProto,verbose) ;
	  _gcList.push_back(cx) ;

	  genDeps.add(*termDeps) ;

	}
      }
      
      delete pdfIter ;

//        cout << "added generator for this term, removing from list" << endl ;

      termList.Remove(term) ;
      depsList.Remove(termDeps) ;
      impDepList.Remove(impDeps) ;
      delete term ;
      delete termDeps ;
      delete impDeps ;
      anyAction=kTRUE ;
    }
  }

  // Check if there are any left over terms that cannot be generated 
  // separately due to cross dependency of observables
  if (termList.GetSize()>0) {

    cxcoutD(Generation) << "RooProdGenContext::ctor() there are left-over terms that need to be generated separately" << endl ;
 
    RooAbsPdf* pdf ;
    RooArgSet* term ;

    // Concatenate remaining terms
    termIter->Reset() ;
    normIter->Reset() ;
    RooArgSet trailerTerm ;
    RooArgSet trailerTermDeps ;
    while((term=(RooArgSet*)termIter->Next())) {
      RooArgSet* termDeps = (RooArgSet*)normIter->Next() ;
      trailerTerm.add(*term) ;
      trailerTermDeps.add(*termDeps) ;
    }

    const std::string name = model.makeRGPPName("PRODGEN_",trailerTerm,RooArgSet(),RooArgSet(),0) ;      
      
    // Construct auxiliary PDF expressing product of composite terms, 
    // following Partial/Full component specification of input model
    RooLinkedList cmdList ;
    RooLinkedList pdfSetList ;
    RooArgSet fullPdfSet ;

    TIterator* pdfIter = trailerTerm.createIterator() ;
    while((pdf=(RooAbsPdf*)pdfIter->Next())) {
	
      RooArgSet* pdfnset = model.findPdfNSet(*pdf) ;
      RooArgSet* pdfSet = new RooArgSet(*pdf) ;
      pdfSetList.Add(pdfSet) ;
      
      if (pdfnset && pdfnset->getSize()>0) {
	// This PDF requires a Conditional() construction
	  cmdList.Add(RooFit::Conditional(*pdfSet,*pdfnset).Clone()) ;
      } else {
	fullPdfSet.add(*pdfSet) ;
      }
      
    }
//     cmdList.Print("v") ;
    RooProdPdf* multiPdf = new RooProdPdf(name.c_str(),name.c_str(),fullPdfSet,cmdList) ;
    cmdList.Delete() ;
    pdfSetList.Delete() ;
   
    multiPdf->setOperMode(RooAbsArg::ADirty,kTRUE) ;
    multiPdf->useDefaultGen(kTRUE) ;
    _ownedMultiProds.addOwned(*multiPdf) ;
    
    cxcoutD(Generation) << "RooProdGenContext(" << model.GetName() << "): creating context for irreducible composite trailer term " 
	 << multiPdf->GetName() << " that generates observables " << trailerTermDeps << endl ;
    RooAbsGenContext* cx = multiPdf->genContext(trailerTermDeps,prototype,auxProto,verbose) ;
    _gcList.push_back(cx) ;    
  }

  // Now check if the are observables in vars that are not generated by any of the above p.d.f.s
  // If not, generate uniform distributions for these using a special context
  _uniObs.add(vars) ;
  _uniObs.remove(genDeps,kTRUE,kTRUE) ;
  if (_uniObs.getSize()>0) {
    _uniIter = _uniObs.createIterator() ;
    coutI(Generation) << "RooProdGenContext(" << model.GetName() << "): generating uniform distribution for non-dependent observable(s) " << _uniObs << endl ;
  }


  delete termIter ;
  delete impIter ;
  delete normIter ;


  // We own contents of lists filled by factorizeProduct() 
  termList.Delete() ;
  depsList.Delete() ;
  impDepList.Delete() ;
  crossDepList.Delete() ;
  intList.Delete() ;

}



////////////////////////////////////////////////////////////////////////////////
/// Destructor. Delete all owned subgenerator contexts

RooProdGenContext::~RooProdGenContext()
{
  delete _uniIter ;
  for (list<RooAbsGenContext*>::iterator iter=_gcList.begin() ; iter!=_gcList.end() ; ++iter) {
    delete (*iter) ;
  }
}


////////////////////////////////////////////////////////////////////////////////
/// Attach generator to given event buffer

void RooProdGenContext::attach(const RooArgSet& args) 
{
  //Forward initGenerator call to all components
  for (list<RooAbsGenContext*>::iterator iter=_gcList.begin() ; iter!=_gcList.end() ; ++iter) {
    (*iter)->attach(args) ;
  }
}


////////////////////////////////////////////////////////////////////////////////
/// One-time initialization of generator context, forward to component generators

void RooProdGenContext::initGenerator(const RooArgSet &theEvent)
{
  // Forward initGenerator call to all components
  for (list<RooAbsGenContext*>::iterator iter=_gcList.begin() ; iter!=_gcList.end() ; ++iter) {
    (*iter)->initGenerator(theEvent) ;
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Generate a single event of the product by generating the components
/// of the products sequentially. The subcontext have been order such
/// that all conditional dependencies are correctly taken into account
/// when processed in sequential order

void RooProdGenContext::generateEvent(RooArgSet &theEvent, Int_t remaining)
{
  // Loop over the component generators

  for (list<RooAbsGenContext*>::iterator iter=_gcList.begin() ; iter!=_gcList.end() ; ++iter) {
    (*iter)->generateEvent(theEvent,remaining) ;
  }

  // Generate uniform variables (non-dependents)  
  if (_uniIter) {
    _uniIter->Reset() ;
    RooAbsArg* uniVar ;
    while((uniVar=(RooAbsArg*)_uniIter->Next())) {
      RooAbsLValue* arglv = dynamic_cast<RooAbsLValue*>(uniVar) ;
      if (arglv) {
	arglv->randomize() ;
      }
    }
    theEvent = _uniObs ;
  }  

}



////////////////////////////////////////////////////////////////////////////////
/// Set the traversal order of the prototype dataset by the
/// given lookup table

void RooProdGenContext::setProtoDataOrder(Int_t* lut)
{
  // Forward call to component generators
  RooAbsGenContext::setProtoDataOrder(lut) ;

  for (list<RooAbsGenContext*>::iterator iter=_gcList.begin() ; iter!=_gcList.end() ; ++iter) {
    (*iter)->setProtoDataOrder(lut) ;
  }

}



////////////////////////////////////////////////////////////////////////////////
/// Detailed printing interface

void RooProdGenContext::printMultiline(ostream &os, Int_t content, Bool_t verbose, TString indent) const 
{
  RooAbsGenContext::printMultiline(os,content,verbose,indent) ;
  os << indent << "--- RooProdGenContext ---" << endl ;
  os << indent << "Using PDF ";
  _pdf->printStream(os,kName|kArgs|kClassName,kSingleLine,indent);
  os << indent << "List of component generators" << endl ;

  TString indent2(indent) ;
  indent2.Append("    ") ;

  for (list<RooAbsGenContext*>::const_iterator iter=_gcList.begin() ; iter!=_gcList.end() ; ++iter) {
    (*iter)->printMultiline(os,content,verbose,indent2) ;
  }
}
