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
\file RooProjectedPdf.cxx
\class RooProjectedPdf
\ingroup Roofitcore

Class RooProjectedPdf is a RooAbsPdf implementation that represent a projection 
of a given input p.d.f and the object returned by RooAbsPdf::createProjection.
The actual projection integral for it value and normalization are
calculated on the fly in getVal() once the normalization observables are known.
Class RooProjectedPdf can cache projected p.d.f.s for multiple normalization
observables simultaneously.
The createProjection() method of RooProjectedPdf is overloaded and will
return a new RooProjectedPdf that performs the projection of itself
and the requested additional projections in one integration step
The performance of <pre>f->createProjection(x)->createProjection(y)</pre>
is therefore identical to that of <pre>f->createProjection(RooArgSet(x,y))</pre>
**/

#include "Riostream.h" 

#include "RooFit.h"
#include "RooProjectedPdf.h" 
#include "RooMsgService.h"
#include "RooAbsReal.h" 
#include "RooRealVar.h"
#include "RooNameReg.h"

using namespace std;

 ClassImp(RooProjectedPdf); 
   ;


////////////////////////////////////////////////////////////////////////////////
/// Default constructor

RooProjectedPdf::RooProjectedPdf()
{
}



////////////////////////////////////////////////////////////////////////////////
/// Construct projection of input pdf '_intpdf' over observables 'intObs'

 RooProjectedPdf::RooProjectedPdf(const char *name, const char *title, RooAbsReal& _intpdf, const RooArgSet& intObs) :
   RooAbsPdf(name,title), 
   intpdf("!IntegratedPdf","intpdf",this,_intpdf,kFALSE,kFALSE),
   intobs("!IntegrationObservables","intobs",this,kFALSE,kFALSE),
   deps("!Dependents","deps",this,kTRUE,kTRUE),
   _cacheMgr(this,10)
 { 
   intobs.add(intObs) ;

   // Add all other dependens of projected p.d.f. directly
   RooArgSet* tmpdeps = _intpdf.getParameters(intObs) ;
   deps.add(*tmpdeps) ;
   delete tmpdeps ;
 } 



////////////////////////////////////////////////////////////////////////////////
/// Copy constructor

 RooProjectedPdf::RooProjectedPdf(const RooProjectedPdf& other, const char* name) :  
   RooAbsPdf(other,name), 
   intpdf("!IntegratedPdf",this,other.intpdf),
   intobs("!IntegrationObservable",this,other.intobs),
   deps("!Dependents",this,other.deps),
   _cacheMgr(other._cacheMgr,this)
{ 
 } 



////////////////////////////////////////////////////////////////////////////////
/// Evaluate projected p.d.f

Double_t RooProjectedPdf::evaluate() const 
{
  // Calculate current unnormalized value of object
  int code ;
  const RooAbsReal* proj = getProjection(&intobs, _normSet, 0, code);
  
  return proj->getVal() ;
}



////////////////////////////////////////////////////////////////////////////////
/// Retrieve object representing projection integral of input p.d.f
/// over observables iset, while normalizing over observables
/// nset. The code argument returned by reference is the unique code
/// defining this particular projection configuration

const RooAbsReal* RooProjectedPdf::getProjection(const RooArgSet* iset, const RooArgSet* nset, const char* rangeName, int& code) const
{

  // Check if this configuration was created before
  Int_t sterileIdx(-1) ;
  CacheElem* cache = (CacheElem*) _cacheMgr.getObj(iset,nset,&sterileIdx,RooNameReg::ptr(rangeName)) ;
  if (cache) {
    code = _cacheMgr.lastIndex() ;
    return static_cast<const RooAbsReal*>(cache->_projection);
  }

  RooArgSet* nset2 =  intpdf.arg().getObservables(*nset) ;

  if (iset) {
    nset2->add(*iset) ;
  }
  RooAbsReal* proj = intpdf.arg().createIntegral(iset?*iset:RooArgSet(),nset2,0,rangeName) ;
  delete nset2 ;

  cache = new CacheElem ;
  cache->_projection = proj ;

  code = _cacheMgr.setObj(iset,nset,(RooAbsCacheElement*)cache,RooNameReg::ptr(rangeName)) ;

  coutI(Integration) << "RooProjectedPdf::getProjection(" << GetName() << ") creating new projection " << proj->GetName() << " with code " << code << endl ;

  return proj ;
}



////////////////////////////////////////////////////////////////////////////////
/// Special version of RooAbsReal::createProjection that deals with
/// projections of projections. Instead of integrating twice, a new
/// RooProjectedPdf is returned that is configured to perform the
/// complete integration in one step

RooAbsPdf* RooProjectedPdf::createProjection(const RooArgSet& iset) 
{
  RooArgSet combiset(iset) ;
  combiset.add(intobs) ;
  return static_cast<RooAbsPdf&>( const_cast<RooAbsReal&>(intpdf.arg()) ).createProjection(combiset) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Force RooRealIntegral to relegate integration of all observables to internal logic

Bool_t RooProjectedPdf::forceAnalyticalInt(const RooAbsArg& /*dep*/) const 
{
  return kTRUE ;
}



////////////////////////////////////////////////////////////////////////////////
/// Mark all requested variables as internally integrated

Int_t RooProjectedPdf::getAnalyticalIntegralWN(RooArgSet& allVars, RooArgSet& analVars, const RooArgSet* normSet, const char* rangeName) const 
{ 
  analVars.add(allVars) ;
  
  // Create the appropriate integral
  int code ;
  RooArgSet allVars2(allVars) ;
  allVars2.add(intobs) ;
  getProjection(&allVars2,normSet,rangeName,code) ;
  
  return code+1 ; 
} 



////////////////////////////////////////////////////////////////////////////////
/// Return analytical integral represent by appropriate element of projection cache

Double_t RooProjectedPdf::analyticalIntegralWN(Int_t code, const RooArgSet* /*normSet*/, const char* rangeName) const 
{ 
  CacheElem *cache = (CacheElem*) _cacheMgr.getObjByIndex(code-1) ;
  
  if (cache) {
    Double_t ret= cache->_projection->getVal() ;
    return ret ;
  } else {
    
    RooArgSet* vars = getParameters(RooArgSet()) ;
    vars->add(intobs) ;
    RooArgSet* iset = _cacheMgr.nameSet1ByIndex(code-1)->select(*vars) ;
    RooArgSet* nset = _cacheMgr.nameSet2ByIndex(code-1)->select(*vars) ;
    
    Int_t code2(-1) ;
    const RooAbsReal* proj = getProjection(iset,nset,rangeName,code2) ;
    
    delete vars ;
    delete nset ;
    delete iset ;
    
    Double_t ret =  proj->getVal() ;
    return ret ;
  } 
  
} 



////////////////////////////////////////////////////////////////////////////////
/// No internal generator is implemented

Int_t RooProjectedPdf::getGenerator(const RooArgSet& /*directVars*/, RooArgSet& /*generateVars*/, Bool_t /*staticInitOK*/) const 
 { 
   return 0 ; 
 } 



////////////////////////////////////////////////////////////////////////////////
/// No internal generator is implemented

void RooProjectedPdf::generateEvent(Int_t /*code*/) 
 { 
   return; 
 } 



////////////////////////////////////////////////////////////////////////////////
/// Intercept a server redirection all and update list of dependents if necessary 
/// Specifically update the set proxy 'deps' which introduces the dependency
/// on server value dirty flags of ourselves

Bool_t RooProjectedPdf::redirectServersHook(const RooAbsCollection& newServerList, Bool_t /*mustReplaceAll*/, 
				       Bool_t /*nameChange*/, Bool_t /*isRecursive*/) 
{
  // Redetermine explicit list of dependents if intPdf is being replaced
  RooAbsArg* newPdf = newServerList.find(intpdf.arg().GetName()) ;
  if (newPdf) {

    // Determine if set of dependens of new p.d.f is different from old p.d.f.
    RooArgSet olddeps(deps) ;
    RooArgSet* newdeps = newPdf->getParameters(intobs) ;
    RooArgSet* common = (RooArgSet*) newdeps->selectCommon(deps) ;    
    newdeps->remove(*common,kTRUE,kTRUE) ;
    olddeps.remove(*common,kTRUE,kTRUE) ;

    // If so, adjust composition of deps Listproxy
    if (newdeps->getSize()>0) {
      deps.add(*newdeps) ;
    }
    if (olddeps.getSize()>0) {
      deps.remove(olddeps,kTRUE,kTRUE) ;
    }

    delete common ;
    delete newdeps ;
  }

  return kFALSE ;
}



////////////////////////////////////////////////////////////////////////////////
/// Return RooAbsArg elements contained in projection cache element.

RooArgList RooProjectedPdf::CacheElem::containedArgs(Action)
{
  RooArgList ret(*_projection) ;  
  return ret ;
}



////////////////////////////////////////////////////////////////////////////////
/// Customized printing of arguments of a RooRealIntegral to more intuitively reflect the contents of the
/// integration operation

void RooProjectedPdf::printMetaArgs(ostream& os) const
{
  os << "Int " << intpdf.arg().GetName() ;
 
  os << " d" ;
  os << intobs ;
  os << " " ;
  
}




////////////////////////////////////////////////////////////////////////////////
/// Print contents of cache when printing self as part of object tree

void RooProjectedPdf::CacheElem::printCompactTreeHook(ostream& os, const char* indent, Int_t curElem, Int_t maxElem) 
{
  if (curElem==0) {
    os << indent << "RooProjectedPdf begin projection cache" << endl ;
  }

  TString indent2(indent) ;
  indent2 += Form("[%d] ",curElem) ;

  _projection->printCompactTree(os,indent2) ;

  if(curElem==maxElem) {
    os << indent << "RooProjectedPdf end projection cache" << endl ;
  }
}


