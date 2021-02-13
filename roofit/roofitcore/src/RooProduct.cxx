/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 * @(#)root/roofitcore:$Id$
 * Authors:                                                                  *
 *   WV, Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu       *
 *   GR, Gerhard Raven,   VU Amsterdan,     graven@nikhef.nl                 *
 *                                                                           *
 * Copyright (c) 2000-2007, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

/**
\file RooProduct.cxx
\class RooProduct
\ingroup Roofitcore

A RooProduct represents the product of a given set of RooAbsReal objects.

**/


#include <cmath>
#include <memory>

#include "RooProduct.h"
#include "RooNameReg.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
#include "RooErrorHandler.h"
#include "RooMsgService.h"
#include "RooTrace.h"

using namespace std ;

ClassImp(RooProduct);
;

class RooProduct::ProdMap : public  std::vector<std::pair<RooArgSet*,RooArgList*> > {} ;

// Namespace with helper functions that have STL stuff that we don't want to expose to CINT
namespace {
  typedef RooProduct::ProdMap::iterator RPPMIter ;
  std::pair<RPPMIter,RPPMIter> findOverlap2nd(RPPMIter i, RPPMIter end)  ;
  void dump_map(ostream& os, RPPMIter i, RPPMIter end) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Default constructor

RooProduct::RooProduct()
{
  TRACE_CREATE
}



////////////////////////////////////////////////////////////////////////////////
/// Destructor

RooProduct::~RooProduct()
{
  TRACE_DESTROY
}



////////////////////////////////////////////////////////////////////////////////
/// Construct function representing the product of functions in prodSet

RooProduct::RooProduct(const char* name, const char* title, const RooArgList& prodSet) :
  RooAbsReal(name, title),
  _compRSet("!compRSet","Set of real product components",this),
  _compCSet("!compCSet","Set of category product components",this),
  _cacheMgr(this,10)
{
  for (auto comp : prodSet) {
    if (dynamic_cast<RooAbsReal*>(comp)) {
      _compRSet.add(*comp) ;
    } else if (dynamic_cast<RooAbsCategory*>(comp)) {
      _compCSet.add(*comp) ;
    } else {
      coutE(InputArguments) << "RooProduct::ctor(" << GetName() << ") ERROR: component " << comp->GetName() 
			    << " is not of type RooAbsReal or RooAbsCategory" << endl ;
      RooErrorHandler::softAbort() ;
    }
  }
  TRACE_CREATE
}



////////////////////////////////////////////////////////////////////////////////
/// Copy constructor

RooProduct::RooProduct(const RooProduct& other, const char* name) :
  RooAbsReal(other, name), 
  _compRSet("!compRSet",this,other._compRSet),
  _compCSet("!compCSet",this,other._compCSet),
  _cacheMgr(other._cacheMgr,this)
{
  TRACE_CREATE
}



////////////////////////////////////////////////////////////////////////////////
/// Force internal handling of integration of given observable if any
/// of the product terms depend on it.

Bool_t RooProduct::forceAnalyticalInt(const RooAbsArg& dep) const
{
  // Force internal handling of integration of given observable if any
  // of the product terms depend on it.

  RooFIter compRIter = _compRSet.fwdIterator() ;
  RooAbsReal* rcomp ;
  Bool_t depends(kFALSE);
  while((rcomp=(RooAbsReal*)compRIter.next())&&!depends) {
        depends = rcomp->dependsOn(dep);
  }
  return depends ;
}



////////////////////////////////////////////////////////////////////////////////
/// Group observables into subsets in which the product factorizes
/// and that can thus be integrated separately

RooProduct::ProdMap* RooProduct::groupProductTerms(const RooArgSet& allVars) const 
{
  ProdMap* map = new ProdMap ;

  // Do we have any terms which do not depend on the
  // on the variables we integrate over?
  RooAbsReal* rcomp ;
  RooFIter compRIter = _compRSet.fwdIterator() ;
  RooArgList *indep = new RooArgList();
  while((rcomp=(RooAbsReal*) compRIter.next())) {
    if( !rcomp->dependsOn(allVars) ) indep->add(*rcomp);
  }
  if (indep->getSize()!=0) {
    map->push_back( std::make_pair(new RooArgSet(),indep) );
  } else {
     delete indep;
  }

  // Map observables -> functions ; start with individual observables
  RooFIter allVarsIter = allVars.fwdIterator() ;
  RooAbsReal* var ;
  while((var=(RooAbsReal*)allVarsIter.next())) {
    RooArgSet *vars  = new RooArgSet(); vars->add(*var);
    RooArgList *comps = new RooArgList();
    RooAbsReal* rcomp2 ; 
    
    compRIter = _compRSet.fwdIterator() ;
    while((rcomp2=(RooAbsReal*) compRIter.next())) {
      if( rcomp2->dependsOn(*var) ) comps->add(*rcomp2);
    }
    map->push_back( std::make_pair(vars,comps) );
  }

  // Merge groups with overlapping dependents
  Bool_t overlap;
  do {
    std::pair<ProdMap::iterator,ProdMap::iterator> i = findOverlap2nd(map->begin(),map->end());
    overlap = (i.first!=i.second);
    if (overlap) {
      i.first->first->add(*i.second->first);

      // In the merging step, make sure not to duplicate
      RooFIter it = i.second->second->fwdIterator() ;
      RooAbsArg* targ ;
      while ((targ = it.next())) {
	if (!i.first->second->find(*targ)) {
	  i.first->second->add(*targ) ;
	}
      }
      //i.first->second->add(*i.second->second);

      delete i.second->first;
      delete i.second->second;
      map->erase(i.second);
    }
  } while (overlap);
  
  // check that we have all variables to be integrated over on the LHS
  // of the map, and all terms in the product do appear on the RHS
  int nVar=0; int nFunc=0;
  for (ProdMap::iterator i = map->begin();i!=map->end();++i) {
    nVar+=i->first->getSize();
    nFunc+=i->second->getSize();
  }
  assert(nVar==allVars.getSize());
  assert(nFunc==_compRSet.getSize());
  return map;
}



////////////////////////////////////////////////////////////////////////////////
/// Return list of (partial) integrals whose product defines the integral of this
/// RooProduct over the observables in iset in range isetRange. If no such list
/// exists, create it now and store it in the cache for future use.

Int_t RooProduct::getPartIntList(const RooArgSet* iset, const char *isetRange) const
{

  // check if we already have integrals for this combination of factors
  Int_t sterileIndex(-1);
  CacheElem* cache = (CacheElem*) _cacheMgr.getObj(iset,iset,&sterileIndex,RooNameReg::ptr(isetRange));
  if (cache!=0) {
    Int_t code = _cacheMgr.lastIndex();
    return code;
  }
  
  ProdMap* map = groupProductTerms(*iset);

  cxcoutD(Integration) << "RooProduct::getPartIntList(" << GetName() << ") groupProductTerms returned map" ;
  if (dologD(Integration)) {
    dump_map(ccoutD(Integration),map->begin(),map->end()); 
    ccoutD(Integration) << endl;
  }
  
  // did we find any factorizable terms?
  if (map->size()<2) {
    
    for (ProdMap::iterator iter = map->begin() ; iter != map->end() ; ++iter) {
      delete iter->first ;
      delete iter->second ;
    }

    delete map ;
    return -1; // RRI caller will zero analVars if return code = 0....
  }
  cache = new CacheElem();

  for (ProdMap::const_iterator i = map->begin();i!=map->end();++i) {
    RooAbsReal *term(0);
    if (i->second->getSize()>1) { // create a RooProd for this subexpression
      const char *name = makeFPName("SUBPROD_",*i->second);
      term = new RooProduct(name,name,*i->second);
      cache->_ownedList.addOwned(*term);
      cxcoutD(Integration) << "RooProduct::getPartIntList(" << GetName() << ") created subexpression " << term->GetName() << endl;
    } else {
      assert(i->second->getSize()==1);
      RooFIter j = i->second->fwdIterator();
      term = (RooAbsReal*)j.next();
    }
    assert(term!=0);
    if (i->first->getSize()==0) { // check whether we need to integrate over this term or not...
      cache->_prodList.add(*term);
      cxcoutD(Integration) << "RooProduct::getPartIntList(" << GetName() << ") adding simple factor " << term->GetName() << endl;
    } else {
      RooAbsReal *integral = term->createIntegral(*i->first,isetRange);
      cache->_prodList.add(*integral);
      cache->_ownedList.addOwned(*integral);
      cxcoutD(Integration) << "RooProduct::getPartIntList(" << GetName() << ") adding integral for " << term->GetName() << " : " << integral->GetName() << endl;
    }
  }
  // add current set-up to cache, and return index..
  Int_t code = _cacheMgr.setObj(iset,iset,(RooAbsCacheElement*)cache,RooNameReg::ptr(isetRange));

  cxcoutD(Integration) << "RooProduct::getPartIntList(" << GetName() << ") created list " << cache->_prodList << " with code " << code+1 << endl
		       << " for iset=" << *iset << " @" << iset << " range: " << (isetRange?isetRange:"<none>") << endl ;

  for (ProdMap::iterator iter = map->begin() ; iter != map->end() ; ++iter) {
    delete iter->first ;
    delete iter->second ;
  }
  delete map ;
  return code;
}


////////////////////////////////////////////////////////////////////////////////
/// Declare that we handle all integrations internally

Int_t RooProduct::getAnalyticalIntegralWN(RooArgSet& allVars, RooArgSet& analVars,
					  const RooArgSet* /*normSet*/,
					  const char* rangeName) const
{
  if (_forceNumInt) return 0 ;

  // Declare that we can analytically integrate all requested observables
  // (basically, we will take care of the problem, and delegate where required)
  //assert(normSet==0);
  assert(analVars.getSize()==0);
  analVars.add(allVars) ;
  Int_t code = getPartIntList(&analVars,rangeName)+1;
  return code ;
}


////////////////////////////////////////////////////////////////////////////////
/// Calculate integral internally from appropriate partial integral cache

Double_t RooProduct::analyticalIntegral(Int_t code, const char* rangeName) const
{
  // note: rangeName implicit encoded in code: see _cacheMgr.setObj in getPartIntList...
  CacheElem *cache = (CacheElem*) _cacheMgr.getObjByIndex(code-1);
  if (cache==0) { 
    // cache got sterilized, trigger repopulation of this slot, then try again...
    std::unique_ptr<RooArgSet> vars( getParameters(RooArgSet()) );
    std::unique_ptr<RooArgSet> iset(  _cacheMgr.nameSet2ByIndex(code-1)->select(*vars) );
    Int_t code2 = getPartIntList(iset.get(),rangeName)+1;
    assert(code==code2); // must have revived the right (sterilized) slot...
    return analyticalIntegral(code2,rangeName);
  }
  assert(cache!=0);
  
  return calculate(cache->_prodList);
}


////////////////////////////////////////////////////////////////////////////////
/// Calculate and return product of partial terms in partIntList

Double_t RooProduct::calculate(const RooArgList& partIntList) const
{
  RooAbsReal *term(0);
  Double_t val=1;
  RooFIter i = partIntList.fwdIterator() ;
  while((term=(RooAbsReal*)i.next())) {
    double x = term->getVal();
    val*= x;
  }
  return val;
}


////////////////////////////////////////////////////////////////////////////////
/// Construct automatic name for internal product terms

const char* RooProduct::makeFPName(const char *pfx,const RooArgSet& terms) const
{
  static TString pname;
  pname = pfx;
  RooFIter i = terms.fwdIterator();
  RooAbsArg *arg;
  Bool_t first(kTRUE);
  while((arg=(RooAbsArg*)i.next())) {
    if (first) { first=kFALSE;}
    else pname.Append("_X_");
    pname.Append(arg->GetName());
  }
  return pname.Data();
}



////////////////////////////////////////////////////////////////////////////////
/// Evaluate product of input functions

Double_t RooProduct::evaluate() const 
{
  Double_t prod(1) ;

  const RooArgSet* nset = _compRSet.nset() ;
  for (const auto item : _compRSet) {
    auto rcomp = static_cast<const RooAbsReal*>(item);

    prod *= rcomp->getVal(nset) ;
  }
  
  for (const auto item : _compCSet) {
    auto ccomp = static_cast<const RooAbsCategory*>(item);

    prod *= ccomp->getCurrentIndex() ;
  }
  
  return prod ;
}



////////////////////////////////////////////////////////////////////////////////
/// Forward the plot sampling hint from the p.d.f. that defines the observable obs  

std::list<Double_t>* RooProduct::binBoundaries(RooAbsRealLValue& obs, Double_t xlo, Double_t xhi) const
{
  for (const auto item : _compRSet) {
    auto func = static_cast<const RooAbsReal*>(item);

    list<Double_t>* binb = func->binBoundaries(obs,xlo,xhi) ;      
    if (binb) {
      return binb ;
    }
  }
  
  return 0 ;  
}


//_____________________________________________________________________________B
Bool_t RooProduct::isBinnedDistribution(const RooArgSet& obs) const 
{
  // If all components that depend on obs are binned that so is the product
  
  for (const auto item : _compRSet) {
    auto func = static_cast<const RooAbsReal*>(item);

    if (func->dependsOn(obs) && !func->isBinnedDistribution(obs)) {
      return kFALSE ;
    }
  }
  
  return kTRUE  ;  
}



////////////////////////////////////////////////////////////////////////////////
/// Forward the plot sampling hint from the p.d.f. that defines the observable obs  

std::list<Double_t>* RooProduct::plotSamplingHint(RooAbsRealLValue& obs, Double_t xlo, Double_t xhi) const
{
  for (const auto item : _compRSet) {
    auto func = static_cast<const RooAbsReal*>(item);

    list<Double_t>* hint = func->plotSamplingHint(obs,xlo,xhi) ;      
    if (hint) {
      return hint ;
    }
  }
  
  return 0 ;
}



////////////////////////////////////////////////////////////////////////////////
/// Destructor

RooProduct::CacheElem::~CacheElem() 
{
}


////////////////////////////////////////////////////////////////////////////////
/// Return list of all RooAbsArgs in cache element

RooArgList RooProduct::CacheElem::containedArgs(Action) 
{
  RooArgList ret(_ownedList) ;
  return ret ;
}




////////////////////////////////////////////////////////////////////////////////
/// Label OK'ed components of a RooProduct with cache-and-track

void RooProduct::setCacheAndTrackHints(RooArgSet& trackNodes) 
{
  RooArgSet comp(components()) ;
  for (const auto parg : comp) {
    if (parg->isDerived()) {
      if (parg->canNodeBeCached()==Always) {
        trackNodes.add(*parg) ;
	//cout << "tracking node RooProduct component " << parg->IsA()->GetName() << "::" << parg->GetName() << endl ;
      }
    }
  }
}							    





////////////////////////////////////////////////////////////////////////////////
/// Customized printing of arguments of a RooProduct to more intuitively reflect the contents of the
/// product operator construction

void RooProduct::printMetaArgs(ostream& os) const 
{
  Bool_t first(kTRUE) ;

  for (const auto rcomp : _compRSet) {
    if (!first) {  os << " * " ; } else {  first = kFALSE ; }
    os << rcomp->GetName() ;
  }
  
  for (const auto item : _compCSet) {
    auto ccomp = static_cast<const RooAbsCategory*>(item);

    if (!first) {  os << " * " ; } else {  first = kFALSE ; }
    os << ccomp->GetName() ;
  }

  os << " " ;    
}





namespace {

std::pair<RPPMIter,RPPMIter> findOverlap2nd(RPPMIter i, RPPMIter end) 
{
  // Utility function finding pairs of overlapping input functions
  for (; i!=end; ++i) for ( RPPMIter j(i+1); j!=end; ++j) {
    if (i->second->overlaps(*j->second)) {
      return std::make_pair(i,j);
    }
  }
  return std::make_pair(end,end);
}

  
void dump_map(ostream& os, RPPMIter i, RPPMIter end) 
{
  // Utility dump function for debugging
  Bool_t first(kTRUE);
  os << " [ " ;
  for(; i!=end;++i) {
    if (first) { first=kFALSE; }
    else { os << " , " ; }
    os << *(i->first) << " -> " << *(i->second) ;
  }
  os << " ] " ; 
}

}




