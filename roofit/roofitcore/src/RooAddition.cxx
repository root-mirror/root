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
\file RooAddition.cxx
\class RooAddition
\ingroup Roofitcore

RooAddition calculates the sum of a set of RooAbsReal terms, or
when constructed with two sets, it sums the product of the terms
in the two sets.
**/


#include "RooFit.h"

#include "Riostream.h"
#include "RooAddition.h"
#include "RooProduct.h"
#include "RooAbsReal.h"
#include "RooErrorHandler.h"
#include "RooArgSet.h"
#include "RooNameReg.h"
#include "RooNLLVar.h"
#include "RooChi2Var.h"
#include "RooMsgService.h"

#include <algorithm>
#include <cmath>

using namespace std;

ClassImp(RooAddition);
;


////////////////////////////////////////////////////////////////////////////////
/// Empty constructor
RooAddition::RooAddition()
{
}



////////////////////////////////////////////////////////////////////////////////
/// Constructor with a single set consisting of RooAbsReal.
/// \param[in] name Name of the PDF
/// \param[in] title Title
/// \param[in] sumSet The value of the function will be the sum of the values in this set
/// \param[in] takeOwnership If true, the RooAddition object will take ownership of the arguments in `sumSet`

RooAddition::RooAddition(const char* name, const char* title, const RooArgList& sumSet, Bool_t takeOwnership) 
  : RooAbsReal(name, title)
  , _set("!set","set of components",this)
  , _cacheMgr(this,10)
{
  for (const auto comp : sumSet) {
    if (!dynamic_cast<RooAbsReal*>(comp)) {
      coutE(InputArguments) << "RooAddition::ctor(" << GetName() << ") ERROR: component " << comp->GetName() 
			    << " is not of type RooAbsReal" << endl ;
      RooErrorHandler::softAbort() ;
    }
    _set.add(*comp) ;
    if (takeOwnership) _ownedList.addOwned(*comp) ;
  }

}



////////////////////////////////////////////////////////////////////////////////
/// Constructor with two sets of RooAbsReals.
///
/// The sum of pair-wise products of elements in the sets will be computed:
/// \f[
///  A = \sum_i \mathrm{Set1}[i] * \mathrm{Set2}[i]
/// \f]
///
/// \param[in] name Name of the PDF
/// \param[in] title Title
/// \param[in] sumSet1 Left-hand element of the pair-wise products
/// \param[in] sumSet2 Right-hand element of the pair-wise products
/// \param[in] takeOwnership If true, the RooAddition object will take ownership of the arguments in the `sumSets`
///
RooAddition::RooAddition(const char* name, const char* title, const RooArgList& sumSet1, const RooArgList& sumSet2, Bool_t takeOwnership) 
    : RooAbsReal(name, title)
    , _set("!set","set of components",this)
    , _cacheMgr(this,10)
{
  if (sumSet1.getSize() != sumSet2.getSize()) {
    coutE(InputArguments) << "RooAddition::ctor(" << GetName() << ") ERROR: input lists should be of equal length" << endl ;
    RooErrorHandler::softAbort() ;    
  }

  for (unsigned int i = 0; i < sumSet1.size(); ++i) {
    const auto comp1 = &sumSet1[i];
    const auto comp2 = &sumSet2[i];

    if (!dynamic_cast<RooAbsReal*>(comp1)) {
      coutE(InputArguments) << "RooAddition::ctor(" << GetName() << ") ERROR: component " << comp1->GetName() 
			    << " in first list is not of type RooAbsReal" << endl ;
      RooErrorHandler::softAbort() ;
    }

    if (!dynamic_cast<RooAbsReal*>(comp2)) {
      coutE(InputArguments) << "RooAddition::ctor(" << GetName() << ") ERROR: component " << comp2->GetName() 
			    << " in first list is not of type RooAbsReal" << endl ;
      RooErrorHandler::softAbort() ;
    }
    // TODO: add flag to RooProduct c'tor to make it assume ownership...
    TString _name(name);
    _name.Append( "_[");
    _name.Append(comp1->GetName());
    _name.Append( "_x_");
    _name.Append(comp2->GetName());
    _name.Append( "]");
    RooProduct  *prod = new RooProduct( _name, _name , RooArgSet(*comp1, *comp2) /*, takeOwnership */ ) ;
    _set.add(*prod);
    _ownedList.addOwned(*prod) ;
    if (takeOwnership) {
        _ownedList.addOwned(*comp1) ;
        _ownedList.addOwned(*comp2) ;
    }
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Copy constructor

RooAddition::RooAddition(const RooAddition& other, const char* name) 
    : RooAbsReal(other, name)
    , _set("!set",this,other._set)
    , _cacheMgr(other._cacheMgr,this)
{
  // Member _ownedList is intentionally not copy-constructed -- ownership is not transferred
}


////////////////////////////////////////////////////////////////////////////////

RooAddition::~RooAddition() 
{ // Destructor

}

////////////////////////////////////////////////////////////////////////////////
/// Calculate and return current value of self

Double_t RooAddition::evaluate() const 
{
  Double_t sum(0);
  const RooArgSet* nset = _set.nset() ;

//   cout << "RooAddition::eval sum = " ;

  for (const auto arg : _set) {
    const auto comp = static_cast<RooAbsReal*>(arg);
    const Double_t tmp = comp->getVal(nset);
//     cout << tmp << " " ;
    sum += tmp ;
  }
//   cout << " = " << sum << endl ;
  return sum ;
}


////////////////////////////////////////////////////////////////////////////////
/// Return the default error level for MINUIT error analysis
/// If the addition contains one or more RooNLLVars and 
/// no RooChi2Vars, return the defaultErrorLevel() of
/// RooNLLVar. If the addition contains one ore more RooChi2Vars
/// and no RooNLLVars, return the defaultErrorLevel() of
/// RooChi2Var. If the addition contains neither or both
/// issue a warning message and return a value of 1

Double_t RooAddition::defaultErrorLevel() const 
{
  RooAbsReal* nllArg(0) ;
  RooAbsReal* chi2Arg(0) ;

  RooAbsArg* arg ;

  RooArgSet* comps = getComponents() ;
  TIterator* iter = comps->createIterator() ;
  while((arg=(RooAbsArg*)iter->Next())) {
    if (dynamic_cast<RooNLLVar*>(arg)) {
      nllArg = (RooAbsReal*)arg ;
    }
    if (dynamic_cast<RooChi2Var*>(arg)) {
      chi2Arg = (RooAbsReal*)arg ;
    }
  }
  delete iter ;
  delete comps ;

  if (nllArg && !chi2Arg) {
    coutI(Fitting) << "RooAddition::defaultErrorLevel(" << GetName() 
		   << ") Summation contains a RooNLLVar, using its error level" << endl ;
    return nllArg->defaultErrorLevel() ;
  } else if (chi2Arg && !nllArg) {
    coutI(Fitting) << "RooAddition::defaultErrorLevel(" << GetName() 
		   << ") Summation contains a RooChi2Var, using its error level" << endl ;
    return chi2Arg->defaultErrorLevel() ;
  } else if (!nllArg && !chi2Arg) {
    coutI(Fitting) << "RooAddition::defaultErrorLevel(" << GetName() << ") WARNING: "
		   << "Summation contains neither RooNLLVar nor RooChi2Var server, using default level of 1.0" << endl ;
  } else {
    coutI(Fitting) << "RooAddition::defaultErrorLevel(" << GetName() << ") WARNING: "
		   << "Summation contains BOTH RooNLLVar and RooChi2Var server, using default level of 1.0" << endl ;
  }

  return 1.0 ;
}


////////////////////////////////////////////////////////////////////////////////

void RooAddition::enableOffsetting(Bool_t flag) 
{
  for (auto arg : _set) {
    static_cast<RooAbsReal*>(arg)->enableOffsetting(flag) ;
  }  
}



////////////////////////////////////////////////////////////////////////////////

Bool_t RooAddition::setData(RooAbsData& data, Bool_t cloneData) 
{
  for (const auto arg : _set) {
    static_cast<RooAbsReal*>(arg)->setData(data,cloneData) ;
  }  
  return kTRUE ;
}



////////////////////////////////////////////////////////////////////////////////

void RooAddition::printMetaArgs(ostream& os) const 
{
  Bool_t first(kTRUE) ;
  for (const auto arg : _set) {
    if (!first) {
      os << " + " ;
    } else {
      first = kFALSE ;
    }
    os << arg->GetName() ; 
  }  
  os << " " ;    
}

////////////////////////////////////////////////////////////////////////////////

Int_t RooAddition::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const
{
  // we always do things ourselves -- actually, always delegate further down the line ;-)
  analVars.add(allVars);

  // check if we already have integrals for this combination of factors
  Int_t sterileIndex(-1);
  CacheElem* cache = (CacheElem*) _cacheMgr.getObj(&analVars,&analVars,&sterileIndex,RooNameReg::ptr(rangeName));
  if (cache!=0) {
    Int_t code = _cacheMgr.lastIndex();
    return code+1;
  }

  // we don't, so we make it right here....
  cache = new CacheElem;
  for (const auto arg : _set) {// checked in c'tor that this will work...
      RooAbsReal *I = static_cast<const RooAbsReal*>(arg)->createIntegral(analVars,rangeName);
      cache->_I.addOwned(*I);
  }

  Int_t code = _cacheMgr.setObj(&analVars,&analVars,(RooAbsCacheElement*)cache,RooNameReg::ptr(rangeName));
  return 1+code;
}

////////////////////////////////////////////////////////////////////////////////
/// Calculate integral internally from appropriate integral cache

Double_t RooAddition::analyticalIntegral(Int_t code, const char* rangeName) const 
{
  // note: rangeName implicit encoded in code: see _cacheMgr.setObj in getPartIntList...
  CacheElem *cache = (CacheElem*) _cacheMgr.getObjByIndex(code-1);
  if (cache==0) {
    // cache got sterilized, trigger repopulation of this slot, then try again...
    std::unique_ptr<RooArgSet> vars( getParameters(RooArgSet()) );
    std::unique_ptr<RooArgSet> iset(  _cacheMgr.nameSet2ByIndex(code-1)->select(*vars) );
    RooArgSet dummy;
    Int_t code2 = getAnalyticalIntegral(*iset,dummy,rangeName);
    assert(code==code2); // must have revived the right (sterilized) slot...
    return analyticalIntegral(code2,rangeName);
  }
  assert(cache!=0);

  // loop over cache, and sum...
  double result(0);
  for (auto I : cache->_I) {
    result += static_cast<const RooAbsReal*>(I)->getVal();
  }
  return result;

}



////////////////////////////////////////////////////////////////////////////////

std::list<Double_t>* RooAddition::binBoundaries(RooAbsRealLValue& obs, Double_t xlo, Double_t xhi) const
{
  std::list<Double_t>* sumBinB = 0 ;
  Bool_t needClean(kFALSE) ;
  
  RooFIter iter = _set.fwdIterator() ;
  RooAbsReal* func ;
  // Loop over components pdf
  while((func=(RooAbsReal*)iter.next())) {

    std::list<Double_t>* funcBinB = func->binBoundaries(obs,xlo,xhi) ;
    
    // Process hint
    if (funcBinB) {
      if (!sumBinB) {
	// If this is the first hint, then just save it
	sumBinB = funcBinB ;
      } else {
	
	std::list<Double_t>* newSumBinB = new std::list<Double_t>(sumBinB->size()+funcBinB->size()) ;

	// Merge hints into temporary array
	merge(funcBinB->begin(),funcBinB->end(),sumBinB->begin(),sumBinB->end(),newSumBinB->begin()) ;
	
	// Copy merged array without duplicates to new sumBinBArrau
	delete sumBinB ;
	delete funcBinB ;
	sumBinB = newSumBinB ;
	needClean = kTRUE ;	
      }
    }
  }

  // Remove consecutive duplicates
  if (needClean) {
    std::list<Double_t>::iterator new_end = unique(sumBinB->begin(),sumBinB->end()) ;
    sumBinB->erase(new_end,sumBinB->end()) ;
  }

  return sumBinB ;
}


//_____________________________________________________________________________B
Bool_t RooAddition::isBinnedDistribution(const RooArgSet& obs) const 
{
  // If all components that depend on obs are binned that so is the product
  
  RooFIter iter = _set.fwdIterator() ;
  RooAbsReal* func ;
  while((func=(RooAbsReal*)iter.next())) {
    if (func->dependsOn(obs) && !func->isBinnedDistribution(obs)) {
      return kFALSE ;
    }
  }
  
  return kTRUE  ;  
}




////////////////////////////////////////////////////////////////////////////////

std::list<Double_t>* RooAddition::plotSamplingHint(RooAbsRealLValue& obs, Double_t xlo, Double_t xhi) const
{
  std::list<Double_t>* sumHint = 0 ;
  Bool_t needClean(kFALSE) ;
  
  RooFIter iter = _set.fwdIterator() ;
  RooAbsReal* func ;
  // Loop over components pdf
  while((func=(RooAbsReal*)iter.next())) {
    
    std::list<Double_t>* funcHint = func->plotSamplingHint(obs,xlo,xhi) ;
    
    // Process hint
    if (funcHint) {
      if (!sumHint) {

	// If this is the first hint, then just save it
	sumHint = funcHint ;

      } else {
	
	std::list<Double_t>* newSumHint = new std::list<Double_t>(sumHint->size()+funcHint->size()) ;
	
	// Merge hints into temporary array
	merge(funcHint->begin(),funcHint->end(),sumHint->begin(),sumHint->end(),newSumHint->begin()) ;

	// Copy merged array without duplicates to new sumHintArrau
	delete sumHint ;
	sumHint = newSumHint ;
	needClean = kTRUE ;	
      }
    }
  }

  // Remove consecutive duplicates
  if (needClean) {
    std::list<Double_t>::iterator new_end = unique(sumHint->begin(),sumHint->end()) ;
    sumHint->erase(new_end,sumHint->end()) ;
  }

  return sumHint ;
}



////////////////////////////////////////////////////////////////////////////////
/// Return list of all RooAbsArgs in cache element

RooArgList RooAddition::CacheElem::containedArgs(Action)
{
  RooArgList ret(_I) ;
  return ret ;
}

RooAddition::CacheElem::~CacheElem()
{
  // Destructor
}


