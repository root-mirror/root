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
\file RooVectorDataStore.cxx
\class RooVectorDataStore
\ingroup Roofitcore

RooVectorDataStore is the abstract base class for data collection that
use a TTree as internal storage mechanism
**/

#include "RooFit.h"
#include "RooMsgService.h"
#include "RooVectorDataStore.h"
#include "RooTreeDataStore.h"

#include "Riostream.h"
#include "TTree.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TROOT.h"
#include "RooFormulaVar.h"
#include "RooRealVar.h"
#include "RooCategory.h"
#include "RooNameSet.h"
#include "RooHistError.h"
#include "RooTrace.h"

#include <iomanip>
using namespace std ;

ClassImp(RooVectorDataStore);
ClassImp(RooVectorDataStore::RealVector);
;





////////////////////////////////////////////////////////////////////////////////

RooVectorDataStore::RooVectorDataStore() :
  _wgtVar(0),
  _nRealF(0),
  _nCat(0),
  _nEntries(0),	 
  _firstRealF(0),
  _firstCat(0),
  _sumWeight(0),
  _sumWeightCarry(0),
  _extWgtArray(0),
  _extWgtErrLoArray(0),
  _extWgtErrHiArray(0),
  _extSumW2Array(0),
  _curWgt(1),
  _curWgtErrLo(0),
  _curWgtErrHi(0),
  _curWgtErr(0),
  _cache(0),
  _cacheOwner(0),
  _forcedUpdate(kFALSE)
{
  TRACE_CREATE
}



////////////////////////////////////////////////////////////////////////////////

RooVectorDataStore::RooVectorDataStore(const char* name, const char* title, const RooArgSet& vars, const char* wgtVarName) :
  RooAbsDataStore(name,title,varsNoWeight(vars,wgtVarName)),
  _varsww(vars),
  _wgtVar(weightVar(vars,wgtVarName)),
  _nRealF(0),
  _nCat(0),
  _nEntries(0),	   
  _firstRealF(0),
  _firstCat(0),
  _sumWeight(0),
  _sumWeightCarry(0),
  _extWgtArray(0),
  _extWgtErrLoArray(0),
  _extWgtErrHiArray(0),
  _extSumW2Array(0),
  _curWgt(1),
  _curWgtErrLo(0),
  _curWgtErrHi(0),
  _curWgtErr(0),
  _cache(0),
  _cacheOwner(0),
  _forcedUpdate(kFALSE)
{
  for (auto arg : _varsww) {
    arg->attachToVStore(*this) ;
  }
  
  setAllBuffersNative() ;
  TRACE_CREATE
}



////////////////////////////////////////////////////////////////////////////////

void RooVectorDataStore::setAllBuffersNative()
{
  for (auto realVec : _realStoreList) {
    realVec->setNativeBuffer();
  }

  for (auto fullVec : _realfStoreList) {
    fullVec->setNativeBuffer();
  }

  for (auto catVec : _catStoreList) {
    catVec->setNativeBuffer();
  }
}




////////////////////////////////////////////////////////////////////////////////
/// Utility function for constructors
/// Return RooArgSet that is copy of allVars minus variable matching wgtName if specified

RooArgSet RooVectorDataStore::varsNoWeight(const RooArgSet& allVars, const char* wgtName) 
{
  RooArgSet ret(allVars) ;
  if(wgtName) {
    RooAbsArg* wgt = allVars.find(wgtName) ;
    if (wgt) {
      ret.remove(*wgt,kTRUE,kTRUE) ;
    }
  }
  return ret ;
}



////////////////////////////////////////////////////////////////////////////////
/// Utility function for constructors
/// Return pointer to weight variable if it is defined

RooRealVar* RooVectorDataStore::weightVar(const RooArgSet& allVars, const char* wgtName) 
{
  if(wgtName) {
    RooRealVar* wgt = dynamic_cast<RooRealVar*>(allVars.find(wgtName)) ;
    return wgt ;
  } 
  return 0 ;
}




////////////////////////////////////////////////////////////////////////////////
/// Regular copy ctor

RooVectorDataStore::RooVectorDataStore(const RooVectorDataStore& other, const char* newname) :
  RooAbsDataStore(other,newname), 
  _varsww(other._varsww),
  _wgtVar(other._wgtVar),
  _nRealF(0),
  _nCat(0),
  _nEntries(other._nEntries),	 
  _sumWeight(other._sumWeight),
  _sumWeightCarry(other._sumWeightCarry),
  _extWgtArray(other._extWgtArray),
  _extWgtErrLoArray(other._extWgtErrLoArray),
  _extWgtErrHiArray(other._extWgtErrHiArray),
  _extSumW2Array(other._extSumW2Array),
  _curWgt(other._curWgt),
  _curWgtErrLo(other._curWgtErrLo),
  _curWgtErrHi(other._curWgtErrHi),
  _curWgtErr(other._curWgtErr),
  _cache(0),
  _cacheOwner(0),
  _forcedUpdate(kFALSE)
{
  for (const auto realVec : other._realStoreList) {
    _realStoreList.push_back(new RealVector(*realVec, (RooAbsReal*)_varsww.find(realVec->_nativeReal->GetName()))) ;
  }

  for (const auto realFullVec : other._realfStoreList) {
    _realfStoreList.push_back(new RealFullVector(*realFullVec, (RooAbsReal*)_varsww.find(realFullVec->_nativeReal->GetName()))) ;
    _nRealF++ ;
  }

  for (const auto catVec : other._catStoreList) {
    _catStoreList.push_back(new CatVector(*catVec, (RooAbsCategory*)_varsww.find(catVec->_cat->GetName()))) ;
    _nCat++ ;
 }

  setAllBuffersNative() ;
  
  _firstRealF = _realfStoreList.size()>0 ? &_realfStoreList.front() : 0 ;
  _firstCat = _catStoreList.size()>0 ? &_catStoreList.front() : 0 ;
  TRACE_CREATE
}


////////////////////////////////////////////////////////////////////////////////

RooVectorDataStore::RooVectorDataStore(const RooTreeDataStore& other, const RooArgSet& vars, const char* newname) :
  RooAbsDataStore(other,varsNoWeight(vars,other._wgtVar?other._wgtVar->GetName():0),newname),
  _varsww(vars),
  _wgtVar(weightVar(vars,other._wgtVar?other._wgtVar->GetName():0)),
  _nRealF(0),
  _nCat(0),
  _nEntries(0),	   
  _firstRealF(0),
  _firstCat(0),
  _sumWeight(0),
  _sumWeightCarry(0),
  _extWgtArray(0),
  _extWgtErrLoArray(0),
  _extWgtErrHiArray(0),
  _extSumW2Array(0),
  _curWgt(1),
  _curWgtErrLo(0),
  _curWgtErrHi(0),
  _curWgtErr(0),
  _cache(0),
  _cacheOwner(0),
  _forcedUpdate(kFALSE)
{
  TIterator* iter = _varsww.createIterator() ;
  RooAbsArg* arg ;
  while((arg=(RooAbsArg*)iter->Next())) {
    arg->attachToVStore(*this) ;
  }
  delete iter ;

  setAllBuffersNative() ;
  
  // now copy contents of tree storage here
  reserve(other.numEntries());
  for (Int_t i=0 ; i<other.numEntries() ; i++) {
    other.get(i) ;
    _varsww = other._varsww ;
    fill() ;
  }
  TRACE_CREATE
  
}


////////////////////////////////////////////////////////////////////////////////
/// Clone ctor, must connect internal storage to given new external set of vars

RooVectorDataStore::RooVectorDataStore(const RooVectorDataStore& other, const RooArgSet& vars, const char* newname) :
  RooAbsDataStore(other,varsNoWeight(vars,other._wgtVar?other._wgtVar->GetName():0),newname),
  _varsww(vars),
  _wgtVar(other._wgtVar?weightVar(vars,other._wgtVar->GetName()):0),
  _nRealF(0),
  _nCat(0),
  _nEntries(other._nEntries),	 
  _sumWeight(other._sumWeight),
  _sumWeightCarry(other._sumWeightCarry),
  _extWgtArray(other._extWgtArray),
  _extWgtErrLoArray(other._extWgtErrLoArray),
  _extWgtErrHiArray(other._extWgtErrHiArray),
  _extSumW2Array(other._extSumW2Array),
  _curWgt(other._curWgt),
  _curWgtErrLo(other._curWgtErrLo),
  _curWgtErrHi(other._curWgtErrHi),
  _curWgtErr(other._curWgtErr),
  _cache(0),
  _forcedUpdate(kFALSE)
{
  for (const auto realVec : other._realStoreList) {
    auto real = static_cast<RooAbsReal*>(vars.find(realVec->bufArg()->GetName()));
    if (real) {
      // Clone vector
      _realStoreList.push_back(new RealVector(*realVec, real)) ;
      // Adjust buffer pointer
      real->attachToVStore(*this) ;
    }
  }
  
  vector<RealFullVector*>::const_iterator fiter = other._realfStoreList.begin() ;
  for (; fiter!=other._realfStoreList.end() ; ++fiter) {
    RooAbsReal* real = (RooAbsReal*) vars.find((*fiter)->bufArg()->GetName()) ;
    if (real) {
      // Clone vector
      _realfStoreList.push_back(new RealFullVector(**fiter,real)) ;
      // Adjust buffer pointer
      real->attachToVStore(*this) ;
      _nRealF++ ;
    }
  }

  vector<CatVector*>::const_iterator citer = other._catStoreList.begin() ;
  for (; citer!=other._catStoreList.end() ; ++citer) {
    RooAbsCategory* cat = (RooAbsCategory*) vars.find((*citer)->bufArg()->GetName()) ;
    if (cat) {
      // Clone vector
      _catStoreList.push_back(new CatVector(**citer,cat)) ;
      // Adjust buffer pointer
      cat->attachToVStore(*this) ;
      _nCat++ ;
    }
  }

  setAllBuffersNative() ;

  _firstRealF = _realfStoreList.size()>0 ? &_realfStoreList.front() : 0 ;
  _firstCat = _catStoreList.size()>0 ? &_catStoreList.front() : 0 ;
  TRACE_CREATE

}





////////////////////////////////////////////////////////////////////////////////

RooVectorDataStore::RooVectorDataStore(const char *name, const char *title, RooAbsDataStore& tds, 
			 const RooArgSet& vars, const RooFormulaVar* cutVar, const char* cutRange,
			 Int_t nStart, Int_t nStop, Bool_t /*copyCache*/, const char* wgtVarName) :

  RooAbsDataStore(name,title,varsNoWeight(vars,wgtVarName)),
  _varsww(vars),
  _wgtVar(weightVar(vars,wgtVarName)),
  _nRealF(0),
  _nCat(0),
  _nEntries(0),	   
  _firstRealF(0),
  _firstCat(0),
  _sumWeight(0),
  _sumWeightCarry(0),
  _extWgtArray(0),
  _extWgtErrLoArray(0),
  _extWgtErrHiArray(0),
  _extSumW2Array(0),
  _curWgt(1),
  _curWgtErrLo(0),
  _curWgtErrHi(0),
  _curWgtErr(0),
  _cache(0),
  _forcedUpdate(kFALSE)
{
  TIterator* iter = _varsww.createIterator() ;
  RooAbsArg* arg ;
  while((arg=(RooAbsArg*)iter->Next())) {
    arg->attachToVStore(*this) ;
  }
  delete iter ;

  setAllBuffersNative() ;

  // Deep clone cutVar and attach clone to this dataset
  RooFormulaVar* cloneVar = 0;
  if (cutVar) {    
    cloneVar = (RooFormulaVar*) cutVar->cloneTree() ;
    cloneVar->attachDataStore(tds) ;
  }

  RooVectorDataStore* vds = dynamic_cast<RooVectorDataStore*>(&tds) ;
  if (vds && vds->_cache) {    
    _cache = new RooVectorDataStore(*vds->_cache) ;
  }
  
  loadValues(&tds,cloneVar,cutRange,nStart,nStop);

  delete cloneVar ;
  TRACE_CREATE
}






////////////////////////////////////////////////////////////////////////////////
/// Destructor

RooVectorDataStore::~RooVectorDataStore()
{
  for (auto elm : _realStoreList) {
    delete elm;
  }

  for (auto elm : _realfStoreList) {
    delete elm;
  }

  for (auto elm : _catStoreList) {
    delete elm;
  }

  delete _cache ;
  TRACE_DESTROY
}




////////////////////////////////////////////////////////////////////////////////
/// Return true if currently loaded coordinate is considered valid within
/// the current range definitions of all observables

Bool_t RooVectorDataStore::valid() const 
{
  return kTRUE ;
}




////////////////////////////////////////////////////////////////////////////////
/// Interface function to TTree::Fill

Int_t RooVectorDataStore::fill()
{
  for (auto realVec : _realStoreList) {
    realVec->fill() ;
  }

  for (auto fullVec : _realfStoreList) {
    fullVec->fill() ;
  }

  for (auto catVec : _catStoreList) {
    catVec->fill() ;
  }
  // use Kahan's algorithm to sum up weights to avoid loss of precision
  Double_t y = (_wgtVar ? _wgtVar->getVal() : 1.) - _sumWeightCarry;
  Double_t t = _sumWeight + y;
  _sumWeightCarry = (t - _sumWeight) - y;
  _sumWeight = t;
  _nEntries++ ;  

  return 0 ;
}
 


////////////////////////////////////////////////////////////////////////////////
/// Load the n-th data point (n='index') in memory
/// and return a pointer to the internal RooArgSet
/// holding its coordinates.

const RooArgSet* RooVectorDataStore::get(Int_t index) const 
{
  if (index>=_nEntries) return 0 ;
    
  for (const auto realV : _realStoreList) {
    realV->get(index);
  }

  if (_nRealF>0) {
    for (Int_t i=0 ; i<_nRealF ; i++) {
      (*(_firstRealF+i))->get(index) ;
    }
  }

  if (_nCat>0) {
    for (Int_t i=0 ; i<_nCat ; i++) {
      (*(_firstCat+i))->get(index) ;
    }
  }

  if (_doDirtyProp) {
    // Raise all dirty flags 
    for (auto var : _vars) {
      var->setValueDirty(); // This triggers recalculation of all clients
    }     
  }
  
  // Update current weight cache
  if (_extWgtArray) {

    // If external array is specified use that  
    _curWgt = _extWgtArray[index] ;
    _curWgtErrLo = _extWgtErrLoArray[index] ;
    _curWgtErrHi = _extWgtErrHiArray[index] ;
    _curWgtErr   = sqrt(_extSumW2Array[index]) ;

  } else if (_wgtVar) {

    // Otherwise look for weight variable
    _curWgt = _wgtVar->getVal() ;
    _curWgtErrLo = _wgtVar->getAsymErrorLo() ;
    _curWgtErrHi = _wgtVar->getAsymErrorHi() ;
    _curWgtErr   = _wgtVar->hasAsymError() ? ((_wgtVar->getAsymErrorHi() - _wgtVar->getAsymErrorLo())/2)  : _wgtVar->getError() ;

  } // else {

//     // Otherwise return 1 
//     _curWgt=1.0 ;
//     _curWgtErrLo = 0 ;
//     _curWgtErrHi = 0 ;
//     _curWgtErr = 0 ;
    
//   }

  if (_cache) {
    _cache->get(index) ;
  }

  return &_vars;
}


////////////////////////////////////////////////////////////////////////////////
/// Load the n-th data point (n='index') in memory
/// and return a pointer to the internal RooArgSet
/// holding its coordinates.

const RooArgSet* RooVectorDataStore::getNative(Int_t index) const 
{
  if (index>=_nEntries) return 0 ;
    
  for (const auto realV : _realStoreList) {
    realV->getNative(index) ;
  }

  if (_nRealF>0) {
    for (Int_t i=0 ; i<_nRealF ; i++) {
      (*(_firstRealF+i))->getNative(index) ;
    }
  }

  if (_nCat>0) {
    for (Int_t i=0 ; i<_nCat ; i++) {
      (*(_firstCat+i))->getNative(index) ;
    }
  }

  if (_doDirtyProp) {
    // Raise all dirty flags 
    for (auto var : _vars) {
      var->setValueDirty() ; // This triggers recalculation of all clients
    }     
  }
  
  // Update current weight cache
  if (_extWgtArray) {

    // If external array is specified use that  
    _curWgt = _extWgtArray[index] ;
    _curWgtErrLo = _extWgtErrLoArray[index] ;
    _curWgtErrHi = _extWgtErrHiArray[index] ;
    _curWgtErr   = sqrt(_extSumW2Array[index]) ;

  } else if (_wgtVar) {

    // Otherwise look for weight variable
    _curWgt = _wgtVar->getVal() ;
    _curWgtErrLo = _wgtVar->getAsymErrorLo() ;
    _curWgtErrHi = _wgtVar->getAsymErrorHi() ;
    _curWgtErr   = _wgtVar->hasAsymError() ? ((_wgtVar->getAsymErrorHi() - _wgtVar->getAsymErrorLo())/2)  : _wgtVar->getError() ;

  } else {

    // Otherwise return 1 
    _curWgt=1.0 ;
    _curWgtErrLo = 0 ;
    _curWgtErrHi = 0 ;
    _curWgtErr = 0 ;
    
  }

  if (_cache) {
    _cache->getNative(index) ;
  }

  return &_vars;
}



////////////////////////////////////////////////////////////////////////////////
/// Return the weight of the n-th data point (n='index') in memory

Double_t RooVectorDataStore::weight(Int_t index) const 
{
  get(index) ;
  return weight() ;
}



////////////////////////////////////////////////////////////////////////////////
/// Return the weight of the n-th data point (n='index') in memory

Double_t RooVectorDataStore::weight() const 
{
  return _curWgt ;
}



////////////////////////////////////////////////////////////////////////////////
/// Return the error of the current weight.
/// @param[in] etype Switch between simple Poisson or sum-of-weights statistics

Double_t RooVectorDataStore::weightError(RooAbsData::ErrorType etype) const 
{
  if (_extWgtArray) {

    // We have a weight array, use that info

    // Return symmetric error on current bin calculated either from Poisson statistics or from SumOfWeights
    Double_t lo,hi ;
    weightError(lo,hi,etype) ;
    return (lo+hi)/2 ;

   } else if (_wgtVar) {

    // We have a a weight variable, use that info
    if (_wgtVar->hasAsymError()) {
      return ( _wgtVar->getAsymErrorHi() - _wgtVar->getAsymErrorLo() ) / 2 ;
    } else if (_wgtVar->hasError(kFALSE)) {
      return _wgtVar->getError() ;    
    } else {
      return 0 ;
    }

  } else {

    // We have no weights
    return 0 ;

  }
}



////////////////////////////////////////////////////////////////////////////////

void RooVectorDataStore::weightError(Double_t& lo, Double_t& hi, RooAbsData::ErrorType etype) const
{
  if (_extWgtArray) {
    
    // We have a weight array, use that info
    switch (etype) {
      
    case RooAbsData::Auto:
      throw string(Form("RooDataHist::weightError(%s) error type Auto not allowed here",GetName())) ;
      break ;
      
    case RooAbsData::Expected:
      throw string(Form("RooDataHist::weightError(%s) error type Expected not allowed here",GetName())) ;
      break ;
      
    case RooAbsData::Poisson:
      // Weight may be preset or precalculated    
      if (_curWgtErrLo>=0) {
	lo = _curWgtErrLo ;
	hi = _curWgtErrHi ;
	return ;
      }
      
      // Otherwise Calculate poisson errors
      Double_t ym,yp ;  
      RooHistError::instance().getPoissonInterval(Int_t(weight()+0.5),ym,yp,1) ;
      lo = weight()-ym ;
      hi = yp-weight() ;
      return ;
      
    case RooAbsData::SumW2:
      lo = _curWgtErr ;
      hi = _curWgtErr ;
      return ;
      
    case RooAbsData::None:
      lo = 0 ;
      hi = 0 ;
      return ;
    }    
    
  } else if (_wgtVar) {

    // We have a a weight variable, use that info
    if (_wgtVar->hasAsymError()) {
      hi = _wgtVar->getAsymErrorHi() ;
      lo = _wgtVar->getAsymErrorLo() ;
    } else {
      hi = _wgtVar->getError() ;
      lo = _wgtVar->getError() ;
    }  

  } else {

    // We are unweighted
    lo=0 ;
    hi=0 ;

  }
}



////////////////////////////////////////////////////////////////////////////////
///

void RooVectorDataStore::loadValues(const RooAbsDataStore *ads, const RooFormulaVar* select, const char* rangeName, Int_t nStart, Int_t nStop) 
{
  // Load values from dataset 't' into this data collection, optionally
  // selecting events using 'select' RooFormulaVar
  //

  // Redirect formula servers to source data row
  RooFormulaVar* selectClone(0) ;
  if (select) {
    selectClone = (RooFormulaVar*) select->cloneTree() ;
    selectClone->recursiveRedirectServers(*ads->get()) ;
    selectClone->setOperMode(RooAbsArg::ADirty,kTRUE) ;
  }

  // Force DS internal initialization
  ads->get(0) ;

  // Loop over events in source tree   
  RooAbsArg* arg = 0;
  TIterator* destIter = _varsww.createIterator() ;
  Int_t nevent = nStop < ads->numEntries() ? nStop : ads->numEntries() ;
  Bool_t allValid ;

  Bool_t isTDS = dynamic_cast<const RooTreeDataStore*>(ads) ;
  Bool_t isVDS = dynamic_cast<const RooVectorDataStore*>(ads) ;

  // Check if weight is being renamed - if so set flag to enable special handling in copy loop
  Bool_t weightRename(kFALSE) ;
  Bool_t newWeightVar = _wgtVar ? _wgtVar->getAttribute("NewWeight") : kFALSE ;

  if (_wgtVar && isVDS && ((RooVectorDataStore*)(ads))->_wgtVar) {
    if (string(_wgtVar->GetName())!=((RooVectorDataStore*)(ads))->_wgtVar->GetName() && !newWeightVar) {
      weightRename=kTRUE ;
    }
  }
  if (_wgtVar && isTDS && ((RooTreeDataStore*)(ads))->_wgtVar) {
    if (string(_wgtVar->GetName())!=((RooTreeDataStore*)(ads))->_wgtVar->GetName() && !newWeightVar) {
      weightRename=kTRUE ;
    }
  }

  reserve(numEntries() + (nevent - nStart));
  for(Int_t i=nStart; i < nevent ; ++i) {
    ads->get(i) ;
    
    // Does this event pass the cuts?
    if (selectClone && selectClone->getVal()==0) {
      continue ; 
    }


    if (isTDS) {
      _varsww.assignValueOnly(((RooTreeDataStore*)ads)->_varsww) ;
      if (weightRename) {
	_wgtVar->setVal(((RooTreeDataStore*)ads)->_wgtVar->getVal()) ;
      }
    } else if (isVDS) {
      _varsww.assignValueOnly(((RooVectorDataStore*)ads)->_varsww) ;
      if (weightRename) {
	_wgtVar->setVal(((RooVectorDataStore*)ads)->_wgtVar->getVal()) ;
      }
    } else {
      _varsww.assignValueOnly(*ads->get()) ;
    }

    destIter->Reset() ;
    // Check that all copied values are valid
    allValid=kTRUE ;
    while((arg=(RooAbsArg*)destIter->Next())) {
      if (!arg->isValid() || (rangeName && !arg->inRange(rangeName))) {
	allValid=kFALSE ;
	break ;
      }
    }
    if (!allValid) {
      continue ;
    }
    
    //_cachedVars = ((RooTreeDataStore*)ads)->_cachedVars ;
    fill() ;
   }

  delete destIter ;  
  delete selectClone ;
  
  SetTitle(ads->GetTitle());
}





////////////////////////////////////////////////////////////////////////////////

Bool_t RooVectorDataStore::changeObservableName(const char* /*from*/, const char* /*to*/) 
{
  return kFALSE ;
}

  

////////////////////////////////////////////////////////////////////////////////
/// Add a new column to the data set which holds the pre-calculated values
/// of 'newVar'. This operation is only meaningful if 'newVar' is a derived
/// value.
///
/// The return value points to the added element holding 'newVar's value
/// in the data collection. The element is always the corresponding fundamental
/// type of 'newVar' (e.g. a RooRealVar if 'newVar' is a RooFormulaVar)
///
/// Note: This function is explicitly NOT intended as a speed optimization
///       opportunity for the user. Components of complex PDFs that can be
///       precalculated with the dataset are automatically identified as such
///       and will be precalculated when fitting to a dataset
/// 
///       By forcibly precalculating functions with non-trivial Jacobians,
///       or functions of multiple variables occurring in the data set,
///       using addColumn(), you may alter the outcome of the fit. 
///
///       Only in cases where such a modification of fit behaviour is intentional, 
///       this function should be used. 

RooAbsArg* RooVectorDataStore::addColumn(RooAbsArg& newVar, Bool_t /*adjustRange*/)
{
  // Create a fundamental object of the right type to hold newVar values
  RooAbsArg* valHolder= newVar.createFundamental();
  // Sanity check that the holder really is fundamental
  if(!valHolder->isFundamental()) {
    coutE(InputArguments) << GetName() << "::addColumn: holder argument is not fundamental: \""
	 << valHolder->GetName() << "\"" << endl;
    return 0;
  }

  // Clone variable and attach to cloned tree 
  RooAbsArg* newVarClone = newVar.cloneTree() ;
  newVarClone->recursiveRedirectServers(_vars,kFALSE) ;

  // Attach value place holder to this tree
  valHolder->attachToVStore(*this) ;
  _vars.add(*valHolder) ;
  _varsww.add(*valHolder) ;

  // Fill values of of placeholder
  RealVector* rv(0) ;
  CatVector* cv(0) ;
  if (dynamic_cast<RooAbsReal*>(valHolder)) {
    rv = addReal((RooAbsReal*)valHolder); 
    rv->resize(numEntries()) ;
  } else if (dynamic_cast<RooAbsCategory*>((RooAbsCategory*)valHolder)) {
    cv = addCategory((RooAbsCategory*)valHolder) ;
    cv->resize(numEntries()) ;
  } 

  for (int i=0 ; i<numEntries() ; i++) {
    getNative(i) ;

    newVarClone->syncCache(&_vars) ;
    valHolder->copyCache(newVarClone) ;

    if (rv) rv->write(i) ;
    if (cv) cv->write(i) ;
  }

  delete newVarClone ;  
  return valHolder ;

}



////////////////////////////////////////////////////////////////////////////////
/// Utility function to add multiple columns in one call
/// See addColumn() for details

RooArgSet* RooVectorDataStore::addColumns(const RooArgList& varList)
{
  TIterator* vIter = varList.createIterator() ;
  RooAbsArg* var ;

  checkInit() ;

  TList cloneSetList ;
  RooArgSet cloneSet ;
  RooArgSet* holderSet = new RooArgSet ;

  while((var=(RooAbsArg*)vIter->Next())) {
    // Create a fundamental object of the right type to hold newVar values
    RooAbsArg* valHolder= var->createFundamental();
    holderSet->add(*valHolder) ;

    // Sanity check that the holder really is fundamental
    if(!valHolder->isFundamental()) {
      coutE(InputArguments) << GetName() << "::addColumn: holder argument is not fundamental: \""
	   << valHolder->GetName() << "\"" << endl;
      return 0;
    }
    
    // Clone variable and attach to cloned tree 
    RooArgSet* newVarCloneList = (RooArgSet*) RooArgSet(*var).snapshot() ;  
    if (!newVarCloneList) {
      coutE(InputArguments) << "RooTreeDataStore::RooTreeData(" << GetName() 
			    << ") Couldn't deep-clone variable " << var->GetName() << ", abort." << endl ;
      return 0 ;
    }
    RooAbsArg* newVarClone = newVarCloneList->find(var->GetName()) ;   
    newVarClone->recursiveRedirectServers(_vars,kFALSE) ;
    newVarClone->recursiveRedirectServers(*holderSet,kFALSE) ;

    cloneSetList.Add(newVarCloneList) ;
    cloneSet.add(*newVarClone) ;

    // Attach value place holder to this tree
    valHolder->attachToVStore(*this) ;
    _vars.add(*valHolder) ;
  }
  delete vIter ;


  TIterator* cIter = cloneSet.createIterator() ;
  TIterator* hIter = holderSet->createIterator() ;
  RooAbsArg *cloneArg, *holder ;

  // Dimension storage area for new vectors
  while((holder = (RooAbsArg*)hIter->Next())) {
      if (dynamic_cast<RooAbsReal*>(holder)) {
	addReal((RooAbsReal*)holder)->resize(numEntries()) ;
      } else {
	addCategory((RooAbsCategory*)holder)->resize(numEntries()) ;
      }
    }

  // Fill values of of placeholder
  for (int i=0 ; i<numEntries() ; i++) {
    getNative(i) ;

    cIter->Reset() ;
    hIter->Reset() ;
    while((cloneArg=(RooAbsArg*)cIter->Next())) {
      holder = (RooAbsArg*)hIter->Next() ;

      cloneArg->syncCache(&_vars) ;

      holder->copyCache(cloneArg) ;

      if (dynamic_cast<RooAbsReal*>(holder)) {
	addReal((RooAbsReal*)holder)->write(i) ;
      } else {
	addCategory((RooAbsCategory*)holder)->write(i) ;
      }
    }
  }
  
  delete cIter ;
  delete hIter ;

  cloneSetList.Delete() ;
  return holderSet ;
}




////////////////////////////////////////////////////////////////////////////////
/// Merge columns of supplied data set(s) with this data set.  All
/// data sets must have equal number of entries.  In case of
/// duplicate columns the column of the last dataset in the list
/// prevails

RooAbsDataStore* RooVectorDataStore::merge(const RooArgSet& allVars, list<RooAbsDataStore*> dstoreList)
{
  RooVectorDataStore* mergedStore = new RooVectorDataStore("merged","merged",allVars) ;

  Int_t nevt = dstoreList.front()->numEntries() ;
  mergedStore->reserve(nevt);
  for (int i=0 ; i<nevt ; i++) {

    // Copy data from self
    mergedStore->_vars = *get(i) ;
      
    // Copy variables from merge sets
    for (list<RooAbsDataStore*>::iterator iter = dstoreList.begin() ; iter!=dstoreList.end() ; ++iter) {
      const RooArgSet* partSet = (*iter)->get(i) ;
      mergedStore->_vars = *partSet ;
    }

    mergedStore->fill() ;
  }
  return mergedStore ;
}



void RooVectorDataStore::reserve(Int_t nEvts)
{
  for (auto elm : _realStoreList) {
    elm->reserve(nEvts);
  }

  for (auto elm : _realfStoreList) {
    elm->reserve(nEvts);
  }

  for (auto elm : _catStoreList) {
    elm->reserve(nEvts);
  }
}

////////////////////////////////////////////////////////////////////////////////

void RooVectorDataStore::append(RooAbsDataStore& other) 
{
  Int_t nevt = other.numEntries() ;
  reserve(nevt + numEntries());
  for (int i=0 ; i<nevt ; i++) {  
    _vars = *other.get(i) ;
    if (_wgtVar) {
      _wgtVar->setVal(other.weight()) ;
    }
    
    fill() ;
  }
}



////////////////////////////////////////////////////////////////////////////////

Int_t RooVectorDataStore::numEntries() const 
{
  return _nEntries ;
}



////////////////////////////////////////////////////////////////////////////////

void RooVectorDataStore::reset() 
{
  _nEntries=0 ;
  _sumWeight=_sumWeightCarry=0 ;

  for (auto elm : _realStoreList) {
    elm->reset() ;
  }  

  for (auto elm : _realfStoreList) {
    elm->reset() ;
  }

  for (auto elm : _catStoreList) {
    elm->reset() ;
  }  

}

////////////////////////////////////////////////////////////////////////////////
/// Cache given RooAbsArgs: The tree is
/// given direct write access of the args internal cache
/// the args values is pre-calculated for all data points
/// in this data collection. Upon a get() call, the
/// internal cache of 'newVar' will be loaded with the
/// precalculated value and it's dirty flag will be cleared.

void RooVectorDataStore::cacheArgs(const RooAbsArg* owner, RooArgSet& newVarSet, const RooArgSet* nset, Bool_t skipZeroWeights) 
{
  // Delete previous cache, if any
  delete _cache ;
  _cache = 0 ;

  // Reorder cached elements. First constant nodes, then tracked nodes in order of dependence

  // Step 1 - split in constant and tracked
  RooArgSet newVarSetCopy(newVarSet);
  RooArgSet orderedArgs;
  vector<RooAbsArg*> trackArgs;
  for (const auto arg : newVarSetCopy) {
    if (arg->getAttribute("ConstantExpression") && !arg->getAttribute("NOCacheAndTrack")) {
      orderedArgs.add(*arg) ;
    } else {

      // Explicitly check that arg depends on any of the observables, if this
      // is not the case, skip it, as inclusion would result in repeated
      // calculation of a function that has the same value for every event
      // in the likelihood
      if (arg->dependsOn(_vars) && !arg->getAttribute("NOCacheAndTrack")) {
        trackArgs.push_back(arg) ;
      } else {
        newVarSet.remove(*arg) ;
      }
    }
  }

  // Step 2 - reorder tracked nodes
  std::sort(trackArgs.begin(), trackArgs.end(), [](RooAbsArg* left, RooAbsArg* right){
    return right->dependsOn(*left);
  });

  // Step 3 - put back together
  for (const auto trackedArg : trackArgs) {
    orderedArgs.add(*trackedArg);
  }
  
  // WVE need to prune tracking entries _below_ constant nodes as the're not needed
//   cout << "Number of Cache-and-Tracked args are " << trackArgs.size() << endl ;
//   cout << "Compound ordered cache parameters = " << endl ;
//   orderedArgs.Print("v") ;
  
  checkInit() ;
  
  std::vector<RooArgSet*> vlist;
  RooArgList cloneSet;

  for (const auto var : orderedArgs) {

    // Clone variable and attach to cloned tree 
    RooArgSet* newVarCloneList = (RooArgSet*) RooArgSet(*var).snapshot() ;  
    RooAbsArg* newVarClone = newVarCloneList->find(var->GetName()) ;   
    newVarClone->recursiveRedirectServers(_vars,kFALSE) ;

    vlist.push_back(newVarCloneList) ;
    cloneSet.add(*newVarClone) ;
  }

  _cacheOwner = (RooAbsArg*) owner ;
  RooVectorDataStore* newCache = new RooVectorDataStore("cache","cache",orderedArgs) ;


  RooAbsArg::setDirtyInhibit(kTRUE) ;

  std::vector<RooArgSet*> nsetList ;
  std::vector<RooArgSet*> argObsList ;

  // Now need to attach branch buffers of clones
  RooArgSet *anset(0), *acset(0) ;
  for (const auto arg : cloneSet) {
    arg->attachToVStore(*newCache) ;
    
    RooArgSet* argObs = nset ? arg->getObservables(*nset) : arg->getVariables() ;
    argObsList.push_back(argObs) ;
    
    RooArgSet* normSet(0) ;
    const char* catNset = arg->getStringAttribute("CATNormSet") ;
    if (catNset) {
//       cout << "RooVectorDataStore::cacheArgs() cached node " << arg->GetName() << " has a normalization set specification CATNormSet = " << catNset << endl ;
      RooNameSet rns ;
      rns.setNameList(catNset) ;
      anset = rns.select(nset?*nset:RooArgSet()) ;
      normSet = (RooArgSet*) anset->selectCommon(*argObs) ;
      
    }
    const char* catCset = arg->getStringAttribute("CATCondSet") ;
    if (catCset) {
//       cout << "RooVectorDataStore::cacheArgs() cached node " << arg->GetName() << " has a conditional observable set specification CATCondSet = " << catCset << endl ;
      RooNameSet rns ;
      rns.setNameList(catCset) ;
      acset = rns.select(nset?*nset:RooArgSet()) ;
      
      argObs->remove(*acset,kTRUE,kTRUE) ;
      normSet = argObs ;
    }

    // now construct normalization set for component from cset/nset spec
//     if (normSet) {
//       cout << "RooVectorDaraStore::cacheArgs() component " << arg->GetName() << " has custom normalization set " << *normSet << endl ;
//     }
    nsetList.push_back(normSet) ;    
  }


  // Fill values of of placeholder
  newCache->reserve(numEntries());
  for (int i=0 ; i<numEntries() ; i++) {
    getNative(i) ;
    if (weight()!=0 || !skipZeroWeights) {    
      for (std::size_t j = 0; j < cloneSet.size(); ++j) {
        auto& cloneArg = cloneSet[j];
        auto argNSet = nsetList[j];
        // WVE need to intervene here for condobs from ProdPdf
        cloneArg.syncCache(argNSet ? argNSet : nset) ;
      }
    }
    newCache->fill() ;
  }

  RooAbsArg::setDirtyInhibit(kFALSE) ;


  // Now need to attach branch buffers of original function objects 
  for (const auto arg : orderedArgs) {
    arg->attachToVStore(*newCache) ;
    
    // Activate change tracking mode, if requested
    if (!arg->getAttribute("ConstantExpression") && dynamic_cast<RooAbsReal*>(arg)) {
      RealVector* rv = newCache->addReal((RooAbsReal*)arg) ;      
      RooArgSet* deps = arg->getParameters(_vars) ;
      rv->setDependents(*deps) ;

      // WV lookup normalization set and associate with RealVector
      // find ordinal number of arg in original list 
      Int_t idx = cloneSet.index(arg->GetName()) ;

      coutI(Optimization) << "RooVectorDataStore::cacheArg() element " << arg->GetName() << " has change tracking enabled on parameters " << *deps << endl ;
      rv->setNset(nsetList[idx]) ;
      delete deps ;
    }
    
  }


  for (auto set : vlist) {
    delete set;
  }  
  for (auto set : argObsList) {
    delete set;
  }  

  _cache = newCache ;
  _cache->setDirtyProp(_doDirtyProp) ;
}


void RooVectorDataStore::forceCacheUpdate()
{
  if (_cache) _forcedUpdate = kTRUE ; 
}



////////////////////////////////////////////////////////////////////////////////

void RooVectorDataStore::recalculateCache( const RooArgSet *projectedArgs, Int_t firstEvent, Int_t lastEvent, Int_t stepSize, Bool_t skipZeroWeights) 
{
  if (!_cache) return ;

  std::vector<RooVectorDataStore::RealVector *> tv;
  tv.reserve(static_cast<std::size_t>(_cache->_realStoreList.size() * 0.7)); // Typically, 30..60% need to be recalculated

  // Check which items need recalculation
  for (const auto realVec : _cache->_realStoreList) {
    if (_forcedUpdate || realVec->needRecalc()) {
       tv.push_back(realVec);
       realVec->_nativeReal->setOperMode(RooAbsArg::ADirty);
       realVec->_nativeReal->_operMode = RooAbsArg::Auto;
    }    
  }
  _forcedUpdate = kFALSE ;

  // If no recalculations are needed stop here
  if (tv.empty()) {
     return;
  }


  // Refill caches of elements that require recalculation
  RooArgSet* ownedNset = 0 ;
  RooArgSet* usedNset = 0 ;
  if (projectedArgs && projectedArgs->getSize()>0) {
    ownedNset = (RooArgSet*) _vars.snapshot(kFALSE) ;
    ownedNset->remove(*projectedArgs,kFALSE,kTRUE);
    usedNset = ownedNset ;
  } else {
    usedNset = &_vars ;
  }


  for (int i=firstEvent ; i<lastEvent ; i+=stepSize) {
    get(i) ;    
    Bool_t zeroWeight = (weight()==0) ;
    if (!zeroWeight || !skipZeroWeights) {
       for (auto realVector : tv) {
          realVector->_nativeReal->_valueDirty = kTRUE;
          realVector->_nativeReal->getValV(realVector->_nset ? realVector->_nset : usedNset);
          realVector->write(i);
      }
    }
  }  
  
  for (auto realVector : tv) {
     realVector->_nativeReal->setOperMode(RooAbsArg::AClean);
  }  

  delete ownedNset ;

}


////////////////////////////////////////////////////////////////////////////////
/// Initialize cache of dataset: attach variables of cache ArgSet
/// to the corresponding TTree branches

void RooVectorDataStore::attachCache(const RooAbsArg* newOwner, const RooArgSet& cachedVarsIn) 
{
  // Only applicable if a cache exists
  if (!_cache) return ;

  // Clone ctor, must connect internal storage to given new external set of vars
  std::vector<RealVector*> cacheElements(_cache->realStoreList());
  cacheElements.insert(cacheElements.end(), _cache->_realfStoreList.begin(), _cache->_realfStoreList.end());

  for (const auto elm : cacheElements) {
    auto real = static_cast<RooAbsReal*>(cachedVarsIn.find(elm->bufArg()->GetName()));
    if (real) {
      // Adjust buffer pointer
      real->attachToVStore(*_cache) ;
    }
  }

  for (const auto catVec : _cache->_catStoreList) {
    auto cat = static_cast<RooAbsCategory*>(cachedVarsIn.find(catVec->bufArg()->GetName()));
    if (cat) {
      // Adjust buffer pointer
      cat->attachToVStore(*_cache) ;
    }
  }

  _cacheOwner = (RooAbsArg*) newOwner ;
}




////////////////////////////////////////////////////////////////////////////////

void RooVectorDataStore::resetCache() 
{
  delete _cache;
  _cache = nullptr;
  _cacheOwner = nullptr;
  return ;
}





////////////////////////////////////////////////////////////////////////////////
/// Disabling of branches is (intentionally) not implemented in vector
/// data stores (as the doesn't result in a net saving of time)

void RooVectorDataStore::setArgStatus(const RooArgSet& /*set*/, Bool_t /*active*/) 
{
  return ;
}




////////////////////////////////////////////////////////////////////////////////

void RooVectorDataStore::attachBuffers(const RooArgSet& extObs) 
{
  for (auto arg : _varsww) {
    RooAbsArg* extArg = extObs.find(arg->GetName()) ;
    if (extArg) {
      extArg->attachToVStore(*this) ;
    }
  }
}



////////////////////////////////////////////////////////////////////////////////

void RooVectorDataStore::resetBuffers() 
{ 
  for (auto arg : _varsww) {
    arg->attachToVStore(*this);
  }
}  



////////////////////////////////////////////////////////////////////////////////

void RooVectorDataStore::dump()
{
  cout << "RooVectorDataStor::dump()" << endl ;
  
  cout << "_varsww = " << endl ; _varsww.Print("v") ;
  cout << "realVector list is" << endl ;

  for (const auto elm : _realStoreList) {
    cout << "RealVector " << elm << " _nativeReal = " << elm->_nativeReal << " = " << elm->_nativeReal->GetName() << " bufptr = " << elm->_buf  << endl ;
    cout << " values : " ;
    Int_t imax = elm->_vec.size()>10 ? 10 : elm->_vec.size() ;
    for (Int_t i=0 ; i<imax ; i++) {
      cout << elm->_vec[i] << " " ;
    }
    cout << endl ;
  }    

  for (const auto elm : _realfStoreList) {
    cout << "RealFullVector " << elm << " _nativeReal = " << elm->_nativeReal << " = " << elm->_nativeReal->GetName()
	 << " bufptr = " << elm->_buf  << " errbufptr = " << elm->_bufE << endl ;

    cout << " values : " ;
    Int_t imax = elm->_vec.size()>10 ? 10 : elm->_vec.size() ;
    for (Int_t i=0 ; i<imax ; i++) {
      cout << elm->_vec[i] << " " ;
    }
    cout << endl ;
    if (elm->_vecE) {
      cout << " errors : " ;
      for (Int_t i=0 ; i<imax ; i++) {
	cout << (*elm->_vecE)[i] << " " ;
      }
      cout << endl ;

    }    
  }
}


////////////////////////////////////////////////////////////////////////////////
/// Stream an object of class RooVectorDataStore.

void RooVectorDataStore::Streamer(TBuffer &R__b)
{
  if (R__b.IsReading()) {
    R__b.ReadClassBuffer(RooVectorDataStore::Class(),this);

    if (_realfStoreList.size() > 0)
      _firstRealF = &_realfStoreList.front() ;
    if (_catStoreList.size() > 0)
      _firstCat = &_catStoreList.front() ;

    for (auto elm : _realStoreList) {
      RooAbsArg* arg = _varsww.find(elm->_nativeReal->GetName()) ;
      arg->attachToVStore(*this) ;
    }
    for (auto elm : _realfStoreList) {
      RooAbsArg* arg = _varsww.find(elm->_nativeReal->GetName()) ;
      arg->attachToVStore(*this) ;
    }
    for (auto elm : _catStoreList) {
      RooAbsArg* arg = _varsww.find(elm->_cat->GetName()) ;
      arg->attachToVStore(*this) ;
    }

  } else {
    R__b.WriteClassBuffer(RooVectorDataStore::Class(),this);
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Stream an object of class RooVectorDataStore::RealVector.

void RooVectorDataStore::RealVector::Streamer(TBuffer &R__b)
{
   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(RooVectorDataStore::RealVector::Class(),this);
   } else {
      R__b.WriteClassBuffer(RooVectorDataStore::RealVector::Class(),this);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Stream an object of class RooVectorDataStore::RealFullVector.

void RooVectorDataStore::RealFullVector::Streamer(TBuffer &R__b)
{
   if (R__b.IsReading()) {
     R__b.ReadClassBuffer(RooVectorDataStore::RealFullVector::Class(),this);

     // WVE - It seems that ROOT persistence turns null pointers to vectors into pointers to null-sized vectors 
     //       Intervene here to remove those null-sized vectors and replace with null pointers to not break
     //       assumptions made elsewhere in this class
     if (_vecE  && _vecE->empty()) { delete _vecE   ; _vecE = 0 ; }
     if (_vecEL && _vecEL->empty()) { delete _vecEL ; _vecEL = 0 ; }
     if (_vecEH && _vecEH->empty()) { delete _vecEH ; _vecEH = 0 ; }
   } else {
     R__b.WriteClassBuffer(RooVectorDataStore::RealFullVector::Class(),this);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Stream an object of class RooVectorDataStore::CatVector.

void RooVectorDataStore::CatVector::Streamer(TBuffer &R__b)
{
   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(RooVectorDataStore::CatVector::Class(),this);
      _vec0 = _vec.size()>0 ? &_vec.front() : 0 ;
   } else {
      R__b.WriteClassBuffer(RooVectorDataStore::CatVector::Class(),this);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return a batch of the data columns for all events in [firstEvent, lastEvent[.

std::vector<RooSpan<const double>> RooVectorDataStore::getBatch(std::size_t firstEvent, std::size_t lastEvent) const
{
  std::vector<RooSpan<const double>> ret;

  ret.reserve(_realStoreList.size());

  for (const auto realVec : _realStoreList) {
    ret.emplace_back(realVec->getRange(firstEvent, lastEvent));
  }

  if (_cache) {
    ret.reserve(ret.size() + _cache->_realStoreList.size());

    for (const auto realVec : _cache->_realStoreList) {
      ret.emplace_back(realVec->getRange(firstEvent, lastEvent));
    }
  }

  return ret;
}


////////////////////////////////////////////////////////////////////////////////
/// Return the weights of all events in [first, last[.

RooSpan<const double> RooVectorDataStore::getWeightBatch(std::size_t first, std::size_t last) const
{
  if (_extWgtArray) {
    return RooSpan<const double>(_extWgtArray + first, _extWgtArray + last);
  }


  if (_wgtVar) {
    return _wgtVar->getValBatch(first, last);
  }

  //TODO FIXME!
  static double dummyWeight = 1.;
  return RooSpan<const double>(&dummyWeight, 1);
}



