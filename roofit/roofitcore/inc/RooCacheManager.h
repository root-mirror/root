/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id$
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
#ifndef ROO_CACHE_MANAGER
#define ROO_CACHE_MANAGER

#include "Rtypes.h"

#include "RooMsgService.h"
#include "RooNormSetCache.h"
#include "RooAbsReal.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooAbsCache.h"
#include "RooAbsCacheElement.h"
#include "RooNameReg.h"
#include <vector>

class RooNameSet ;


template<class T>
class RooCacheManager : public RooAbsCache {

public:

  RooCacheManager(Int_t maxSize=2) ;
  RooCacheManager(RooAbsArg* owner, Int_t maxSize=2) ;
  RooCacheManager(const RooCacheManager& other, RooAbsArg* owner=0) ;
  virtual ~RooCacheManager() ;

  T* getObj(const RooArgSet* nset, Int_t* sterileIndex=0, const TNamed* isetRangeName=0) {
    // Getter function without integration set
    return getObj(nset,0,sterileIndex,isetRangeName) ;
  }

  Int_t setObj(const RooArgSet* nset, T* obj, const TNamed* isetRangeName=0) {
    // Setter function without integration set
    return setObj(nset,0,obj,isetRangeName) ;
  }

  inline T* getObj(const RooArgSet* nset, const RooArgSet* iset, Int_t* sterileIdx, const char* isetRangeName)  {
    if (_wired) return _object[0] ;
    return getObj(nset,iset,sterileIdx,RooNameReg::ptr(isetRangeName)) ;
  }

  T* getObj(const RooArgSet* nset, const RooArgSet* iset, Int_t* sterileIndex=0, const TNamed* isetRangeName=0) ;
  Int_t setObj(const RooArgSet* nset, const RooArgSet* iset, T* obj, const TNamed* isetRangeName=0) ;

  void reset() ;
  virtual void sterilize() ;

  Int_t lastIndex() const {
    // Return index of slot used in last get or set operation
    return _lastIndex ;
  }
  Int_t cacheSize() const {
    // Return size of cache
    return _size ;
  }

  virtual Bool_t redirectServersHook(const RooAbsCollection& /*newServerList*/, Bool_t /*mustReplaceAll*/,
                 Bool_t /*nameChange*/, Bool_t /*isRecursive*/) {
    // Interface function to intercept server redirects
    return kFALSE ;
  }
  virtual void operModeHook() {
    // Interface function to intercept cache operation mode changes
  }
  virtual void printCompactTreeHook(std::ostream&, const char *) {
    // Interface function to cache add contents to output in tree printing mode
  }

  T* getObjByIndex(Int_t index) const ;
  const RooNameSet* nameSet1ByIndex(Int_t index) const ;
  const RooNameSet* nameSet2ByIndex(Int_t index) const ;

  virtual void insertObjectHook(T&) {
    // Interface function to perform post-insert operations on cached object
  }

  void wireCache() {
    if (_size==0) {
      oocoutI(_owner,Optimization) << "RooCacheManager::wireCache(" << _owner->GetName() << ") no cached elements!" << std::endl ;
    } else if (_size==1) {
      oocoutI(_owner,Optimization) << "RooCacheManager::wireCache(" << _owner->GetName() << ") now wiring cache" << std::endl ;
      _wired=kTRUE ;
    } else if (_size>1) {
      oocoutI(_owner,Optimization) << "RooCacheManager::wireCache(" << _owner->GetName() << ") cache cannot be wired because it contains more than one element" << std::endl ;
    }
  }

protected:

  Int_t _maxSize ;    //! Maximum size
  Int_t _size ;       //! Actual use
  Int_t _lastIndex ;  //! Last slot accessed

  std::vector<RooNormSetCache> _nsetCache ; //! Normalization/Integration set manager
  std::vector<T*> _object ;                 //! Payload
  Bool_t _wired ;               //! In wired mode, there is a single payload which is returned always

  ClassDef(RooCacheManager,2) // Cache Manager class generic objects
} ;


/// Constructor for simple caches without RooAbsArg payload. A cache
/// made with this constructor is not registered with its owner
/// and will not receive information on server redirects and
/// cache operation mode changes.
template<class T>
RooCacheManager<T>::RooCacheManager(Int_t maxSize) : RooAbsCache(0)
{
  _maxSize = maxSize ;
  _nsetCache.resize(_maxSize) ; // = new RooNormSetCache[maxSize] ;
  _object.resize(_maxSize,0) ; // = new T*[maxSize] ;
  _wired = kFALSE ;
}


/// Constructor for simple caches with RooAbsArg derived payload. A cache
/// made with this constructor is registered with its owner
/// and will receive information on server redirects and
/// cache operation mode changes.
template<class T>
RooCacheManager<T>::RooCacheManager(RooAbsArg* owner, Int_t maxSize) : RooAbsCache(owner)
{
  _maxSize = maxSize ;
  _size = 0 ;

  _nsetCache.resize(_maxSize) ; // = new RooNormSetCache[maxSize] ;
  _object.resize(_maxSize,0) ; // = new T*[maxSize] ;
  _wired = kFALSE ;
  _lastIndex = -1 ;

  Int_t i ;
  for (i=0 ; i<_maxSize ; i++) {
    _object[i]=0 ;
  }

}

/// Copy constructor.
template<class T>
RooCacheManager<T>::RooCacheManager(const RooCacheManager& other, RooAbsArg* owner) : RooAbsCache(other,owner)
{
  _maxSize = other._maxSize ;
  _size = other._size ;

  _nsetCache.resize(_maxSize) ; // = new RooNormSetCache[_maxSize] ;
  _object.resize(_maxSize,0) ; // = new T*[_maxSize] ;
  _wired = kFALSE ;
  _lastIndex = -1 ;

  //std::cout << "RooCacheManager:cctor(" << this << ") other = " << &other << " _size=" << _size << " _maxSize = " << _maxSize << std::endl ;

  Int_t i ;
  for (i=0 ; i<other._size ; i++) {
    _nsetCache[i].initialize(other._nsetCache[i]) ;
    _object[i] = 0 ;
  }

  for (i=other._size ; i<_maxSize ; i++) {
    _object[i] = 0 ;
  }
}

  /// Destructor
template<class T>
RooCacheManager<T>::~RooCacheManager()
{
  for (int i=0 ; i<_size ; i++) {
    delete _object[i] ;
  }
}


  /// Clear the cache
template<class T>
void RooCacheManager<T>::reset()
{
  for (int i=0 ; i<_maxSize ; i++) {
    delete _object[i] ;
    _object[i]=0 ;
    _nsetCache[i].clear() ;
  }
  _lastIndex = -1 ;
  _size = 0 ;
}


/// Clear the cache payload but retain slot mapping w.r.t to
/// normalization and integration sets.
template<class T>
void RooCacheManager<T>::sterilize()
{
  Int_t i ;
  for (i=0 ; i<_maxSize ; i++) {
    delete _object[i] ;
    _object[i]=0 ;
  }
}


/// Insert payload object 'obj' in cache indexed on nset,iset and isetRangeName.
template<class T>
Int_t RooCacheManager<T>::setObj(const RooArgSet* nset, const RooArgSet* iset, T* obj, const TNamed* isetRangeName)
{
  // Check if object is already registered
  Int_t sterileIdx(-1) ;
  if (getObj(nset,iset,&sterileIdx,isetRangeName)) {
    delete obj; // important! do not forget to cleanup memory
    return lastIndex() ;
  }


  if (sterileIdx>=0) {
    // Found sterile slot that can should be recycled [ sterileIndex only set if isetRangeName matches ]

    if (sterileIdx>=_maxSize) {
      //cout << "RooCacheManager<T>::setObj()/SI increasing object cache size from " << _maxSize << " to " << sterileIdx+4 << endl ;
      _maxSize = sterileIdx+4;
      _object.resize(_maxSize,0) ;
      _nsetCache.resize(_maxSize) ;
    }


    _object[sterileIdx] = obj ;

    // Allow optional post-processing of object inserted in cache
    insertObjectHook(*obj) ;

    return lastIndex() ;
  }

  if (_size>=_maxSize-1) {
    //cout << "RooCacheManager<T>::setObj() increasing object cache size from " << _maxSize << " to " << _maxSize*2 << endl ;
    _maxSize *=2 ;
    _object.resize(_maxSize,0) ;
    _nsetCache.resize(_maxSize) ;
  }

  //cout << "RooCacheManager::setObj<T>(" << this << ") _size = " << _size << " _maxSize = " << _maxSize << endl ;
  _nsetCache[_size].autoCache(_owner,nset,iset,isetRangeName,kTRUE) ;
  if (_object[_size]) {
    delete _object[_size] ;
  }

  _object[_size] = obj ;
  _size++ ;

  // Allow optional post-processing of object inserted in cache
  insertObjectHook(*obj) ;

  // Unwire cache in case it was wired
  _wired = kFALSE ;

  return _size-1 ;
}


/// Retrieve payload object indexed on nset,uset and isetRangeName
/// If sterileIdx is not null, it is set to the index of the sterile
/// slot in case such a slot is recycled.
template<class T>
T* RooCacheManager<T>::getObj(const RooArgSet* nset, const RooArgSet* iset, Int_t* sterileIdx, const TNamed* isetRangeName)
{
  // Fast-track for wired mode
  if (_wired) {
    if(_object[0]==0 && sterileIdx) *sterileIdx=0 ;
    return _object[0] ;
  }

  Int_t i ;
  for (i=0 ; i<_size ; i++) {
    if (_nsetCache[i].contains(nset,iset,isetRangeName)==kTRUE) {
      _lastIndex = i ;
      if(_object[i]==0 && sterileIdx) *sterileIdx=i ;
      return _object[i] ;
    }
  }

  for (i=0 ; i<_size ; i++) {
    if (_nsetCache[i].autoCache(_owner,nset,iset,isetRangeName,kFALSE)==kFALSE) {
      _lastIndex = i ;
      if(_object[i]==0 && sterileIdx) *sterileIdx=i ;
      return _object[i] ;
    }
  }

  return 0 ;
}


/// Retrieve payload object by slot index.
template<class T>
T* RooCacheManager<T>::getObjByIndex(Int_t index) const
{
  if (index<0||index>=_size) {
    oocoutE(_owner,ObjectHandling) << "RooCacheManager::getNormListByIndex: ERROR index ("
               << index << ") out of range [0," << _size-1 << "]" << std::endl ;
    return 0 ;
  }
  return _object[index] ;
}


/// Retrieve RooNameSet associated with slot at given index.
template<class T>
const RooNameSet* RooCacheManager<T>::nameSet1ByIndex(Int_t index) const
{
  if (index<0||index>=_size) {
    oocoutE(_owner,ObjectHandling) << "RooCacheManager::getNormListByIndex: ERROR index ("
               << index << ") out of range [0," << _size-1 << "]" << std::endl ;
    return 0 ;
  }
  return &_nsetCache[index].nameSet1() ;
}


/// Retrieve RooNameSet associated with slot at given index.
template<class T>
const RooNameSet* RooCacheManager<T>::nameSet2ByIndex(Int_t index) const
{
  if (index<0||index>=_size) {
    oocoutE(_owner,ObjectHandling) << "RooCacheManager::getNormListByIndex: ERROR index ("
               << index << ") out of range [0," << _size-1 << "]" << std::endl ;
    return 0 ;
  }
  return &_nsetCache[index].nameSet2() ;
}


#endif
