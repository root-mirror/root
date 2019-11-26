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
#ifndef ROO_VECTOR_DATA_STORE
#define ROO_VECTOR_DATA_STORE

#include <list>
#include <vector>
#include <string>
#include <algorithm>
#include "RooAbsDataStore.h" 
#include "TString.h"
#include "RooCatType.h"
#include "RooAbsCategory.h"
#include "RooAbsReal.h"
#include "RooChangeTracker.h"

#define VECTOR_BUFFER_SIZE 1024

class RooAbsArg ;
class RooArgList ;
class TTree ;
class RooFormulaVar ;
class RooArgSet ;

class RooVectorDataStore : public RooAbsDataStore {
public:

  RooVectorDataStore() ; 

  // Empty ctor
  RooVectorDataStore(const char* name, const char* title, const RooArgSet& vars, const char* wgtVarName=0) ;
  virtual RooAbsDataStore* clone(const char* newname=0) const override { return new RooVectorDataStore(*this,newname) ; }
  virtual RooAbsDataStore* clone(const RooArgSet& vars, const char* newname=0) const override { return new RooVectorDataStore(*this,vars,newname) ; }

  RooVectorDataStore(const RooVectorDataStore& other, const char* newname=0) ;
  RooVectorDataStore(const RooTreeDataStore& other, const RooArgSet& vars, const char* newname=0) ;
  RooVectorDataStore(const RooVectorDataStore& other, const RooArgSet& vars, const char* newname=0) ;


  RooVectorDataStore(const char *name, const char *title, RooAbsDataStore& tds, 
		     const RooArgSet& vars, const RooFormulaVar* cutVar, const char* cutRange,
		     Int_t nStart, Int_t nStop, Bool_t /*copyCache*/, const char* wgtVarName=0) ;

  virtual ~RooVectorDataStore() ;

private:
  RooArgSet varsNoWeight(const RooArgSet& allVars, const char* wgtName);
  RooRealVar* weightVar(const RooArgSet& allVars, const char* wgtName);

  // reserve storage for nEvt entries
  void reserve(Int_t nEvt);

public:
  // Write current row
  virtual Int_t fill() override;

  // Retrieve a row
  using RooAbsDataStore::get;
  virtual const RooArgSet* get(Int_t index) const override;

  virtual const RooArgSet* getNative(Int_t index) const;

  virtual Double_t weight() const override;
  virtual Double_t weightError(RooAbsData::ErrorType etype=RooAbsData::Poisson) const override;
  virtual void weightError(Double_t& lo, Double_t& hi, RooAbsData::ErrorType etype=RooAbsData::Poisson) const override;
  virtual Double_t weight(Int_t index) const override;
  virtual Bool_t isWeighted() const override { return (_wgtVar!=0||_extWgtArray!=0) ; }

  virtual std::vector<RooSpan<const double>> getBatch(std::size_t first, std::size_t last) const override;
  virtual RooSpan<const double> getWeightBatch(std::size_t first, std::size_t last) const override;

  // Change observable name
  virtual Bool_t changeObservableName(const char* from, const char* to) override;
  
  // Add one or more columns
  virtual RooAbsArg* addColumn(RooAbsArg& var, Bool_t adjustRange=kTRUE) override;
  virtual RooArgSet* addColumns(const RooArgList& varList) override;

  // Merge column-wise
  RooAbsDataStore* merge(const RooArgSet& allvars, std::list<RooAbsDataStore*> dstoreList) override;

  // Add rows 
  virtual void append(RooAbsDataStore& other) override;

  // General & bookkeeping methods
  virtual Bool_t valid() const override;
  virtual Int_t numEntries() const override;
  virtual Double_t sumEntries() const override { return _sumWeight ; }
  virtual void reset() override;

  // Buffer redirection routines used in inside RooAbsOptTestStatistics
  virtual void attachBuffers(const RooArgSet& extObs) override;
  virtual void resetBuffers() override;
  
  
  // Constant term  optimizer interface
  virtual const RooAbsArg* cacheOwner() override { return _cacheOwner ; }
  virtual void cacheArgs(const RooAbsArg* owner, RooArgSet& varSet, const RooArgSet* nset=0, Bool_t skipZeroWeights=kTRUE) override;
  virtual void attachCache(const RooAbsArg* newOwner, const RooArgSet& cachedVars) override;
  virtual void resetCache() override;
  virtual void recalculateCache(const RooArgSet* /*proj*/, Int_t firstEvent, Int_t lastEvent, Int_t stepSize, Bool_t skipZeroWeights) override;

  virtual void setArgStatus(const RooArgSet& set, Bool_t active) override;

  const RooVectorDataStore* cache() const { return _cache ; }

  void loadValues(const RooAbsDataStore *tds, const RooFormulaVar* select=0, const char* rangeName=0, Int_t nStart=0, Int_t nStop=2000000000) override;
  
  void dump() override;

  void setExternalWeightArray(const Double_t* arrayWgt, const Double_t* arrayWgtErrLo,
      const Double_t* arrayWgtErrHi, const Double_t* arraySumW2) override {
    _extWgtArray = arrayWgt ; 
    _extWgtErrLoArray = arrayWgtErrLo ;
    _extWgtErrHiArray = arrayWgtErrHi ;
    _extSumW2Array = arraySumW2 ;
  }

  virtual void setDirtyProp(Bool_t flag) override {
    _doDirtyProp = flag ; 
    if (_cache) {
      _cache->setDirtyProp(flag) ;
    }
  }

  const RooArgSet& row() { return _varsww ; }

  class RealVector {
  public:

    RealVector(UInt_t initialCapacity=(VECTOR_BUFFER_SIZE / sizeof(Double_t))) : 
      _nativeReal(0), _real(0), _buf(0), _nativeBuf(0), _tracker(0), _nset(0) {
      _vec.reserve(initialCapacity);
    }

    RealVector(RooAbsReal* arg, UInt_t initialCapacity=(VECTOR_BUFFER_SIZE / sizeof(Double_t))) : 
      _nativeReal(arg), _real(0), _buf(0), _nativeBuf(0), _tracker(0), _nset(0) {
      _vec.reserve(initialCapacity);
    }

    virtual ~RealVector() {
      delete _tracker;
      if (_nset) delete _nset ;
    }

    RealVector(const RealVector& other, RooAbsReal* real=0) : 
      _vec(other._vec), _nativeReal(real?real:other._nativeReal), _real(real?real:other._real), _buf(other._buf), _nativeBuf(other._nativeBuf), _nset(0) {
      if (other._tracker) {
        _tracker = new RooChangeTracker(Form("track_%s",_nativeReal->GetName()),"tracker",other._tracker->parameters()) ;
      } else {
        _tracker = 0 ;
      }
      if (other._nset) {
        _nset = new RooArgSet(*other._nset) ;
      }
    }

    RealVector& operator=(const RealVector& other) {
      if (&other==this) return *this;
      _nativeReal = other._nativeReal;
      _real = other._real;
      _buf = other._buf;
      _nativeBuf = other._nativeBuf;
      if (other._vec.size() <= _vec.capacity() / 2 && _vec.capacity() > (VECTOR_BUFFER_SIZE / sizeof(Double_t))) {
        std::vector<Double_t> tmp;
        tmp.reserve(std::max(other._vec.size(), VECTOR_BUFFER_SIZE / sizeof(Double_t)));
        tmp.assign(other._vec.begin(), other._vec.end());
        _vec.swap(tmp);
      } else {
        _vec = other._vec;
      }

      return *this;
    }

    void setNset(RooArgSet* newNset) { _nset = newNset ? new RooArgSet(*newNset) : 0 ; }

    RooArgSet* nset() const { return _nset ; }

    void setBufArg(RooAbsReal* arg) { _nativeReal = arg ; }
    const RooAbsReal* bufArg() const { return _nativeReal ; }

    void setBuffer(RooAbsReal* real, Double_t* newBuf) { 
      _real = real ;
      _buf = newBuf ; 
      if (_nativeBuf==0) {
        _nativeBuf=newBuf ;
      }
    }

    void setNativeBuffer(Double_t* newBuf=0) {       
      _nativeBuf = newBuf ? newBuf : _buf ; 
    }

    void setDependents(const RooArgSet& deps) {
      if (_tracker) {
        delete _tracker ;
      }
      _tracker = new RooChangeTracker(Form("track_%s",_nativeReal->GetName()),"tracker",deps) ;
    }

    Bool_t needRecalc() {
      if (!_tracker) return kFALSE ;
      return _tracker->hasChanged(kTRUE) ;
    }

    void fill() { 
      _vec.push_back(*_buf);
    }

    void write(Int_t i) {
      _vec[i] = *_buf ;
    }
    
    void reset() { 
      // make sure the vector releases the underlying memory
      std::vector<Double_t> tmp;
      _vec.swap(tmp);
    }

    inline void get(Int_t idx) const { 
      *_buf = *(_vec.begin() + idx) ;
    }

    RooSpan<const double> getRange(std::size_t first, std::size_t last) const {
      auto beg = std::min(_vec.cbegin() + first, _vec.cend());
      auto end = std::min(_vec.cbegin() + last,  _vec.cend());

      return RooSpan<const double>(beg, end);
    }

    inline void getNative(Int_t idx) const { 
      *_nativeBuf = *(_vec.begin() + idx) ;
    }

    Int_t size() const { return _vec.size() ; }

    void resize(Int_t siz) {
      if (siz < Int_t(_vec.capacity()) / 2 && _vec.capacity() > (VECTOR_BUFFER_SIZE / sizeof(Double_t))) {
        // do an expensive copy, if we save at least a factor 2 in size
        std::vector<Double_t> tmp;
        tmp.reserve(std::max(siz, Int_t(VECTOR_BUFFER_SIZE / sizeof(Double_t))));
        if (!_vec.empty())
          tmp.assign(_vec.begin(), std::min(_vec.end(), _vec.begin() + siz));
        if (Int_t(tmp.size()) != siz)
          tmp.resize(siz);
        _vec.swap(tmp);
      } else {
        _vec.resize(siz);
      }
    }

    void reserve(Int_t siz) {
      _vec.reserve(siz);
    }

    const std::vector<double>& data() const {
      return _vec;
    }

  protected:
    std::vector<double> _vec;

  private:
    friend class RooVectorDataStore ;
    RooAbsReal* _nativeReal ;
    RooAbsReal* _real ;
    Double_t* _buf ; //!
    Double_t* _nativeBuf ; //!
    RooChangeTracker* _tracker ; //
    RooArgSet* _nset ; //! 
    ClassDef(RealVector,1) // STL-vector-based Data Storage class
  } ;
  

  class RealFullVector : public RealVector {
  public:
    RealFullVector(UInt_t initialCapacity=(VECTOR_BUFFER_SIZE / sizeof(Double_t))) : RealVector(initialCapacity),
      _bufE(0), _bufEL(0), _bufEH(0), 
      _nativeBufE(0), _nativeBufEL(0), _nativeBufEH(0), 
      _vecE(0), _vecEL(0), _vecEH(0) { 
    }

    RealFullVector(RooAbsReal* arg, UInt_t initialCapacity=(VECTOR_BUFFER_SIZE / sizeof(Double_t))) : 
      RealVector(arg,initialCapacity), 
      _bufE(0), _bufEL(0), _bufEH(0), 
      _nativeBufE(0), _nativeBufEL(0), _nativeBufEH(0), 
      _vecE(0), _vecEL(0), _vecEH(0) { 
    }

    virtual ~RealFullVector() {
      if (_vecE) delete _vecE ;
      if (_vecEL) delete _vecEL ;
      if (_vecEH) delete _vecEH ;
    }
    
    RealFullVector(const RealFullVector& other, RooAbsReal* real=0) : RealVector(other,real),
      _bufE(other._bufE), _bufEL(other._bufEL), _bufEH(other._bufEH),
      _nativeBufE(other._nativeBufE), _nativeBufEL(other._nativeBufEL), _nativeBufEH(other._nativeBufEH) {
      _vecE = (other._vecE) ? new std::vector<Double_t>(*other._vecE) : 0 ;
      _vecEL = (other._vecEL) ? new std::vector<Double_t>(*other._vecEL) : 0 ;
      _vecEH = (other._vecEH) ? new std::vector<Double_t>(*other._vecEH) : 0 ;
    }

    RealFullVector(const RealVector& other, RooAbsReal* real=0) : RealVector(other,real), 
      _bufE(0), _bufEL(0), _bufEH(0),
      _nativeBufE(0), _nativeBufEL(0), _nativeBufEH(0) {
      _vecE = 0 ;
      _vecEL = 0 ; 
      _vecEH = 0 ; 
    }

    RealFullVector& operator=(const RealFullVector& other) {
      if (&other==this) return *this;
      RealVector::operator=(other);
      _bufE = other._bufE;
      _bufEL = other._bufEL;
      _bufEH = other._bufEH;
      _nativeBufE = other._nativeBufE;
      _nativeBufEL = other._nativeBufEL;
      _nativeBufEH = other._nativeBufEH;
      std::vector<Double_t>* src[3] = { other._vecE, other._vecEL, other._vecEH };
      std::vector<Double_t>* dst[3] = { _vecE, _vecEL, _vecEH };
      for (unsigned i = 0; i < 3; ++i) {
        if (src[i]) {
          if (dst[i]) {
            if (dst[i]->size() <= src[i]->capacity() / 2 &&
                src[i]->capacity() > (VECTOR_BUFFER_SIZE / sizeof(Double_t))) {
              std::vector<Double_t> tmp;
              tmp.reserve(std::max(src[i]->size(), VECTOR_BUFFER_SIZE / sizeof(Double_t)));
              tmp.assign(src[i]->begin(), src[i]->end());
              dst[i]->swap(tmp);
            } else {
              *dst[i] = *src[i];
            }
          } else {
            dst[i] = new std::vector<Double_t>(*src[i]);
          }
        } else {
          delete dst[i];
          dst[i] = 0;
        }
      }
      return *this;
    }

    void setErrorBuffer(Double_t* newBuf) { 
      /*       std::cout << "setErrorBuffer(" << _nativeReal->GetName() << ") newBuf = " << newBuf << std::endl ; */
      _bufE = newBuf ; 
      if (!_vecE) _vecE = new std::vector<Double_t> ;
      _vecE->reserve(_vec.capacity()) ;
      if (!_nativeBufE) _nativeBufE = _bufE ;
    }
    void setAsymErrorBuffer(Double_t* newBufL, Double_t* newBufH) { 
      _bufEL = newBufL ; _bufEH = newBufH ; 
      if (!_vecEL) {
        _vecEL = new std::vector<Double_t> ;
        _vecEH = new std::vector<Double_t> ;
        _vecEL->reserve(_vec.capacity()) ;
        _vecEH->reserve(_vec.capacity()) ;
      }
      if (!_nativeBufEL) {
        _nativeBufEL = _bufEL ;
        _nativeBufEH = _bufEH ;
      }
    }

    inline void getNative(Int_t idx) const { 
      RealVector::getNative(idx) ;
      if (_vecE) {
        *_nativeBufE = (*_vecE)[idx] ;
      }
      if (_vecEL) {
        *_nativeBufEL = (*_vecEL)[idx] ;
        *_nativeBufEH = (*_vecEH)[idx] ;
      }
    }
        
    void fill() { 
      RealVector::fill() ;
      if (_vecE) _vecE->push_back(*_bufE) ;
      if (_vecEL) _vecEL->push_back(*_bufEL) ;
      if (_vecEH) _vecEH->push_back(*_bufEH) ;
    } ;

    void write(Int_t i) {
      RealVector::write(i) ;
      if (_vecE) (*_vecE)[i] = *_bufE ;
      if (_vecEL) (*_vecEL)[i] = *_bufEL ;
      if (_vecEH) (*_vecEH)[i] = *_bufEH ;
    }
    
    void reset() { 
      RealVector::reset();
      if (_vecE) {
        std::vector<Double_t> tmp;
        _vecE->swap(tmp);
      }
      if (_vecEL) {
        std::vector<Double_t> tmp;
        _vecEL->swap(tmp);
      }
      if (_vecEH) {
        std::vector<Double_t> tmp;
        _vecEH->swap(tmp);
      }
    }

    inline void get(Int_t idx) const { 
      RealVector::get(idx) ;
      if (_vecE) *_bufE = (*_vecE)[idx];
      if (_vecEL) *_bufEL = (*_vecEL)[idx] ;
      if (_vecEH) *_bufEH = (*_vecEH)[idx] ;
    }

    void resize(Int_t siz) {
      RealVector::resize(siz);
      std::vector<Double_t>* vlist[3] = { _vecE, _vecEL, _vecEH };
      for (unsigned i = 0; i < 3; ++i) {
        if (!vlist[i]) continue;
        if (vlist[i]) {
          if (siz < Int_t(vlist[i]->capacity()) / 2 && vlist[i]->capacity() > (VECTOR_BUFFER_SIZE / sizeof(Double_t))) {
            // if we gain a factor of 2 in memory, we copy and swap
            std::vector<Double_t> tmp;
            tmp.reserve(std::max(siz, Int_t(VECTOR_BUFFER_SIZE / sizeof(Double_t))));
            if (!vlist[i]->empty())
              tmp.assign(vlist[i]->begin(),
                  std::min(_vec.end(), _vec.begin() + siz));
            if (Int_t(tmp.size()) != siz)
              tmp.resize(siz);
            vlist[i]->swap(tmp);
          } else {
            vlist[i]->resize(siz);
          }
        }
      }
    }

    void reserve(Int_t siz) {
      RealVector::reserve(siz);
      if (_vecE) _vecE->reserve(siz);
      if (_vecEL) _vecEL->reserve(siz);
      if (_vecEH) _vecEH->reserve(siz);
    }

  private:
    friend class RooVectorDataStore ;
    Double_t *_bufE ; //!
    Double_t *_bufEL ; //!
    Double_t *_bufEH ; //!
    Double_t *_nativeBufE ; //!
    Double_t *_nativeBufEL ; //! 
    Double_t *_nativeBufEH ; //!
    std::vector<double> *_vecE, *_vecEL, *_vecEH ;
    ClassDef(RealFullVector,1) // STL-vector-based Data Storage class
  } ;
  

  class CatVector {
  public:
    CatVector(UInt_t initialCapacity=(VECTOR_BUFFER_SIZE / sizeof(RooCatType))) : 
      _cat(0), _buf(0), _nativeBuf(0), _vec0(0)
    {
      _vec.reserve(initialCapacity);
    }

    CatVector(RooAbsCategory* cat, UInt_t initialCapacity=(VECTOR_BUFFER_SIZE / sizeof(RooCatType))) : 
      _cat(cat), _buf(0), _nativeBuf(0), _vec0(0)
    {
      _vec.reserve(initialCapacity);
    }

    virtual ~CatVector() {
    }

    CatVector(const CatVector& other, RooAbsCategory* cat=0) : 
      _cat(cat?cat:other._cat), _buf(other._buf), _nativeBuf(other._nativeBuf), _vec(other._vec) 
      {
	_vec0 = _vec.size()>0 ? &_vec.front() : 0 ;
      }

    CatVector& operator=(const CatVector& other) {
      if (&other==this) return *this;
      _cat = other._cat;
      _buf = other._buf;
      _nativeBuf = other._nativeBuf;
      if (other._vec.size() <= _vec.capacity() / 2 && _vec.capacity() > (VECTOR_BUFFER_SIZE / sizeof(RooCatType))) {
	std::vector<RooCatType> tmp;
	tmp.reserve(std::max(other._vec.size(), VECTOR_BUFFER_SIZE / sizeof(RooCatType)));
	tmp.assign(other._vec.begin(), other._vec.end());
	_vec.swap(tmp);
      } else {
	_vec = other._vec;
      }
      _vec0 = _vec.size()>0 ? &_vec.front() : 0;
      return *this;
    }

    void setBuffer(RooCatType* newBuf) { 
      _buf = newBuf ; 
      if (_nativeBuf==0) _nativeBuf=newBuf ;
    }

    void setNativeBuffer(RooCatType* newBuf=0) {       
      _nativeBuf = newBuf ? newBuf : _buf ; 
    }
    
    void fill() { 
      _vec.push_back(*_buf) ; 
      _vec0 = &_vec.front() ;
    } ;
    void write(Int_t i) { 
      _vec[i]=*_buf ; 
    } ;
    void reset() { 
      // make sure the vector releases the underlying memory
      std::vector<RooCatType> tmp;
      _vec.swap(tmp);
      _vec0 = 0;
    }
    inline void get(Int_t idx) const { 
      _buf->assignFast(*(_vec0+idx)) ;
    }
    inline void getNative(Int_t idx) const { 
      _nativeBuf->assignFast(*(_vec0+idx)) ;
    }
    Int_t size() const { return _vec.size() ; }

    void resize(Int_t siz) {
      if (siz < Int_t(_vec.capacity()) / 2 && _vec.capacity() > (VECTOR_BUFFER_SIZE / sizeof(RooCatType))) {
	// do an expensive copy, if we save at least a factor 2 in size
	std::vector<RooCatType> tmp;
	tmp.reserve(std::max(siz, Int_t(VECTOR_BUFFER_SIZE / sizeof(RooCatType))));
	if (!_vec.empty())
	    tmp.assign(_vec.begin(), std::min(_vec.end(), _vec.begin() + siz));
	if (Int_t(tmp.size()) != siz) 
	    tmp.resize(siz);
	_vec.swap(tmp);
      } else {
	_vec.resize(siz);
      }
      _vec0 = _vec.size() > 0 ? &_vec.front() : 0;
    }

    void reserve(Int_t siz) {
      _vec.reserve(siz);
      _vec0 = _vec.size() > 0 ? &_vec.front() : 0;
    }

    void setBufArg(RooAbsCategory* arg) { _cat = arg; }
    const RooAbsCategory* bufArg() const { return _cat; }

  private:
    friend class RooVectorDataStore ;
    RooAbsCategory* _cat ;
    RooCatType* _buf ;  //!
    RooCatType* _nativeBuf ;  //!
    std::vector<RooCatType> _vec ;
    RooCatType* _vec0 ; //!
    ClassDef(CatVector,1) // STL-vector-based Data Storage class
  } ;
  

 protected:

  friend class RooAbsReal ;
  friend class RooAbsCategory ;
  friend class RooRealVar ;
  std::vector<RealVector*>& realStoreList() { return _realStoreList ; }
  std::vector<RealFullVector*>& realfStoreList() { return _realfStoreList ; }
  std::vector<CatVector*>& catStoreList() { return _catStoreList ; }

  CatVector* addCategory(RooAbsCategory* cat) {

    // First try a match by name
    for (auto catVec : _catStoreList) {
      if (std::string(catVec->bufArg()->GetName())==cat->GetName()) {
        return catVec;
      }
    }

    // If nothing found this will make an entry
    _catStoreList.push_back(new CatVector(cat)) ;
    _nCat++ ;

    // Update cached ptr to first element as push_back may have reallocated
    _firstCat = &_catStoreList.front() ;

    return _catStoreList.back() ;
  }

  RealVector* addReal(RooAbsReal* real) {
    
    // First try a match by name
    for (auto realVec : _realStoreList) {
      if (realVec->bufArg()->namePtr()==real->namePtr()) {
        return realVec;
      }
    }    

    // Then check if an entry already exists for a full real    
    for (auto fullVec : _realfStoreList) {
      if (fullVec->bufArg()->namePtr()==real->namePtr()) {
        // Return full vector as RealVector base class here
        return fullVec;
      }
    }    

    // If nothing found this will make an entry
    _realStoreList.push_back(new RealVector(real)) ;

    return _realStoreList.back() ;
  }

  Bool_t isFullReal(RooAbsReal* real) {
    
    // First try a match by name
    for (auto fullVec : _realfStoreList) {
      if (std::string(fullVec->bufArg()->GetName())==real->GetName()) {
        return kTRUE ;
      }
    }        
    return kFALSE ;
  }

  Bool_t hasError(RooAbsReal* real) {
    
    // First try a match by name
    for (auto fullVec : _realfStoreList) {
      if (std::string(fullVec->bufArg()->GetName())==real->GetName()) {
        return fullVec->_vecE ? kTRUE : kFALSE ;
      }
    }        
    return kFALSE ;
  }

  Bool_t hasAsymError(RooAbsReal* real) {
    
    // First try a match by name
    for (auto fullVec : _realfStoreList) {
      if (std::string(fullVec->bufArg()->GetName())==real->GetName()) {
        return fullVec->_vecEL ? kTRUE : kFALSE ;
      }
    }        
    return kFALSE ;
  }

  RealFullVector* addRealFull(RooAbsReal* real) {
    
    // First try a match by name
    for (auto fullVec : _realfStoreList) {
      if (std::string(fullVec->bufArg()->GetName())==real->GetName()) {
	    return fullVec;
      }
    }    

    // Then check if an entry already exists for a bare real    
    for (auto realVec : _realStoreList) {
      if (std::string(realVec->bufArg()->GetName())==real->GetName()) {

        // Convert element to full and add to full list
        _realfStoreList.push_back(new RealFullVector(*realVec,real)) ;
        _nRealF++ ;
        _firstRealF = &_realfStoreList.front() ;

        // Delete bare element
        _realStoreList.erase(std::find(_realStoreList.begin(), _realStoreList.end(), realVec));
        delete realVec;

        return _realfStoreList.back() ;
      }
    }    

    // If nothing found this will make an entry
    _realfStoreList.push_back(new RealFullVector(real)) ;
    _nRealF++ ;

    // Update cached ptr to first element as push_back may have reallocated
    _firstRealF = &_realfStoreList.front() ;


    return _realfStoreList.back() ;
  }

  virtual Bool_t hasFilledCache() const override { return _cache ? kTRUE : kFALSE ; }

  void forceCacheUpdate() override;

 private:
  RooArgSet _varsww ;
  RooRealVar* _wgtVar ;     // Pointer to weight variable (if set)

  std::vector<RealVector*> _realStoreList ;
  std::vector<RealFullVector*> _realfStoreList ;
  std::vector<CatVector*> _catStoreList ;
  std::vector<double> _weights;

  void setAllBuffersNative() ;

  Int_t _nRealF ;
  Int_t _nCat ;
  Int_t _nEntries ;
  RealFullVector** _firstRealF ; //! do not persist
  CatVector** _firstCat ; //! do not persist
  Double_t _sumWeight ; 
  Double_t _sumWeightCarry;

  const Double_t* _extWgtArray ;         //! External weight array
  const Double_t* _extWgtErrLoArray ;    //! External weight array - low error
  const Double_t* _extWgtErrHiArray ;    //! External weight array - high error
  const Double_t* _extSumW2Array ;       //! External sum of weights array

  mutable Double_t  _curWgt ;      // Weight of current event
  mutable Double_t  _curWgtErrLo ; // Weight of current event
  mutable Double_t  _curWgtErrHi ; // Weight of current event
  mutable Double_t  _curWgtErr ;   // Weight of current event

  RooVectorDataStore* _cache ; //! Optimization cache
  RooAbsArg* _cacheOwner ; //! Cache owner

  Bool_t _forcedUpdate ; //! Request for forced cache update 

  ClassDefOverride(RooVectorDataStore,3) // STL-vector-based Data Storage class
};


#endif
