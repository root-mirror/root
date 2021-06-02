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

#include "RooAbsDataStore.h"
#include "RooAbsCategory.h"
#include "RooAbsReal.h"
#include "RooChangeTracker.h"
#include "RooRealVar.h"

#include "ROOT/RStringView.hxx"

#include <list>
#include <vector>
#include <algorithm>

class RooAbsArg ;
class RooArgList ;
class TTree ;
class RooFormulaVar ;
class RooArgSet ;

#define VECTOR_BUFFER_SIZE 1024

class RooVectorDataStore : public RooAbsDataStore {
public:

  RooVectorDataStore() ; 

  // Empty ctor
  RooVectorDataStore(std::string_view name, std::string_view title, const RooArgSet& vars, const char* wgtVarName=0) ;

  WRITE_TSTRING_COMPATIBLE_CONSTRUCTOR(RooVectorDataStore)

  virtual RooAbsDataStore* clone(const char* newname=0) const override { return new RooVectorDataStore(*this,newname) ; }
  virtual RooAbsDataStore* clone(const RooArgSet& vars, const char* newname=0) const override { return new RooVectorDataStore(*this,vars,newname) ; }

  RooVectorDataStore(const RooVectorDataStore& other, const char* newname=0) ;
  RooVectorDataStore(const RooTreeDataStore& other, const RooArgSet& vars, const char* newname=0) ;
  RooVectorDataStore(const RooVectorDataStore& other, const RooArgSet& vars, const char* newname=0) ;


  RooVectorDataStore(std::string_view name, std::string_view title, RooAbsDataStore& tds, 
		     const RooArgSet& vars, const RooFormulaVar* cutVar, const char* cutRange,
		     std::size_t nStart, std::size_t nStop, Bool_t /*copyCache*/, const char* wgtVarName=0) ;

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

  using RooAbsDataStore::weight ;
  /// Return the weight of the last-retrieved data point.
  Double_t weight() const override
  {
    if (_extWgtArray)
      return _extWgtArray[_currentWeightIndex];
    if (_wgtVar)
      return _wgtVar->getVal();

    return 1.0;
  }
  virtual Double_t weightError(RooAbsData::ErrorType etype=RooAbsData::Poisson) const override;
  virtual void weightError(Double_t& lo, Double_t& hi, RooAbsData::ErrorType etype=RooAbsData::Poisson) const override;
  virtual Bool_t isWeighted() const override { return _wgtVar || _extWgtArray; }

  RooBatchCompute::RunContext getBatches(std::size_t first, std::size_t len) const override;
  virtual RooSpan<const double> getWeightBatch(std::size_t first, std::size_t len) const override;

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
  virtual Int_t numEntries() const override { return static_cast<int>(size()); }
  virtual Double_t sumEntries() const override { return _sumWeight ; }
  /// Get size of stored dataset.
  std::size_t size() const {
    if (!_realStoreList.empty()) {
      return _realStoreList.front()->size();
    } else if (!_realfStoreList.empty()) {
      return _realfStoreList.front()->size();
    } else if (!_catStoreList.empty()) {
      return _catStoreList.front()->size();
    }

    return 0;
  }
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

  void loadValues(const RooAbsDataStore *tds, const RooFormulaVar* select=0, const char* rangeName=0, std::size_t nStart=0, std::size_t nStop = std::numeric_limits<std::size_t>::max()) override;
  
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
      assert(static_cast<std::size_t>(i) < _vec.size());
      _vec[i] = *_buf ;
    }
    
    void reset() { 
      _vec.clear();
    }

    inline void load(std::size_t idx) const {
      assert(idx < _vec.size());
      *_buf = *(_vec.begin() + idx) ;
    }

    RooSpan<const double> getRange(std::size_t first, std::size_t last) const {
      auto beg = std::min(_vec.cbegin() + first, _vec.cend());
      auto end = std::min(_vec.cbegin() + last,  _vec.cend());

      return RooSpan<const double>(beg, end);
    }

    inline void loadToNative(std::size_t idx) const {
      *_nativeBuf = *(_vec.begin() + idx) ;
    }

    std::size_t size() const { return _vec.size() ; }

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
    RooAbsReal* _nativeReal ; // Instance which our data belongs to. This is the variable in the dataset.
    RooAbsReal* _real ; // Instance where we should write data into when load() is called.
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

    inline void loadToNative(Int_t idx) const { 
      RealVector::loadToNative(idx) ;
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
      RealVector::load(idx) ;
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
    CatVector(UInt_t initialCapacity = VECTOR_BUFFER_SIZE) :
      _cat(nullptr), _buf(nullptr), _nativeBuf(nullptr)
    {
      _vec.reserve(initialCapacity);
    }

    CatVector(RooAbsCategory* cat, UInt_t initialCapacity = VECTOR_BUFFER_SIZE) :
      _cat(cat), _buf(nullptr), _nativeBuf(nullptr)
    {
      _vec.reserve(initialCapacity);
    }

    virtual ~CatVector() {
    }

    CatVector(const CatVector& other, RooAbsCategory* cat = nullptr) :
      _cat(cat?cat:other._cat), _buf(other._buf), _nativeBuf(other._nativeBuf), _vec(other._vec) 
    {

    }

    CatVector& operator=(const CatVector& other) {
      if (&other==this) return *this;
      _cat = other._cat;
      _buf = other._buf;
      _nativeBuf = other._nativeBuf;
      if (other._vec.size() <= _vec.capacity() / 2 && _vec.capacity() > VECTOR_BUFFER_SIZE) {
        std::vector<RooAbsCategory::value_type> tmp;
        tmp.reserve(std::max(other._vec.size(), std::size_t(VECTOR_BUFFER_SIZE)));
        tmp.assign(other._vec.begin(), other._vec.end());
        _vec.swap(tmp);
      } else {
        _vec = other._vec;
      }

      return *this;
    }

    void setBuffer(RooAbsCategory::value_type* newBuf) {
      _buf = newBuf ; 
      if (!_nativeBuf) _nativeBuf = newBuf;
    }

    void setNativeBuffer(RooAbsCategory::value_type* newBuf = nullptr) {
      _nativeBuf = newBuf ? newBuf : _buf;
    }
    
    void fill() { 
      _vec.push_back(*_buf) ; 
    }

    void write(std::size_t i) {
      _vec[i] = *_buf;
    }

    void reset() { 
      // make sure the vector releases the underlying memory
      std::vector<RooAbsCategory::value_type> tmp;
      _vec.swap(tmp);
    }

    inline void load(std::size_t idx) const {
      *_buf = _vec[idx];
    }

    inline void loadToNative(std::size_t idx) const {
      *_nativeBuf = _vec[idx];
    }

    std::size_t size() const { return _vec.size() ; }

    void resize(Int_t siz) {
      if (siz < Int_t(_vec.capacity()) / 2 && _vec.capacity() > VECTOR_BUFFER_SIZE) {
        // do an expensive copy, if we save at least a factor 2 in size
        std::vector<RooAbsCategory::value_type> tmp;
        tmp.reserve(std::max(siz, VECTOR_BUFFER_SIZE));
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

    void setBufArg(RooAbsCategory* arg) { _cat = arg; }
    const RooAbsCategory* bufArg() const { return _cat; }

  private:
    friend class RooVectorDataStore ;
    RooAbsCategory* _cat;
    RooAbsCategory::value_type* _buf;  //!
    RooAbsCategory::value_type* _nativeBuf;  //!
    std::vector<RooAbsCategory::value_type> _vec;
    ClassDef(CatVector,2) // STL-vector-based Data Storage class
  } ;
  

 protected:

  friend class RooAbsReal ;
  friend class RooAbsCategory ;
  friend class RooRealVar ;
  std::vector<RealVector*>& realStoreList() { return _realStoreList ; }
  std::vector<RealFullVector*>& realfStoreList() { return _realfStoreList ; }
  std::vector<CatVector*>& catStoreList() { return _catStoreList ; }

  CatVector* addCategory(RooAbsCategory* cat);

  RealVector* addReal(RooAbsReal* real);

  Bool_t isFullReal(RooAbsReal* real);

  Bool_t hasError(RooAbsReal* real);

  Bool_t hasAsymError(RooAbsReal* real);

  RealFullVector* addRealFull(RooAbsReal* real);

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

  Double_t _sumWeight ; 
  Double_t _sumWeightCarry;

  const Double_t* _extWgtArray ;         //! External weight array
  const Double_t* _extWgtErrLoArray ;    //! External weight array - low error
  const Double_t* _extWgtErrHiArray ;    //! External weight array - high error
  const Double_t* _extSumW2Array ;       //! External sum of weights array

  mutable std::size_t _currentWeightIndex{0}; //

  RooVectorDataStore* _cache ; //! Optimization cache
  RooAbsArg* _cacheOwner ; //! Cache owner

  Bool_t _forcedUpdate ; //! Request for forced cache update 

  ClassDefOverride(RooVectorDataStore, 5) // STL-vector-based Data Storage class
};


#endif
