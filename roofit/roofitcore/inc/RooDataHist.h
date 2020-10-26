/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooDataHist.h,v 1.37 2007/05/11 09:11:30 verkerke Exp $
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
#ifndef ROO_DATA_HIST
#define ROO_DATA_HIST

#include "RooAbsData.h"
#include "RooDirItem.h"
#include "RooArgSet.h"
#include "RooNameSet.h"
#include "RooCacheManager.h"

#include <map>
#include <vector>
#include <string>

class TObject ;
class RooAbsArg;
class RooAbsReal ;
class RooAbsCategory ;
class Roo1DTable ;
class RooPlot;
class RooArgSet ;
class RooLinkedList ;
class RooAbsLValue ;

class RooDataHist : public RooAbsData, public RooDirItem {
public:

  // Constructors, factory methods etc.
  RooDataHist() ; 
  RooDataHist(const char *name, const char *title, const RooArgSet& vars, const char* binningName=0) ;
  RooDataHist(const char *name, const char *title, const RooArgSet& vars, const RooAbsData& data, Double_t initWgt=1.0) ;
  RooDataHist(const char *name, const char *title, const RooArgList& vars, const TH1* hist, Double_t initWgt=1.0) ;
  RooDataHist(const char *name, const char *title, const RooArgList& vars, RooCategory& indexCat, std::map<std::string,TH1*> histMap, Double_t initWgt=1.0) ;
  RooDataHist(const char *name, const char *title, const RooArgList& vars, RooCategory& indexCat, std::map<std::string,RooDataHist*> dhistMap, Double_t wgt=1.0) ;
  //RooDataHist(const char *name, const char *title, const RooArgList& vars, Double_t initWgt=1.0) ;
  RooDataHist(const char *name, const char *title, const RooArgList& vars, const RooCmdArg& arg1, const RooCmdArg& arg2=RooCmdArg(), const RooCmdArg& arg3=RooCmdArg(),
        const RooCmdArg& arg4=RooCmdArg(),const RooCmdArg& arg5=RooCmdArg(),const RooCmdArg& arg6=RooCmdArg(),const RooCmdArg& arg7=RooCmdArg(),const RooCmdArg& arg8=RooCmdArg()) ;
  RooDataHist& operator=(const RooDataHist&) = delete;

  RooDataHist(const RooDataHist& other, const char* newname = 0) ;
  TObject* Clone(const char* newname="") const override {
    return new RooDataHist(*this, newname && newname[0] != '\0' ? newname : GetName());
  }
  ~RooDataHist() override ;

  /// Return empty clone of this RooDataHist.
  RooAbsData* emptyClone(const char* newName=0, const char* newTitle=0, const RooArgSet*vars=0, const char* /*wgtVarName*/=0) const override {
    return new RooDataHist(newName?newName:GetName(),newTitle?newTitle:GetTitle(),vars?*vars:*get()) ; 
  }

  /// Add `wgt` to the bin content enclosed by the coordinates passed in `row`.
  virtual void add(const RooArgSet& row, Double_t wgt=1.0) { add(row,wgt,-1.); }
  void add(const RooArgSet& row, Double_t weight, Double_t sumw2) override ;
  void set(std::size_t binNumber, double weight, double wgtErr);
  void set(const RooArgSet& row, Double_t weight, Double_t wgtErr=-1.) ;
  void set(const RooArgSet& row, Double_t weight, Double_t wgtErrLo, Double_t wgtErrHi) ;

  void add(const RooAbsData& dset, const RooFormulaVar* cutVar=0, Double_t weight=1.0 ) ;
  void add(const RooAbsData& dset, const char* cut, Double_t weight=1.0 ) ;

  /// Get bin centre of current bin.
  const RooArgSet* get() const override { return &_vars; }
  const RooArgSet* get(Int_t binNumber) const override;
  virtual const RooArgSet* get(const RooArgSet& coord) const;
  Int_t numEntries() const override;
  Double_t sumEntries() const override;
  Double_t sumEntries(const char* cutSpec, const char* cutRange=0) const override;

  /// Always returns true as all histograms have in principle events weight != 1.
  Bool_t isWeighted() const override { return true; }
  Bool_t isNonPoissonWeighted() const override ;

  RooSpan<const double> getWeightBatch(std::size_t first, std::size_t len) const override;
  void getBatches(BatchHelpers::RunContext& evalData, std::size_t begin, std::size_t len) const override;

  Double_t sum(Bool_t correctForBinSize, Bool_t inverseCorr=kFALSE) const ;
  Double_t sum(const RooArgSet& sumSet, const RooArgSet& sliceSet, Bool_t correctForBinSize, Bool_t inverseCorr=kFALSE) ;
  Double_t sum(const RooArgSet& sumSet, const RooArgSet& sliceSet, Bool_t correctForBinSize, Bool_t inverseCorr, const std::map<const RooAbsArg*, std::pair<Double_t, Double_t> >& ranges);

  /// Return weight of i-th bin. \see getIndex()
  double weight(std::size_t i) const { return _wgtVec[i]; }
  Double_t weight(const RooArgSet& bin, Int_t intOrder=1, Bool_t correctForBinSize=kFALSE, Bool_t cdfBoundaries=kFALSE, Bool_t oneSafe=kFALSE);
  /// Return squared weight sum of i-th bin. \see getIndex()
  double weightSquared(std::size_t i) const { return get_sumw2(i); }
  /// Return bin volume of i-th bin. \see getIndex()
  double binVolume(std::size_t i) const { return _binvVec[i]; }
  double binVolume(const RooArgSet& bin) const; 
  /// Return true if bin `i` is considered valid within the current range definitions of all observables. \see getIndex()
  bool valid(std::size_t i) const { return i <= _wgtVec.size() && (_maskedWeights.empty() || _maskedWeights[i] != 0.);}

  TIterator* sliceIterator(RooAbsArg& sliceArg, const RooArgSet& otherArgs) ;

  void weightError(Double_t& lo, Double_t& hi, ErrorType etype=Poisson) const override;
  /// Return the error of the weight of the last-retrieved entry. See also weightError(Double_t&,Double_t&,ErrorType) const.
  Double_t weightError(ErrorType etype=Poisson) const override {
    // Return symmetric error on current bin calculated either from Poisson statistics or from SumOfWeights
    Double_t lo,hi ;
    weightError(lo,hi,etype) ;
    return (lo+hi)/2 ;
  }

  using RooAbsData::plotOn ;
  RooPlot *plotOn(RooPlot *frame, PlotOpt o) const override;

  void reset() override;

  virtual void printMultiline(std::ostream& os, Int_t content, Bool_t verbose=kFALSE, TString indent="") const override;
  virtual void printArgs(std::ostream& os) const override;
  virtual void printValue(std::ostream& os) const override;

  void SetName(const char *name) override;
  void SetNameTitle(const char *name, const char* title) override;

  Int_t getIndex(const RooArgSet& coord, Bool_t fast = false) const;

  void removeSelfFromDir() { removeFromDir(this) ; }

  // A shortcut function only for RooAbsOptTestStatistic.
  void cacheValidEntries();


  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// @name Deprecated functions
  /// These functions rely on the fact that an event has been loaded before they are called. It is advised
  /// to switch to their counterparts that take bin numbers as arguments. In this way, code like,
  /// ```
  ///   const RooArgSet* coordinates = dataHist.get(i); // Need this to achieve side effect on next call of weight() - bad.
  ///   const double weight = dataHist.weight();
  ///   processEvent(coordinates, weight);
  /// ```
  /// becomes
  /// ```
  ///   processEvent(dataHist.get(i), dataHist.weight(i));
  /// ```
  /// The index of a set of coordinates can be computed using getIndex().
  /// @{

  /// Return weight of last bin that was requested with get().
  /// \deprecated Use the safer weight(std::size_t) const.
  Double_t weight() const override { return get_curWeight(); }
  /// Return squared weight of last bin that was requested with get().
  /// \deprecated Use the safer weightSquared(std::size_t) const.
  Double_t weightSquared() const override { return get_curSumW2(); }
  /// Return volume of current bin. \deprecated Use binVolume(std::size_t) const.
  Double_t binVolume() const { return _binvVec[_curIndex]; }
  /// Write `weight` into current bin. \deprecated Use set(std::size_t,double,double)
  void set(Double_t weight, Double_t wgtErr=-1);

  /// Return true if currently loaded coordinate is considered valid within
  /// the current range definitions of all observables.
  /// \deprecated Use the safer valid(std::size_t) const.
  bool valid() const override { return _curIndex <= _wgtVec.size() && (_maskedWeights.empty() || _maskedWeights[_curIndex] != 0.);}

  void dump2();

  ///@}
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////

protected:

  friend class RooAbsCachedPdf ;
  friend class RooAbsCachedReal ;
  friend class RooDataHistSliceIter ;

  std::size_t calcTreeIndex(const RooArgSet& coords, bool fast) const;

  void setAllWeights(Double_t value) ;
 
  void initialize(const char* binningName=0,Bool_t fillTree=kTRUE) ;
  RooDataHist(const char* name, const char* title, RooDataHist* h, const RooArgSet& varSubset, 
        const RooFormulaVar* cutVar, const char* cutRange, Int_t nStart, Int_t nStop, Bool_t copyCache) ;
  RooAbsData* reduceEng(const RooArgSet& varSubset, const RooFormulaVar* cutVar, const char* cutRange=0, 
                  std::size_t nStart=0, std::size_t nStop=std::numeric_limits<std::size_t>::max(), Bool_t copyCache=kTRUE) override;
  Double_t interpolateDim(RooRealVar& dim, const RooAbsBinning* binning, Double_t xval, Int_t intOrder, Bool_t correctForBinSize, Bool_t cdfBoundaries) ;
  void calculatePartialBinVolume(const RooArgSet& dimSet) const ;
  void checkBinBounds() const;

  void adjustBinning(const RooArgList& vars, const TH1& href, Int_t* offset=0) ;
  void importTH1(const RooArgList& vars, const TH1& histo, Double_t initWgt, Bool_t doDensityCorrection) ;
  void importTH1Set(const RooArgList& vars, RooCategory& indexCat, std::map<std::string,TH1*> hmap, Double_t initWgt, Bool_t doDensityCorrection) ;
  void importDHistSet(const RooArgList& vars, RooCategory& indexCat, std::map<std::string,RooDataHist*> dmap, Double_t initWgt) ;

  RooAbsData* cacheClone(const RooAbsArg* newCacheOwner, const RooArgSet* newCacheVars, const char* newName=0) override ;

  Double_t get_wgt(std::size_t idx)   const { return _wgtVec[idx]; }
  Double_t get_errLo(std::size_t idx) const { return _errLoVec.empty() ? -1 : _errLoVec[idx]; }
  Double_t get_errHi(std::size_t idx) const { return _errHiVec.empty() ? -1 : _errHiVec[idx]; }
  /// Return sumw^2 of bin `i`. If this is not being tracked, assume that all fill operations
  /// had a weight of 1.
  Double_t get_sumw2(std::size_t idx) const { return _sumw2Vec.empty() ? _wgtVec[idx] : _sumw2Vec[idx];}

  Double_t get_curWeight()   const { return get_wgt(_curIndex); }
  Double_t get_curWgtErrLo() const { return get_errLo(_curIndex); }
  Double_t get_curWgtErrHi() const { return get_errHi(_curIndex); }
  Double_t get_curSumW2()    const { return get_sumw2(_curIndex); }

  Int_t get_curIndex() const { return _curIndex; }

  Int_t _arrSize; /// For backward-compatible I/O.
  std::vector<Int_t> _idxMult ; // Multiplier jump table for index calculation

  double* _wgt{nullptr};   //[_arrSize] For backward-compatible I/O. Used to be [_arrSize] Weight array
  double* _errLo{nullptr}; //[_arrSize] For backward-compatible I/O. Used to be [_arrSize] Low-side error on weight array
  double* _errHi{nullptr}; //[_arrSize] For backward-compatible I/O. Used to be [_arrSize] High-side error on weight array
  double* _sumw2{nullptr}; //[_arrSize] For backward-compatible I/O. Used to be [_arrSize] Sum of weights^2
  double* _binv{nullptr};  //[_arrSize] For backward-compatible I/O. Used to be [_arrSize] Bin volume array

  std::vector<double> _wgtVec;   /// Weight of each bin
  std::vector<double> _errLoVec; /// Low-side error of weight
  std::vector<double> _errHiVec; /// High-side error of weight
  std::vector<double> _sumw2Vec; /// Sum of weights^2
  std::vector<double> _binvVec;  /// Bin volumes

  RooArgSet  _realVars ; // Real dimensions of the dataset 
  mutable std::vector<double> _maskedWeights; //! Copy of _wgtVec, but masked events have a weight of zero.
 
  mutable std::size_t _curIndex{std::numeric_limits<std::size_t>::max()}; // Current index

  mutable std::vector<Double_t>* _pbinv ; //! Partial bin volume array
  mutable RooCacheManager<std::vector<Double_t> > _pbinvCacheMgr ; //! Cache manager for arrays of partial bin volumes
  std::vector<RooAbsLValue*> _lvvars ; //! List of observables casted as RooAbsLValue
  std::vector<const RooAbsBinning*> _lvbins ; //! List of used binnings associated with lvalues
  mutable std::vector<std::vector<Double_t> > _binbounds; //! list of bin bounds per dimension

  enum CacheSumState_t{kInvalid, kValid, kValidCorrectForBinSize, kValidInvBinCorr};
  mutable CacheSumState_t _cache_sum_valid{kInvalid}; //! Is cache sum valid
  mutable Double_t _cache_sum{0.}; //! Cache for sum of entries ;


private:
  void _adjustBinning(RooRealVar &theirVar, const TAxis &axis, RooRealVar *ourVar, Int_t *offset);
  void registerWeightArraysToDataStore() const;

  ClassDefOverride(RooDataHist, 5) // Binned data set
};

#endif

