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
  virtual TObject* Clone(const char* newname="") const {
    return new RooDataHist(*this, newname && newname[0] != '\0' ? newname : GetName());
  }
  virtual ~RooDataHist() ;

  virtual RooAbsData* emptyClone(const char* newName=0, const char* newTitle=0, const RooArgSet*vars=0, const char* /*wgtVarName*/=0) const {
    // Return empty clone of this RooDataHist
    return new RooDataHist(newName?newName:GetName(),newTitle?newTitle:GetTitle(),vars?*vars:*get()) ; 
  }

  // Add one ore more rows of data
  virtual void add(const RooArgSet& row, Double_t wgt=1.0) { 
    // Increment weight of bin enclosing coordinate stored in row by wgt
    add(row,wgt,-1.) ; 
  }
  virtual void add(const RooArgSet& row, Double_t weight, Double_t sumw2) ;
  void set(Double_t weight, Double_t wgtErr=-1) ;
  void set(const RooArgSet& row, Double_t weight, Double_t wgtErr=-1) ;
  void set(const RooArgSet& row, Double_t weight, Double_t wgtErrLo, Double_t wgtErrHi) ;

  void add(const RooAbsData& dset, const RooFormulaVar* cutVar=0, Double_t weight=1.0 ) ;
  void add(const RooAbsData& dset, const char* cut, Double_t weight=1.0 ) ;

  virtual const RooArgSet* get() const { 
    // Return set with coordinates of center of current bin
    return &_vars ; 
  } 
  virtual const RooArgSet* get(Int_t masterIdx) const ;
  virtual const RooArgSet* get(const RooArgSet& coord) const ;
  virtual Int_t numEntries() const ; 
  virtual Double_t sumEntries() const  ;
  virtual Double_t sumEntries(const char* cutSpec, const char* cutRange=0) const ;
  virtual Bool_t isWeighted() const { 
    // Return true as all histograms have in principle events weight != 1
    return kTRUE ;     
  }
  virtual Bool_t isNonPoissonWeighted() const ;

  virtual RooSpan<const double> getWeightBatch(std::size_t first, std::size_t len) const;
  void getBatches(BatchHelpers::RunContext& evalData, std::size_t begin, std::size_t len) const;

  Double_t sum(Bool_t correctForBinSize, Bool_t inverseCorr=kFALSE) const ;
  Double_t sum(const RooArgSet& sumSet, const RooArgSet& sliceSet, Bool_t correctForBinSize, Bool_t inverseCorr=kFALSE) ;
  Double_t sum(const RooArgSet& sumSet, const RooArgSet& sliceSet, Bool_t correctForBinSize, Bool_t inverseCorr, const std::map<const RooAbsArg*, std::pair<Double_t, Double_t> >& ranges);

  virtual Double_t weight() const { 
    // Return weight of current bin
    return get_curWeight();
  }
  Double_t weightSquared() const ;
  Double_t weight(const RooArgSet& bin, Int_t intOrder=1, Bool_t correctForBinSize=kFALSE, Bool_t cdfBoundaries=kFALSE, Bool_t oneSafe=kFALSE) ;   
  Double_t binVolume() const { return _curVolume ; }
  Double_t binVolume(const RooArgSet& bin) ; 
  virtual Bool_t valid() const ;

  TIterator* sliceIterator(RooAbsArg& sliceArg, const RooArgSet& otherArgs) ;
  
  virtual void weightError(Double_t& lo, Double_t& hi, ErrorType etype=Poisson) const ;
  virtual Double_t weightError(ErrorType etype=Poisson) const { 
    // Return symmetric error on current bin calculated either from Poisson statistics or from SumOfWeights
    Double_t lo,hi ;
    weightError(lo,hi,etype) ;
    return (lo+hi)/2 ;
  }

  using RooAbsData::plotOn ;
  virtual RooPlot *plotOn(RooPlot *frame, PlotOpt o) const;
  
  virtual void reset() ;
  void dump2() ;

  virtual void printMultiline(std::ostream& os, Int_t content, Bool_t verbose=kFALSE, TString indent="") const ;
  virtual void printArgs(std::ostream& os) const ;
  virtual void printValue(std::ostream& os) const ;

  void SetName(const char *name) ;
  void SetNameTitle(const char *name, const char* title) ;

  Int_t getIndex(const RooArgSet& coord, Bool_t fast=kFALSE) ;

  void removeSelfFromDir() { removeFromDir(this) ; }

protected:

  friend class RooAbsCachedPdf ;
  friend class RooAbsCachedReal ;
  friend class RooDataHistSliceIter ;
  friend class RooAbsOptTestStatistic ;

  Int_t calcTreeIndex() const ;
  void cacheValidEntries() ;

  void setAllWeights(Double_t value) ;
 
  void initialize(const char* binningName=0,Bool_t fillTree=kTRUE) ;
  RooDataHist(const char* name, const char* title, RooDataHist* h, const RooArgSet& varSubset, 
	      const RooFormulaVar* cutVar, const char* cutRange, Int_t nStart, Int_t nStop, Bool_t copyCache) ;
  RooAbsData* reduceEng(const RooArgSet& varSubset, const RooFormulaVar* cutVar, const char* cutRange=0, 
	                std::size_t nStart=0, std::size_t nStop=std::numeric_limits<std::size_t>::max(), Bool_t copyCache=kTRUE) ;
  Double_t interpolateDim(RooRealVar& dim, const RooAbsBinning* binning, Double_t xval, Int_t intOrder, Bool_t correctForBinSize, Bool_t cdfBoundaries) ;
  void calculatePartialBinVolume(const RooArgSet& dimSet) const ;
  void checkBinBounds() const;

  void adjustBinning(const RooArgList& vars, const TH1& href, Int_t* offset=0) ;
  void importTH1(const RooArgList& vars, const TH1& histo, Double_t initWgt, Bool_t doDensityCorrection) ;
  void importTH1Set(const RooArgList& vars, RooCategory& indexCat, std::map<std::string,TH1*> hmap, Double_t initWgt, Bool_t doDensityCorrection) ;
  void importDHistSet(const RooArgList& vars, RooCategory& indexCat, std::map<std::string,RooDataHist*> dmap, Double_t initWgt) ;

  virtual RooAbsData* cacheClone(const RooAbsArg* newCacheOwner, const RooArgSet* newCacheVars, const char* newName=0) ;

  Double_t get_wgt(const Int_t &idx) const { return _wgt[idx]; }
  Double_t get_errLo(const Int_t &idx) const { return _errLo[idx]; }
  Double_t get_errHi(const Int_t &idx) const { return _errHi[idx]; }
  Double_t get_sumw2(const Int_t &idx) const { return _sumw2[idx]; }

  Double_t get_curWeight() const { return _curWeight; }
  Double_t get_curWgtErrLo() const { return _curWgtErrLo; }
  Double_t get_curWgtErrHi() const { return _curWgtErrHi; }
  Double_t get_curSumW2() const { return _curSumW2; }

  Int_t get_curIndex() const { return _curIndex; }

  Int_t       _arrSize ; //  Size of the weight array
  std::vector<Int_t> _idxMult ; // Multiplier jump table for index calculation

  Double_t*       _wgt ; //[_arrSize] Weight array
  Double_t*     _errLo ; //[_arrSize] Low-side error on weight array
  Double_t*     _errHi ; //[_arrSize] High-side error on weight array
  Double_t*     _sumw2 ; //[_arrSize] Sum of weights^2
  Double_t*      _binv ; //[_arrSize] Bin volume array  

  RooArgSet  _realVars ; // Real dimensions of the dataset 
  Bool_t*    _binValid ; //! Valid bins with current range definition
 
  mutable Double_t _curWeight{0.}; // Weight associated with the current coordinate
  mutable Double_t _curWgtErrLo{0.}; // Error on weight associated with the current coordinate
  mutable Double_t _curWgtErrHi{0.}; // Error on weight associated with the current coordinate
  mutable Double_t _curSumW2{0.}; // Current sum of weights^2
  mutable Double_t _curVolume{0.}; // Volume of bin enclosing current coordinate
  mutable Int_t    _curIndex{0}; // Current index

  mutable std::vector<Double_t>* _pbinv ; //! Partial bin volume array
  mutable RooCacheManager<std::vector<Double_t> > _pbinvCacheMgr ; //! Cache manager for arrays of partial bin volumes
  std::vector<RooAbsLValue*> _lvvars ; //! List of observables casted as RooAbsLValue
  std::vector<const RooAbsBinning*> _lvbins ; //! List of used binnings associated with lvalues
  mutable std::vector<std::vector<Double_t> > _binbounds; //! list of bin bounds per dimension

  mutable Int_t _cache_sum_valid{0}; //! Is cache sum valid
  mutable Double_t _cache_sum{0.}; //! Cache for sum of entries ;


private:
  void _adjustBinning(RooRealVar &theirVar, const TAxis &axis, RooRealVar *ourVar, Int_t *offset);

  ClassDef(RooDataHist,4) // Binned data set
};

#endif

