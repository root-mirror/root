/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooAbsOptGoodnessOfFit.h,v 1.15 2007/05/11 09:11:30 verkerke Exp $
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
#ifndef ROO_ABS_OPT_TEST_STATISTIC
#define ROO_ABS_OPT_TEST_STATISTIC

#include "RooAbsTestStatistic.h"
#include "RooSetProxy.h"
#include "RooCategoryProxy.h"
#include "TString.h"

class RooArgSet ;
class RooAbsData ;
class RooAbsReal ;

class RooAbsOptTestStatistic : public RooAbsTestStatistic {
public:

  // Constructors, assignment etc
  RooAbsOptTestStatistic() ;
  RooAbsOptTestStatistic(const char *name, const char *title, RooAbsReal& real, RooAbsData& data,
                         const RooArgSet& projDeps,
                         RooAbsTestStatistic::Configuration const& cfg);
  RooAbsOptTestStatistic(const RooAbsOptTestStatistic& other, const char* name=0);
  virtual ~RooAbsOptTestStatistic();

  virtual Double_t combinedValue(RooAbsReal** gofArray, Int_t nVal) const ;

  RooAbsReal& function() { return *_funcClone ; }
  const RooAbsReal& function() const { return *_funcClone ; }

  RooAbsData& data() ;
  const RooAbsData& data() const ;


  virtual const char* cacheUniqueSuffix() const { return Form("_%lx", (ULong_t)_dataClone) ; }

  // Override this to be always true to force calculation of likelihood without parameters
  virtual Bool_t isDerived() const { return kTRUE ; }

  void seal(const char* notice="") { _sealed = kTRUE ; _sealNotice = notice ; }
  Bool_t isSealed() const { return _sealed ; }
  const char* sealNotice() const { return _sealNotice.Data() ; }

private:
  void setUpBinSampling();

protected:

  Bool_t setDataSlave(RooAbsData& data, Bool_t cloneData=kTRUE, Bool_t ownNewDataAnyway=kFALSE) ;
  void initSlave(RooAbsReal& real, RooAbsData& indata, const RooArgSet& projDeps, const char* rangeName,
       const char* addCoefRangeName)  ;

  friend class RooAbsReal ;
  friend class RooAbsTestStatistic ;

  virtual Bool_t allowFunctionCache() { return kTRUE ;  }
  void constOptimizeTestStatistic(ConstOpCode opcode, Bool_t doAlsoTrackingOpt=kTRUE) ;

  virtual Bool_t redirectServersHook(const RooAbsCollection& newServerList, Bool_t mustReplaceAll, Bool_t nameChange, Bool_t isRecursive) ;
  virtual void printCompactTreeHook(std::ostream& os, const char* indent="") ;
  virtual RooArgSet requiredExtraObservables() const { return RooArgSet() ; }
  void optimizeCaching() ;
  void optimizeConstantTerms(Bool_t,Bool_t=kTRUE) ;

  RooArgSet*  _normSet ; // Pointer to set with observables used for normalization
  RooArgSet*  _funcCloneSet ; // Set owning all components of internal clone of input function
  RooAbsData* _dataClone{nullptr}; // Pointer to internal clone if input data
  RooAbsReal* _funcClone ; // Pointer to internal clone of input function
  RooArgSet*  _projDeps ; // Set of projected observable
  Bool_t      _ownData  ; // Do we own the dataset
  Bool_t      _sealed ; // Is test statistic sealed -- i.e. no access to data
  TString     _sealNotice ; // User-defined notice shown when reading a sealed likelihood
  RooArgSet*  _funcObsSet ; // List of observables in the pdf expression
  RooArgSet   _cachedNodes ; //! List of nodes that are cached as constant expressions

  RooAbsReal* _origFunc ; // Original function
  RooAbsData* _origData ; // Original data
  Bool_t      _optimized ; //!
  double      _integrateBinsPrecision{-1.}; // Precision for finer sampling of bins.

  ClassDef(RooAbsOptTestStatistic,5) // Abstract base class for optimized test statistics
};

#endif
