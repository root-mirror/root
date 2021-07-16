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
#ifndef ROO_ABS_DATA_STORE
#define ROO_ABS_DATA_STORE

#include "Rtypes.h"
#include "RooArgSet.h"
#include "RooAbsData.h"

#include "ROOT/RStringView.hxx"
#include "TNamed.h"

#include <list>
#include <vector>

class RooAbsArg ;
class RooArgList ;
class TTree ;
namespace RooBatchCompute {
struct RunContext;
}


class RooAbsDataStore : public TNamed, public RooPrintable {
public:

  RooAbsDataStore() {}
  RooAbsDataStore(std::string_view name, std::string_view title, const RooArgSet& vars)
    : TNamed(TString{name},TString{title}), _vars{vars} {}
  RooAbsDataStore(const RooAbsDataStore& other, const char* newname=0)
    : RooAbsDataStore(other, other._vars, newname) {}
  RooAbsDataStore(const RooAbsDataStore& other, const RooArgSet& vars, const char* newname=0)
    : TNamed(other), RooPrintable(other), _vars{vars}, _doDirtyProp{other._doDirtyProp}
  {
    if(newname) SetName(newname);
  }

  WRITE_TSTRING_COMPATIBLE_CONSTRUCTOR(RooAbsDataStore)

  virtual RooAbsDataStore* clone(const char* newname=0) const = 0 ;
  virtual RooAbsDataStore* clone(const RooArgSet& vars, const char* newname=0) const = 0 ;

  // Write current row
  virtual Int_t fill() = 0 ;
  
  // Retrieve a row
  virtual const RooArgSet* get(Int_t index) const = 0 ;
  virtual const RooArgSet* get() const { return &_vars ; } 
  virtual Double_t weight() const = 0 ;

  virtual Double_t weightError(RooAbsData::ErrorType etype=RooAbsData::Poisson) const = 0 ;
  virtual void weightError(Double_t& lo, Double_t& hi, RooAbsData::ErrorType etype=RooAbsData::Poisson) const = 0 ; 

  double weight(Int_t index) const {
    get(index) ;
    return weight() ;
  }
  virtual Bool_t isWeighted() const = 0 ;

  /// Retrieve batches for all observables in this data store.
  virtual RooBatchCompute::RunContext getBatches(std::size_t first, std::size_t len) const = 0;
  virtual RooSpan<const double> getWeightBatch(std::size_t first, std::size_t len) const = 0;

  // Change observable name
  virtual Bool_t changeObservableName(const char* from, const char* to) =0 ;
  
  // Add one or more columns
  virtual RooAbsArg* addColumn(RooAbsArg& var, Bool_t adjustRange=kTRUE) = 0 ;
  virtual RooArgSet* addColumns(const RooArgList& varList) = 0 ;

  // Merge column-wise
  virtual RooAbsDataStore* merge(const RooArgSet& allvars, std::list<RooAbsDataStore*> dstoreList) = 0 ;

  // Add rows 
  virtual void append(RooAbsDataStore& other)= 0 ;

  // General & bookkeeping methods
  virtual Int_t numEntries() const = 0 ;
  virtual Double_t sumEntries() const { return 0 ; } ;
  virtual void reset() = 0 ;

  // Buffer redirection routines used in inside RooAbsOptTestStatistics
  virtual void attachBuffers(const RooArgSet& extObs) = 0 ; 
  virtual void resetBuffers() = 0 ;

  virtual void setExternalWeightArray(const Double_t* /*arrayWgt*/, const Double_t* /*arrayWgtErrLo*/, const Double_t* /*arrayWgtErrHi*/, const Double_t* /*arraySumW2*/) {} ;

  // Printing interface (human readable)
  inline void Print(Option_t *options= 0) const override {
    // Print contents on stdout
    printStream(defaultPrintStream(),defaultPrintContents(options),defaultPrintStyle(options));
  }

  /// Print name of dataset
  void printName(std::ostream& os) const override { os << GetName(); }
  /// Print title of dataset
  void printTitle(std::ostream& os) const override { os << GetTitle(); }
  void printClassName(std::ostream& os) const override ;
  void printArgs(std::ostream& os) const override;
  /// Print value of the dataset, i.e. the sum of weights contained in the dataset
  void printValue(std::ostream& os) const override { os << numEntries() << " entries" ; }
  void printMultiline(std::ostream& os, Int_t content, Bool_t verbose, TString indent) const override;

  /// Define default print options, for a given print style
  virtual int defaultPrintContents(Option_t* /*opt*/) const override { return kName|kClassName|kArgs|kValue ; }
   

  // Constant term  optimizer interface
  virtual void cacheArgs(const RooAbsArg* cacheOwner, RooArgSet& varSet, const RooArgSet* nset=0, Bool_t skipZeroWeights=kFALSE) = 0 ;
  virtual const RooAbsArg* cacheOwner() = 0 ;
  virtual void attachCache(const RooAbsArg* newOwner, const RooArgSet& cachedVars) = 0 ;
  virtual void setArgStatus(const RooArgSet& set, Bool_t active) = 0 ;
  const RooArgSet& cachedVars() const { return _cachedVars ; }
  virtual void resetCache() = 0 ;
  virtual void recalculateCache(const RooArgSet* /*proj*/, Int_t /*firstEvent*/, Int_t /*lastEvent*/, Int_t /*stepSize*/, Bool_t /* skipZeroWeights*/) {} ;

  virtual void setDirtyProp(Bool_t flag) { _doDirtyProp = flag ; }
  Bool_t dirtyProp() const { return _doDirtyProp ; }

  virtual void checkInit() const {} ;
  
  virtual Bool_t hasFilledCache() const { return kFALSE ; }
  
  virtual const TTree* tree() const { return 0 ; }
  virtual void dump() {} 

  virtual void loadValues(const RooAbsDataStore *tds, const RooFormulaVar* select=0, const char* rangeName=0,
      std::size_t nStart=0, std::size_t nStop = std::numeric_limits<std::size_t>::max()) = 0 ;

  virtual void forceCacheUpdate() {} ;
  
 protected:

  RooArgSet _vars;
  RooArgSet _cachedVars;

  Bool_t _doDirtyProp = true; // Switch do (de)activate dirty state propagation when loading a data point

  ClassDefOverride(RooAbsDataStore,1) // Abstract Data Storage class
};


#endif
