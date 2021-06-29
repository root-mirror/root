/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooRealVar.h,v 1.54 2007/05/11 09:11:30 verkerke Exp $
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
#ifndef ROO_REAL_VAR
#define ROO_REAL_VAR

#include "RooAbsRealLValue.h"

#include "TString.h"

#include <list>
#include <string>
#include <map>
#include <memory>
#include <unordered_map>


class RooArgSet ;
class RooErrorVar ;
class RooVectorDataStore ;
class RooExpensiveObjectCache ;
class RooRealVarSharedProperties;
namespace RooBatchCompute{
struct RunContext;
}

class RooRealVar : public RooAbsRealLValue {
public:
  // Constructors, assignment etc.
  RooRealVar() ;
  RooRealVar(const char *name, const char *title,
  	   Double_t value, const char *unit= "") ;
  RooRealVar(const char *name, const char *title, Double_t minValue, 
	   Double_t maxValue, const char *unit= "");
  RooRealVar(const char *name, const char *title, Double_t value, 
	   Double_t minValue, Double_t maxValue, const char *unit= "") ;
  RooRealVar(const RooRealVar& other, const char* name=0);
  RooRealVar& operator=(const RooRealVar& other);
  virtual TObject* clone(const char* newname) const { return new RooRealVar(*this,newname); }
  virtual ~RooRealVar();
  
  // Parameter value and error accessors
  virtual Double_t getValV(const RooArgSet* nset=0) const ;
  RooSpan<const double> getValues(RooBatchCompute::RunContext& inputData, const RooArgSet* = nullptr) const final;

  virtual void setVal(Double_t value);
  virtual void setVal(Double_t value, const char* rangeName);
  inline Double_t getError() const { return _error>=0?_error:0. ; }
  inline Bool_t hasError(Bool_t allowZero=kTRUE) const { return allowZero ? (_error>=0) : (_error>0) ; }
  inline void setError(Double_t value) { _error= value ; }
  inline void removeError() { _error = -1 ; }
  inline Double_t getAsymErrorLo() const { return _asymErrLo<=0?_asymErrLo:0. ; }
  inline Double_t getAsymErrorHi() const { return _asymErrHi>=0?_asymErrHi:0. ; }
  inline Bool_t hasAsymError(Bool_t allowZero=kTRUE) const { return allowZero ? ((_asymErrHi>=0 && _asymErrLo<=0)) :  ((_asymErrHi>0 && _asymErrLo<0)) ; }
  inline void removeAsymError() { _asymErrLo = 1 ; _asymErrHi = -1 ; }
  inline void setAsymError(Double_t lo, Double_t hi) { _asymErrLo = lo ; _asymErrHi = hi ; }
  inline Double_t getErrorLo() const { return _asymErrLo<=0?_asymErrLo:-1*_error ; }
  inline Double_t getErrorHi() const { return _asymErrHi>=0?_asymErrHi:_error ; }
  
  RooErrorVar* errorVar() const ;

  // Set/get finite fit range limits
  void setMin(const char* name, Double_t value) ;
  void setMax(const char* name, Double_t value) ;
  void setRange(const char* name, Double_t min, Double_t max) ;
  void setRange(const char* name, RooAbsReal& min, RooAbsReal& max) ;
  inline void setMin(Double_t value) { setMin(0,value) ; }
  inline void setMax(Double_t value) { setMax(0,value) ; }
  /// Set the limits of the default range.
  inline void setRange(Double_t min, Double_t max) { setRange(0,min,max) ; }
  /// Set parameterised limits of the default range. See setRange(const char*, RooAbsReal&, RooAbsReal&).
  inline void setRange(RooAbsReal& min, RooAbsReal& max) { setRange(0,min,max) ; }

  void setBins(Int_t nBins, const char* name=0);
  void setBinning(const RooAbsBinning& binning, const char* name=0) ;

  // RooAbsRealLValue implementation
  Bool_t hasBinning(const char* name) const ;
  const RooAbsBinning& getBinning(const char* name=0, Bool_t verbose=kTRUE, Bool_t createOnTheFly=kFALSE) const ;
  RooAbsBinning& getBinning(const char* name=0, Bool_t verbose=kTRUE, Bool_t createOnTheFly=kFALSE) ; 
  std::list<std::string> getBinningNames() const ;

  // Set infinite fit range limits
  /// Remove lower range limit for binning with given name. Empty name means default range.
  void removeMin(const char* name=0);
  /// Remove upper range limit for binning with given name. Empty name means default range.
  void removeMax(const char* name=0);
  /// Remove range limits for binning with given name. Empty name means default range.
  void removeRange(const char* name=0);
 
  // I/O streaming interface (machine readable)
  virtual Bool_t readFromStream(std::istream& is, Bool_t compact, Bool_t verbose=kFALSE) ;
  virtual void writeToStream(std::ostream& os, Bool_t compact) const ;

  // We implement a fundamental type of AbsArg that can be stored in a dataset
  inline virtual Bool_t isFundamental() const { return kTRUE; }

  // Force to be a leaf-node of any expression tree, even if we have (shape) servers
  virtual Bool_t isDerived() const { 
    // Does value or shape of this arg depend on any other arg?
    return !_serverList.empty() || _proxyList.GetEntries()>0;
  }

  // Printing interface (human readable)
  virtual void printValue(std::ostream& os) const ;
  virtual void printExtras(std::ostream& os) const ;
  virtual void printMultiline(std::ostream& os, Int_t contents, Bool_t verbose=kFALSE, TString indent="") const ;
  virtual Int_t defaultPrintContents(Option_t* opt) const ;


  TString* format(const RooCmdArg& formatArg) const ;
  TString* format(Int_t sigDigits, const char *options) const ;

  static void printScientific(Bool_t flag=kFALSE) ;
  static void printSigDigits(Int_t ndig=5) ;

  using RooAbsRealLValue::operator= ;

  void deleteSharedProperties() ;

  void copyCacheFast(const RooRealVar& other, Bool_t setValDirty=kTRUE) { _value = other._value ; if (setValDirty) setValueDirty() ; }

  protected:

  static Bool_t _printScientific ;
  static Int_t  _printSigDigits ;

  friend class RooAbsRealLValue ;
  virtual void setValFast(Double_t value) { _value = value ; setValueDirty() ; }


  virtual Double_t evaluate() const { return _value ; } // dummy because we overloaded getVal()
  virtual void copyCache(const RooAbsArg* source, Bool_t valueOnly=kFALSE, Bool_t setValDirty=kTRUE) ;
  virtual void attachToTree(TTree& t, Int_t bufSize=32000) ;
  virtual void attachToVStore(RooVectorDataStore& vstore) ;
  virtual void fillTreeBranch(TTree& t) ;

  Double_t chopAt(Double_t what, Int_t where) const ;

  Double_t _error;      // Symmetric error associated with current value
  Double_t _asymErrLo ; // Low side of asymmetric error associated with current value
  Double_t _asymErrHi ; // High side of asymmetric error associated with current value
  std::unique_ptr<RooAbsBinning> _binning;
  std::unordered_map<std::string,std::unique_ptr<RooAbsBinning>> _altNonSharedBinning ; //! Non-shareable alternative binnings

  std::shared_ptr<RooRealVarSharedProperties> sharedProp() const;
  void installSharedProp(std::shared_ptr<RooRealVarSharedProperties>&& prop);

  virtual void setExpensiveObjectCache(RooExpensiveObjectCache&) { ; } // variables don't need caches 
  static RooRealVarSharedProperties& _nullProp(); // Null property
  static std::map<std::string,std::weak_ptr<RooRealVarSharedProperties>>& _sharedPropList(); // List of properties shared among clones of a variable
  
  std::shared_ptr<RooRealVarSharedProperties> _sharedProp; //! Shared binnings associated with this instance

  ClassDef(RooRealVar,7) // Real-valued variable
};




#endif
