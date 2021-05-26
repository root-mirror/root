/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooRealIntegral.h,v 1.44 2007/05/11 09:11:30 verkerke Exp $
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
#ifndef ROO_REAL_INTEGRAL
#define ROO_REAL_INTEGRAL

#include "RooAbsReal.h"
#include "RooArgSet.h"
#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooSetProxy.h"
#include "RooListProxy.h"
#include <list>

class RooArgSet ;
class TH1F ;
class RooAbsCategory ;
class RooRealVar ;
class RooAbsIntegrator ;
class RooNumIntConfig ;

class RooRealIntegral : public RooAbsReal {
public:

  // Constructors, assignment etc
  RooRealIntegral() ;
  RooRealIntegral(const char *name, const char *title, const RooAbsReal& function, const RooArgSet& depList,
        const RooArgSet* funcNormSet=0, const RooNumIntConfig* config=0, const char* rangeName=0) ;
  RooRealIntegral(const RooRealIntegral& other, const char* name=0);
  virtual TObject* clone(const char* newname) const { return new RooRealIntegral(*this,newname); }
  virtual ~RooRealIntegral();

  virtual Double_t getValV(const RooArgSet* set=0) const ;

  Bool_t isValid() const { return _valid; }

  void printMultiline(std::ostream& os, Int_t contents, Bool_t verbose=kFALSE, TString indent="") const ;
  void printMetaArgs(std::ostream& os) const ;

  const RooArgSet& numIntCatVars() const { return _sumList ; }
  const RooArgSet& numIntRealVars() const { return _intList ; }
  const RooArgSet& anaIntVars() const { return _anaList ; }

  RooArgSet intVars() const { RooArgSet tmp(_sumList) ; tmp.add(_intList) ; tmp.add(_anaList) ; tmp.add(_facList) ; return tmp ; }
  const char* intRange() { return _rangeName ? _rangeName->GetName() : 0 ; }
  const RooAbsReal& integrand() const { return _function.arg() ; }

  void setCacheNumeric(Bool_t flag) {
    // If true, value of this integral is cached if it is (partially numeric)
    _cacheNum = flag ;
  }

  Bool_t getCacheNumeric() {
    // If true, value of this integral is cached if it is (partially numeric)
    return _cacheNum ;
  }

  static void setCacheAllNumeric(Int_t ndim) ;

  static Int_t getCacheAllNumeric() ;

  virtual std::list<Double_t>* plotSamplingHint(RooAbsRealLValue& obs, Double_t xlo, Double_t xhi) const {
    // Forward plot sampling hint of integrand
    return _function.arg().plotSamplingHint(obs,xlo,xhi) ;
  }

  virtual RooAbsReal* createIntegral(const RooArgSet& iset, const RooArgSet* nset=0, const RooNumIntConfig* cfg=0, const char* rangeName=0) const ;

  void setAllowComponentSelection(Bool_t allow);
  Bool_t getAllowComponentSelection() const;

protected:

  mutable Bool_t _valid;
  Bool_t _respectCompSelect;

  const RooArgSet& parameters() const ;

  enum IntOperMode { Hybrid, Analytic, PassThrough } ;
  //friend class RooAbsPdf ;

  Bool_t initNumIntegrator() const;
  void autoSelectDirtyMode() ;

  virtual Double_t sum() const ;
  virtual Double_t integrate() const ;
  virtual Double_t jacobianProduct() const ;

  // Evaluation and validation implementation
  Double_t evaluate() const ;
  virtual Bool_t isValidReal(Double_t value, Bool_t printError=kFALSE) const ;
  Bool_t servesExclusively(const RooAbsArg* server,const RooArgSet& exclLVBranches, const RooArgSet& allBranches) const ;


  virtual Bool_t redirectServersHook(const RooAbsCollection& newServerList,
                 Bool_t mustReplaceAll, Bool_t nameChange, Bool_t isRecursive) ;

  // Function pointer and integrands list
  mutable RooSetProxy _sumList ; // Set of discrete observable over which is summed numerically
  mutable RooSetProxy _intList ; // Set of continuous observables over which is integrated numerically
  mutable RooSetProxy _anaList ; // Set of observables over which is integrated/summed analytically
  mutable RooSetProxy _jacList ; // Set of lvalue observables over which is analytically integration that have a non-unit Jacobian
  mutable RooSetProxy _facList ; // Set of observables on which function does not depends, which are integrated nevertheless

  mutable RooArgSet   _facListOwned ;  // Owned components in _facList
  RooRealProxy       _function ;     // Function being integration
  RooArgSet*      _funcNormSet ;     // Optional normalization set passed to function

  mutable RooArgSet       _saveInt ; //! do not persist
  mutable RooArgSet       _saveSum ; //! do not persist

  RooNumIntConfig* _iconfig ;

  mutable RooListProxy _sumCat ; //! do not persist

  Int_t _mode ;
  IntOperMode _intOperMode ;   // integration operation mode

  mutable Bool_t _restartNumIntEngine ; //! do not persist
  mutable RooAbsIntegrator* _numIntEngine ;  //! do not persist
  mutable RooAbsFunc *_numIntegrand;         //! do not persist

  TNamed* _rangeName ;

  mutable RooArgSet* _params ; //! cache for set of parameters

  Bool_t _cacheNum ;           // Cache integral if numeric
  static Int_t _cacheAllNDim ; //! Cache all integrals with given numeric dimension

  ClassDef(RooRealIntegral,3) // Real-valued function representing an integral over a RooAbsReal object
};

#endif
