/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: PiecewiseInterpolation.h,v 1.3 2007/05/11 09:11:30 verkerke Exp $
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
#ifndef ROO_PIECEWISEINTERPOLATION
#define ROO_PIECEWISEINTERPOLATION

#include "RooAbsReal.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"

#include "RooObjCacheManager.h"
#include <vector>
#include <list>

class RooRealVar;
class RooArgList;

class PiecewiseInterpolation : public RooAbsReal {
public:

  PiecewiseInterpolation() ;
  PiecewiseInterpolation(const char *name, const char *title, const RooAbsReal& nominal, const RooArgList& lowSet, const RooArgList& highSet, const RooArgList& paramSet, Bool_t takeOwnerShip=kFALSE) ;
  virtual ~PiecewiseInterpolation() ;

  PiecewiseInterpolation(const PiecewiseInterpolation& other, const char* name = 0);
  virtual TObject* clone(const char* newname) const { return new PiecewiseInterpolation(*this, newname); }

  //  virtual Double_t defaultErrorLevel() const ;

  //  void printMetaArgs(std::ostream& os) const ;

  const RooArgList& lowList() const { return _lowSet ; }
  const RooArgList& highList() const { return _highSet ; }
  const RooArgList& paramList() const { return _paramSet ; }

  //virtual Bool_t forceAnalyticalInt(const RooAbsArg&) const { return kTRUE ; }
  Bool_t setBinIntegrator(RooArgSet& allVars) ;

  Int_t getAnalyticalIntegralWN(RooArgSet& allVars, RooArgSet& analVars, const RooArgSet* normSet,const char* rangeName=0) const ;
  Double_t analyticalIntegralWN(Int_t code, const RooArgSet* normSet, const char* rangeName=0) const ;

  void setPositiveDefinite(bool flag=true){_positiveDefinite=flag;}

  void setInterpCode(RooAbsReal& param, int code);
  void setAllInterpCodes(int code);
  void printAllInterpCodes();

  virtual std::list<Double_t>* binBoundaries(RooAbsRealLValue& /*obs*/, Double_t /*xlo*/, Double_t /*xhi*/) const ;
  virtual std::list<Double_t>* plotSamplingHint(RooAbsRealLValue& obs, Double_t xlo, Double_t xhi) const ; 
  virtual Bool_t isBinnedDistribution(const RooArgSet& obs) const ;

protected:

  class CacheElem : public RooAbsCacheElement {
  public:
    CacheElem()  {} ;
    virtual ~CacheElem() {} ; 
    virtual RooArgList containedArgs(Action) { 
      RooArgList ret(_funcIntList) ; 
      ret.add(_lowIntList); 
      ret.add(_highIntList);
      return ret ; 
    }
    RooArgList _funcIntList ;
    RooArgList _lowIntList ;
    RooArgList _highIntList ;
    // will want std::vector<RooRealVar*> for low and high also
  } ;
  mutable RooObjCacheManager _normIntMgr ; // The integration cache manager

  RooRealProxy _nominal;           // The nominal value
  RooArgList   _ownedList ;       // List of owned components
  RooListProxy _lowSet ;            // Low-side variation
  RooListProxy _highSet ;            // High-side varaition
  RooListProxy _paramSet ;            // interpolation parameters
  RooListProxy _normSet ;            // interpolation parameters
  Bool_t _positiveDefinite; // protect against negative and 0 bins.

  std::vector<int> _interpCode;

  Double_t evaluate() const;

  ClassDef(PiecewiseInterpolation,3) // Sum of RooAbsReal objects
};

#endif
