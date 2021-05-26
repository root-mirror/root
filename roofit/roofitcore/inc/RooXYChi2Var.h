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

#ifndef ROO_XY_CHI2_VAR
#define ROO_XY_CHI2_VAR

#include "RooAbsOptTestStatistic.h"
#include "RooCmdArg.h"
#include "RooDataSet.h"
#include "RooAbsPdf.h"
#include "RooNumIntConfig.h"
#include <list>
class RooAbsIntegrator ;


class RooXYChi2Var : public RooAbsOptTestStatistic {
public:

  // Constructors, assignment etc
  RooXYChi2Var() ;
  RooXYChi2Var(const char *name, const char* title, RooAbsReal& func, RooDataSet& data, Bool_t integrate=kFALSE) ;
  RooXYChi2Var(const char *name, const char* title, RooAbsReal& func, RooDataSet& data, RooRealVar& yvar, Bool_t integrate=kFALSE) ;
  RooXYChi2Var(const char *name, const char* title, RooAbsPdf& extPdf, RooDataSet& data, Bool_t integrate=kFALSE) ;
  RooXYChi2Var(const char *name, const char* title, RooAbsPdf& extPdf, RooDataSet& data, RooRealVar& yvar, Bool_t integrate=kFALSE) ;

  RooXYChi2Var(const RooXYChi2Var& other, const char* name=0);
  virtual TObject* clone(const char* newname) const { return new RooXYChi2Var(*this,newname); }

  virtual RooAbsTestStatistic* create(const char *name, const char *title, RooAbsReal& pdf, RooAbsData& adata,
                                      const RooArgSet&, RooAbsTestStatistic::Configuration const&) {
    // Virtual constructor
    return new RooXYChi2Var(name,title,pdf,(RooDataSet&)adata) ;
  }

  virtual ~RooXYChi2Var();

  virtual Double_t defaultErrorLevel() const {
    // The default error level for MINUIT error analysis for a chi^2 is 1.0
    return 1.0 ;
  }

  RooNumIntConfig& binIntegratorConfig() { return _intConfig ; }
  const RooNumIntConfig& binIntegratorConfig() const { return _intConfig ; }

protected:

  Bool_t allowFunctionCache() {
    // Disable function (component) caching if integration is requested as the function
    // will be evaluated at coordinates other than the points in the dataset
    return !_integrate ;
  }

  RooArgSet requiredExtraObservables() const ;

  Double_t fy() const ;

  Bool_t _extended ; // Is the input function and extended p.d.f.
  Bool_t _integrate ; // Is integration over the bin volume requested

  RooRealVar* _yvar ; // Y variable if so designated
  RooArgSet _rrvArgs ; // Set of real-valued observables
  TIterator* _rrvIter ; //! Iterator over set of real-valued observables

  void initialize() ;
  void initIntegrator() ;
  Double_t xErrorContribution(Double_t ydata) const ;

  virtual Double_t evaluatePartition(std::size_t firstEvent, std::size_t lastEvent, std::size_t stepSize) const ;

  RooNumIntConfig   _intConfig ; // Numeric integrator configuration for integration of function over bin
  RooAbsReal*       _funcInt ; //! Function integral
  std::list<RooAbsBinning*> _binList ; //! Bin ranges

  ClassDef(RooXYChi2Var,1) // Chi^2 function of p.d.f w.r.t a unbinned dataset with X and Y values
};


#endif
