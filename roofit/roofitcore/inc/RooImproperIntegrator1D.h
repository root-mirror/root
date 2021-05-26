/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooImproperIntegrator1D.h,v 1.12 2007/05/11 09:11:30 verkerke Exp $
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
#ifndef ROO_IMPROPER_INTEGRATOR_1D
#define ROO_IMPROPER_INTEGRATOR_1D

#include "RooAbsIntegrator.h"
#include "RooNumIntConfig.h"

class RooInvTransform;
class RooIntegrator1D;

class RooImproperIntegrator1D : public RooAbsIntegrator {
public:

  RooImproperIntegrator1D() ;
  RooImproperIntegrator1D(const RooAbsFunc& function);
  RooImproperIntegrator1D(const RooAbsFunc& function, const RooNumIntConfig& config);
  RooImproperIntegrator1D(const RooAbsFunc& function, Double_t xmin, Double_t xmax, const RooNumIntConfig& config);
  virtual RooAbsIntegrator* clone(const RooAbsFunc& function, const RooNumIntConfig& config) const ;
  virtual ~RooImproperIntegrator1D();

  virtual Bool_t checkLimits() const;
  using RooAbsIntegrator::setLimits ;
  Bool_t setLimits(Double_t* xmin, Double_t* xmax);
  virtual Bool_t setUseIntegrandLimits(Bool_t flag) {_useIntegrandLimits = flag ; return kTRUE ; }
  virtual Double_t integral(const Double_t* yvec=0) ;

  virtual Bool_t canIntegrate1D() const { return kTRUE ; }
  virtual Bool_t canIntegrate2D() const { return kFALSE ; }
  virtual Bool_t canIntegrateND() const { return kFALSE ; }
  virtual Bool_t canIntegrateOpenEnded() const { return kTRUE ; }

protected:

  friend class RooNumIntFactory ;
  static void registerIntegrator(RooNumIntFactory& fact) ;

  void initialize(const RooAbsFunc* function=0) ;

  enum LimitsCase { Invalid, ClosedBothEnds, OpenBothEnds, OpenBelowSpansZero, OpenBelow,
          OpenAboveSpansZero, OpenAbove };
  LimitsCase limitsCase() const;
  LimitsCase _case; // Configuration of limits
  mutable Double_t _xmin, _xmax;  // Value of limits
  Bool_t _useIntegrandLimits;  // Use limits in function binding?

  RooAbsFunc*      _origFunc ;  // Original function binding
  RooInvTransform *_function;   // Binding with inverse of function
  RooNumIntConfig  _config ;  // Configuration object
  mutable RooIntegrator1D *_integrator1,*_integrator2,*_integrator3; // Piece integrators

  ClassDef(RooImproperIntegrator1D,0) // 1-dimensional improper integration engine
};

#endif
