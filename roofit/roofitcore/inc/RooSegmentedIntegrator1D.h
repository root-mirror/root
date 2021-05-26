/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooSegmentedIntegrator1D.h,v 1.7 2007/05/11 09:11:30 verkerke Exp $
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
#ifndef ROO_SEGMENTED_INTEGRATOR_1D
#define ROO_SEGMENTED_INTEGRATOR_1D

#include "RooAbsIntegrator.h"
#include "RooIntegrator1D.h"
#include "RooNumIntConfig.h"

class RooSegmentedIntegrator1D : public RooAbsIntegrator {
public:

  // Constructors, assignment etc
  RooSegmentedIntegrator1D() ;
  RooSegmentedIntegrator1D(const RooAbsFunc& function, const RooNumIntConfig& config) ;
  RooSegmentedIntegrator1D(const RooAbsFunc& function, Double_t xmin, Double_t xmax, const RooNumIntConfig& config) ;

  virtual RooAbsIntegrator* clone(const RooAbsFunc& function, const RooNumIntConfig& config) const ;
  virtual ~RooSegmentedIntegrator1D();

  virtual Bool_t checkLimits() const;
  virtual Double_t integral(const Double_t *yvec=0) ;

  using RooAbsIntegrator::setLimits ;
  Bool_t setLimits(Double_t *xmin, Double_t *xmax);
  virtual Bool_t setUseIntegrandLimits(Bool_t flag) { _useIntegrandLimits = flag ; return kTRUE ; }

  virtual Bool_t canIntegrate1D() const { return kTRUE ; }
  virtual Bool_t canIntegrate2D() const { return kFALSE ; }
  virtual Bool_t canIntegrateND() const { return kFALSE ; }
  virtual Bool_t canIntegrateOpenEnded() const { return kFALSE ; }

protected:

  friend class RooNumIntFactory ;
  static void registerIntegrator(RooNumIntFactory& fact) ;

  mutable Double_t _xmin ;
  mutable Double_t _xmax ;
  mutable Double_t _range ;
  Bool_t _valid ;
  Int_t _nseg ; // Number of segments
  Bool_t _useIntegrandLimits ;

  RooNumIntConfig _config ;
  RooIntegrator1D** _array ; // Array of segment integrators

  Bool_t initialize();

  ClassDef(RooSegmentedIntegrator1D,0) // 1-dimensional piece-wise numerical integration engine
};

#endif
