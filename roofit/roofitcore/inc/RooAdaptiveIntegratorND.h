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
#ifndef ROO_ADAPTIVE_INTEGRATOR_ND
#define ROO_ADAPTIVE_INTEGRATOR_ND

#include "RooAbsIntegrator.h"
#include "RooNumIntConfig.h"
#include "TString.h"

namespace ROOT { namespace Math { class AdaptiveIntegratorMultiDim ; } }
class RooMultiGenFunction ;

class RooAdaptiveIntegratorND : public RooAbsIntegrator {
public:

  // Constructors, assignment etc
  RooAdaptiveIntegratorND() ;
  RooAdaptiveIntegratorND(const RooAbsFunc& function, const RooNumIntConfig& config) ;

  virtual RooAbsIntegrator* clone(const RooAbsFunc& function, const RooNumIntConfig& config) const ;
  virtual ~RooAdaptiveIntegratorND();

  virtual Bool_t checkLimits() const;
  virtual Double_t integral(const Double_t *yvec=0) ;

  using RooAbsIntegrator::setLimits ;
  Bool_t setLimits(Double_t* xmin, Double_t* xmax);

  virtual Bool_t canIntegrate1D() const { return kFALSE ; }
  virtual Bool_t canIntegrate2D() const { return kTRUE ; }
  virtual Bool_t canIntegrateND() const { return kTRUE ; }
  virtual Bool_t canIntegrateOpenEnded() const { return kFALSE ; }

  virtual Bool_t setUseIntegrandLimits(Bool_t flag) {_useIntegrandLimits = flag ; return kTRUE ; }

protected:

  RooAdaptiveIntegratorND(const RooAdaptiveIntegratorND&) ;

  Bool_t _useIntegrandLimits;  // If true limits of function binding are used

  mutable Double_t* _xmin ;    // Lower bound in each dimension
  mutable Double_t* _xmax ;    // Upper bound in each dimension
  Double_t _epsRel ;           // Relative precision
  Double_t _epsAbs ;           // Absolute precision
  Int_t    _nmax ;             // Max number of divisions
  Int_t    _nError ;           // Number of error occurrences
  Int_t    _nWarn ;            // Max number of warnings to be issued ;
  RooMultiGenFunction* _func ; //! ROOT::Math multi-parameter function binding
  ROOT::Math::AdaptiveIntegratorMultiDim* _integrator ;
  TString _intName ; // Integrand name

  friend class RooNumIntFactory ;
  static void registerIntegrator(RooNumIntFactory& fact) ;

  ClassDef(RooAdaptiveIntegratorND,0) // N-dimensional adaptive integration (interface to MathCore integrator)
};

#endif
