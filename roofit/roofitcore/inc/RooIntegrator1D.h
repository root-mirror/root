/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooIntegrator1D.h,v 1.21 2007/05/11 09:11:30 verkerke Exp $
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
#ifndef ROO_INTEGRATOR_1D
#define ROO_INTEGRATOR_1D

#include "RooAbsIntegrator.h"
#include "RooNumIntConfig.h"

class RooIntegrator1D : public RooAbsIntegrator {
public:

  // Constructors, assignment etc
  enum SummationRule { Trapezoid, Midpoint };
  RooIntegrator1D() ;

  RooIntegrator1D(const RooAbsFunc& function, SummationRule rule= Trapezoid,
        Int_t maxSteps= 0, Double_t eps= 0) ;
  RooIntegrator1D(const RooAbsFunc& function, Double_t xmin, Double_t xmax,
        SummationRule rule= Trapezoid, Int_t maxSteps= 0, Double_t eps= 0) ;

  RooIntegrator1D(const RooAbsFunc& function, const RooNumIntConfig& config) ;
  RooIntegrator1D(const RooAbsFunc& function, Double_t xmin, Double_t xmax,
        const RooNumIntConfig& config) ;

  virtual RooAbsIntegrator* clone(const RooAbsFunc& function, const RooNumIntConfig& config) const ;
  virtual ~RooIntegrator1D();

  virtual Bool_t checkLimits() const;
  virtual Double_t integral(const Double_t *yvec=0) ;

  using RooAbsIntegrator::setLimits ;
  Bool_t setLimits(Double_t* xmin, Double_t* xmax);
  virtual Bool_t setUseIntegrandLimits(Bool_t flag) {_useIntegrandLimits = flag ; return kTRUE ; }

  virtual Bool_t canIntegrate1D() const { return kTRUE ; }
  virtual Bool_t canIntegrate2D() const { return kFALSE ; }
  virtual Bool_t canIntegrateND() const { return kFALSE ; }
  virtual Bool_t canIntegrateOpenEnded() const { return kFALSE ; }

protected:

  friend class RooNumIntFactory ;
  static void registerIntegrator(RooNumIntFactory& fact) ;

  Bool_t initialize();

  Bool_t _useIntegrandLimits;  // If true limits of function binding are used

  // Integrator configuration
  SummationRule _rule;
  Int_t _maxSteps ;      // Maximum number of steps
  Int_t _minStepsZero ;  // Minimum number of steps to declare convergence to zero
  Int_t _fixSteps ;      // Fixed number of steps
  Double_t _epsAbs ;     // Absolute convergence tolerance
  Double_t _epsRel ;     // Relative convergence tolerance
  Bool_t _doExtrap ;     // Apply conversion step?
  enum { _nPoints = 5 };

  // Numerical integrator support functions
  Double_t addTrapezoids(Int_t n) ;
  Double_t addMidpoints(Int_t n) ;
  void extrapolate(Int_t n) ;

  // Numerical integrator workspace
  Double_t _xmin;              //! Lower integration bound
  Double_t _xmax;              //! Upper integration bound
  Double_t _range;             //! Size of integration range
  Double_t _extrapValue;               //! Extrapolated value
  Double_t _extrapError;               //! Error on extrapolated value
  Double_t *_h ;                       //! Integrator workspace
  Double_t *_s ;                       //! Integrator workspace
  Double_t *_c ;                       //! Integrator workspace
  Double_t *_d ;                       //! Integrator workspace
  Double_t _savedResult;               //! Integrator workspace

  Double_t* xvec(Double_t& xx) { _x[0] = xx ; return _x ; }

  Double_t *_x ; //! do not persist

  ClassDef(RooIntegrator1D,0) // 1-dimensional numerical integration engine
};

#endif
