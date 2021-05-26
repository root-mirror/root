/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooMCIntegrator.h,v 1.15 2007/05/11 09:11:30 verkerke Exp $
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
#ifndef ROO_MC_INTEGRATOR
#define ROO_MC_INTEGRATOR

#include "RooAbsIntegrator.h"
#include "RooGrid.h"
#include "RooNumIntConfig.h"
#include "TStopwatch.h"

class RooMCIntegrator : public RooAbsIntegrator {
public:

  // Constructors, assignment etc
  enum SamplingMode { Importance, ImportanceOnly, Stratified };
  enum GeneratorType { QuasiRandom, PseudoRandom };
  RooMCIntegrator() ;
  RooMCIntegrator(const RooAbsFunc& function, SamplingMode mode= Importance,
        GeneratorType genType= QuasiRandom, Bool_t verbose= kFALSE);
  RooMCIntegrator(const RooAbsFunc& function, const RooNumIntConfig& config);
  virtual RooAbsIntegrator* clone(const RooAbsFunc& function, const RooNumIntConfig& config) const ;
  virtual ~RooMCIntegrator();

  virtual Bool_t checkLimits() const;
  virtual Double_t integral(const Double_t* yvec=0);

  enum Stage { AllStages, ReuseGrid, RefineGrid };
  Double_t vegas(Stage stage, UInt_t calls, UInt_t iterations, Double_t *absError= 0);

  Double_t getAlpha() const { return _alpha;   }
  void setAlpha(Double_t alpha) { _alpha= alpha; }

  GeneratorType getGenType() const { return _genType; }
  void setGenType(GeneratorType type) { _genType= type; }

  const RooGrid &grid() const { return _grid; }

  virtual Bool_t canIntegrate1D() const { return kTRUE ; }
  virtual Bool_t canIntegrate2D() const { return kTRUE ; }
  virtual Bool_t canIntegrateND() const { return kTRUE ; }
  virtual Bool_t canIntegrateOpenEnded() const { return kFALSE ; }

protected:

  friend class RooNumIntFactory ;
  static void registerIntegrator(RooNumIntFactory& fact) ;

  mutable RooGrid _grid;  // Sampling grid definition

  // control variables
  Bool_t _verbose;          // Verbosity control
  Double_t _alpha;          // Grid stiffness parameter
  Int_t _mode;              // Sampling mode
  GeneratorType _genType;   // Generator type
  Int_t _nRefineIter ;      // Number of refinement iterations
  Int_t _nRefinePerDim ;    // Number of refinement samplings (per dim)
  Int_t _nIntegratePerDim ; // Number of integration samplings (per dim)

  TStopwatch _timer;        // Timer

  Double_t _jac,_wtd_int_sum,_sum_wgts,_chi_sum,_chisq,_result,_sigma; // Scratch variables preserved between calls to vegas1/2/2
  UInt_t _it_start,_it_num,_samples,_calls_per_box;                    // Scratch variables preserved between calls to vegas1/2/2

  ClassDef(RooMCIntegrator,0) // VEGAS based multi-dimensional numerical integration engine
};

#endif
