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
#ifndef ROO_FOAM_GENERATOR
#define ROO_FOAM_GENERATOR

#include "RooAbsNumGenerator.h"
#include "RooPrintable.h"
#include "RooArgSet.h"

class RooAbsReal;
class RooRealVar;
class RooDataSet;

class TFoam ;
class RooTFoamBinding ;
class RooNumGenFactory ;

class RooFoamGenerator : public RooAbsNumGenerator {
public:
  RooFoamGenerator() : _binding(0), _tfoam(0), _xmin(0), _range(0), _vec(0) {} ;
  RooFoamGenerator(const RooAbsReal &func, const RooArgSet &genVars, const RooNumGenConfig& config, Bool_t verbose=kFALSE, const RooAbsReal* maxFuncVal=0);
  RooAbsNumGenerator* clone(const RooAbsReal& func, const RooArgSet& genVars, const RooArgSet& /*condVars*/,
             const RooNumGenConfig& config, Bool_t verbose=kFALSE, const RooAbsReal* maxFuncVal=0) const {
    return new RooFoamGenerator(func,genVars,config,verbose,maxFuncVal) ;
  }
  virtual ~RooFoamGenerator();

  const RooArgSet *generateEvent(UInt_t remaining, Double_t& resampleRatio);

  TFoam& engine() { return *_tfoam; }

  virtual Bool_t canSampleConditional() const { return kFALSE ; }
  virtual Bool_t canSampleCategories() const { return kFALSE ; }

protected:

  friend class RooNumGenFactory ;
  static void registerSampler(RooNumGenFactory& fact) ;

  RooTFoamBinding* _binding ; // Binding of RooAbsReal to TFoam function interface
  TFoam*           _tfoam ;   // Instance of TFOAM generator
  Double_t*        _xmin ;    // Lower bound of observables to be generated ;
  Double_t*        _range ;   // Range of observables to be generated ;
  Double_t*        _vec ;     // Transfer array for FOAM output


  ClassDef(RooFoamGenerator,0) // Context for generating a dataset from a PDF using the TFoam class
};

#endif
