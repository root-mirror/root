/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooDLLSignificanceMCSModule.h,v 1.2 2007/05/11 09:11:30 verkerke Exp $
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

#ifndef ROO_DELTA_LL_SIGNIFICANCE_MCS_MODULE
#define ROO_DELTA_LL_SIGNIFICANCE_MCS_MODULE

#include "RooAbsMCStudyModule.h"
#include <string>

class RooDLLSignificanceMCSModule : public RooAbsMCStudyModule {
public:

  RooDLLSignificanceMCSModule(const RooRealVar& param, Double_t nullHypoValue=0) ;
  RooDLLSignificanceMCSModule(const char* parName, Double_t nullHypoValue=0) ;
  RooDLLSignificanceMCSModule(const RooDLLSignificanceMCSModule& other) ;
  virtual ~RooDLLSignificanceMCSModule() ;

  Bool_t initializeInstance() ;

  Bool_t initializeRun(Int_t /*numSamples*/) ;
  RooDataSet* finalizeRun() ;

  Bool_t processAfterFit(Int_t /*sampleNum*/)  ;

private:

  std::string _parName ;  // Name of Nsignal parameter
  RooDataSet* _data ;     // Summary dataset to store results
  RooRealVar* _nll0h ;    // Container variable for NLL result on null hypothesis
  RooRealVar* _dll0h ;    // Container variable for delta NLL
  RooRealVar* _sig0h ;    // Container variable for NLL result with signal
  Double_t    _nullValue ;  // Numeric value of Nsignal parameter representing the null hypothesis

  ClassDef(RooDLLSignificanceMCSModule,0) // MCStudy module to calculate Delta(-logL) significance w.r.t given null hypothesis
} ;


#endif

