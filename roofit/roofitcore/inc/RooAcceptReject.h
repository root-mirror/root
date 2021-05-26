/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooAcceptReject.h,v 1.16 2007/05/11 09:11:30 verkerke Exp $
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
#ifndef ROO_ACCEPT_REJECT
#define ROO_ACCEPT_REJECT

#include "RooAbsNumGenerator.h"
#include "RooPrintable.h"
#include "RooArgSet.h"

class RooAbsReal;
class RooRealVar;
class RooDataSet;
class RooRealBinding;
class RooNumGenFactory ;

class RooAcceptReject : public RooAbsNumGenerator {
public:
  RooAcceptReject() : _nextCatVar(0), _nextRealVar(0) {
    // coverity[UNINIT_CTOR]
  } ;
  RooAcceptReject(const RooAbsReal &func, const RooArgSet &genVars, const RooNumGenConfig& config, Bool_t verbose=kFALSE, const RooAbsReal* maxFuncVal=0);
  RooAbsNumGenerator* clone(const RooAbsReal& func, const RooArgSet& genVars, const RooArgSet& /*condVars*/,
			    const RooNumGenConfig& config, Bool_t verbose=kFALSE, const RooAbsReal* maxFuncVal=0) const {
    return new RooAcceptReject(func,genVars,config,verbose,maxFuncVal) ;
  }
  virtual ~RooAcceptReject();

  const RooArgSet *generateEvent(UInt_t remaining, Double_t& resampleRatio);
  Double_t getFuncMax() ;


  // Advertisement of capabilities
  virtual Bool_t canSampleConditional() const { return kTRUE ; }
  virtual Bool_t canSampleCategories() const { return kTRUE ; }


protected:

  friend class RooNumGenFactory ;
  static void registerSampler(RooNumGenFactory& fact) ;

  void addEventToCache();
  const RooArgSet *nextAcceptedEvent();

  Double_t _maxFuncVal, _funcSum;      // Maximum function value found, and sum of all samples made
  UInt_t _realSampleDim,_catSampleMult;// Number of real and discrete dimensions to be sampled
  UInt_t _minTrials;                   // Minimum number of max.finding trials, total number of samples
  UInt_t _totalEvents;                 // Total number of function samples
  UInt_t _eventsUsed;                  // Accepted number of function samples
  TIterator *_nextCatVar;              // Iterator of categories to be generated
  TIterator *_nextRealVar;             // Iterator over variables to be generated

  UInt_t _minTrialsArray[4];           // Minimum number of trials samples for 1,2,3 dimensional problems

  ClassDef(RooAcceptReject,0) // Context for generating a dataset from a PDF
};

#endif
