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
#ifndef ROO_DATA_WEIGHTED_AVERAGE
#define ROO_DATA_WEIGHTED_AVERAGE

#include "RooAbsOptTestStatistic.h"
#include "RooCmdArg.h"

class RooDataWeightedAverage : public RooAbsOptTestStatistic {
public:

  // Constructors, assignment etc
  RooDataWeightedAverage() {
    // Default constructor
  } ;

  RooDataWeightedAverage(const char *name, const char *title, RooAbsReal& real, RooAbsData& data, const RooArgSet& projDeps,
                         RooAbsTestStatistic::Configuration const& cfg, bool showProgress=false) ;

  RooDataWeightedAverage(const RooDataWeightedAverage& other, const char* name=0);
  virtual TObject* clone(const char* newname) const { return new RooDataWeightedAverage(*this,newname); }

  virtual RooAbsTestStatistic* create(const char *name, const char *title, RooAbsReal& real, RooAbsData& adata,
                                      const RooArgSet& projDeps,
                                      RooAbsTestStatistic::Configuration const& cfg) {
    // Virtual constructor
    return new RooDataWeightedAverage(name,title,real,adata,projDeps,cfg) ;
  }

  virtual Double_t globalNormalization() const ;

  virtual ~RooDataWeightedAverage();


protected:

  Double_t _sumWeight ;  // Global sum of weights needed for normalization
  Bool_t _showProgress ; // Show progress indication during evaluation if true
  virtual Double_t evaluatePartition(std::size_t firstEvent, std::size_t lastEvent, std::size_t stepSize) const ;

  ClassDef(RooDataWeightedAverage,1) // Optimized calculator of data weighted average of a RooAbsReal
};

#endif
