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
#ifndef ROO_STUDY_MANAGER
#define ROO_STUDY_MANAGER

#include "TNamed.h"

class RooAbsPdf;
class RooDataSet ;
class RooAbsData ;
class RooFitResult ;
class RooPlot ;
class RooRealVar ;
class RooWorkspace ;
class RooAbsStudy ;
#include "RooStudyPackage.h"
#include <list>
#include <string>

class RooStudyManager : public TNamed {
public:

  RooStudyManager(RooWorkspace& w) ;
  RooStudyManager(RooWorkspace& w, RooAbsStudy& study) ;
  RooStudyManager(const char* studyPackFileName) ;
  void addStudy(RooAbsStudy& study) ;

  // Interactive running
  void run(Int_t nExperiments) ;

  // PROOF-based parallel running
  void runProof(Int_t nExperiments, const char* proofHost="", Bool_t showGui=kTRUE) ;
  static void closeProof(Option_t *option = "s") ;

  // Batch running
  void prepareBatchInput(const char* studyName, Int_t nExpPerJob, Bool_t unifiedInput) ;
  void processBatchOutput(const char* filePat) ;

  RooWorkspace& wspace() { return _pkg->wspace() ; }
  std::list<RooAbsStudy*>& studies() { return _pkg->studies() ; }

protected:

  void aggregateData(TList* olist) ;
  void expandWildCardSpec(const char* spec, std::list<std::string>& result) ;

  RooStudyPackage* _pkg ;

  RooStudyManager(const RooStudyManager&) ;

  ClassDef(RooStudyManager,1) // A general purpose workspace oriented parallelizing study manager
} ;


#endif

