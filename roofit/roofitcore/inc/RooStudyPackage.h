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
#ifndef ROO_STUDY_PACKAGE
#define ROO_STUDY_PACKAGE

#include "TNamed.h"

class RooAbsPdf;
class RooDataSet ;
class RooAbsData ;
class RooFitResult ;
class RooPlot ;
class RooRealVar ;
class RooWorkspace ;
class RooAbsStudy ;
#include <list>

class RooStudyPackage : public TNamed {
public:

  RooStudyPackage() ;
  RooStudyPackage(RooWorkspace& w) ;
  RooStudyPackage(const RooStudyPackage&) ;
  void addStudy(RooAbsStudy& study) ;
  TObject* Clone(const char* /*newname*/="") const { return new RooStudyPackage(*this) ; }

  RooWorkspace& wspace() { return *_ws ; }
  std::list<RooAbsStudy*>& studies() { return _studies ; }

  void driver(Int_t nExperiments) ;

  Int_t initRandom() ;
  void initialize() ;
  void runOne() ;
  void run(Int_t nExperiments) ;
  void finalize() ;

  void exportData(TList* olist, Int_t seqno) ;

  static void processFile(const char* infile, Int_t nexp) ;

protected:

  RooWorkspace* _ws ;
  std::list<RooAbsStudy*> _studies ;


  ClassDef(RooStudyPackage,1) // A general purpose workspace oriented parallelizing study manager
} ;


#endif

