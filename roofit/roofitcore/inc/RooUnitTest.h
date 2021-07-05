/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id$
 * Authors:                                                                  *
 *   WV, Wouter Verkerke, NIKHEF, verkerke@nikhef.nl                         *
 *                                                                           *
 * Copyright (c) 2000-2011, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/
#ifndef ROO_UNIT_TEST
#define ROO_UNIT_TEST

#include "TNamed.h"
#include "RooTable.h"
#include "RooWorkspace.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "TFile.h"
#include "TH1.h"
#include <list>
#include <string>
#include <utility>

/*
 * The tolerance for the curve test is put to 0.4 instead of 0.2 to take into
 * account the small variations in the values of the likelihood which can occur
 * in presence of a different treatment of floating point numbers.
 */

class RooUnitTest : public TNamed {
public:
  RooUnitTest(const char* name, TFile* refFile, Bool_t writeRef, Int_t verbose, int batchMode=0) ;
  ~RooUnitTest() ;
  
  void setDebug(Bool_t flag) { _debug = flag ; }
  void setSilentMode() ;
  void clearSilentMode() ;
  void regPlot(RooPlot* frame, const char* refName) ;  
  void regResult(RooFitResult* r, const char* refName) ;
  void regValue(Double_t value, const char* refName) ;
  void regTable(RooTable* t, const char* refName) ;
  void regWS(RooWorkspace* ws, const char* refName) ;
  void regTH(TH1* h, const char* refName) ;
  RooWorkspace* getWS(const char* refName) ;
  Bool_t runTest() ;
  Bool_t runCompTests() ;
  Bool_t areTHidentical(TH1* htest, TH1* href) ;

  virtual Bool_t isTestAvailable() { return kTRUE ; }
  virtual Bool_t testCode() = 0 ;  

  virtual Double_t htol() { return 5e-4 ; } // histogram test tolerance (KS dist != prob)
#ifdef R__FAST_MATH
  virtual Double_t ctol() { return 2e-3 ; } // curve test tolerance
#else
  virtual Double_t ctol() { return 4e-3 ; } // curve test tolerance
#endif
  virtual Double_t fptol() { return 1e-5 ; } // fit parameter test tolerance
  virtual Double_t fctol() { return 1e-4 ; } // fit correlation test tolerance
  virtual Double_t vtol() { return 1e-3 ; } // value test tolerance

  static void setMemDir(TDirectory* memDir);

protected:

  static TDirectory* gMemDir ;

  TFile* _refFile ;
  Bool_t _debug ;
  Bool_t _write ;
  Int_t _verb ;
  int _batchMode=0;
   std::list<std::pair<RooPlot*, std::string> > _regPlots ;
   std::list<std::pair<RooFitResult*, std::string> > _regResults ;
   std::list<std::pair<Double_t, std::string> > _regValues ;
   std::list<std::pair<RooTable*,std::string> > _regTables ;
   std::list<std::pair<RooWorkspace*,std::string> > _regWS ;
   std::list<std::pair<TH1*,std::string> > _regTH ;

  ClassDef(RooUnitTest,1) ; // Abstract base class for RooFit/RooStats unit regression tests
} ;
#endif
