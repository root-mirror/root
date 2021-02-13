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
#ifndef ROO_GEN_FIT_STUDY
#define ROO_GEN_FIT_STUDY

#include "RooAbsStudy.h"

class RooAbsPdf;
class RooDataSet ;
class RooAbsData ;
class RooFitResult ;
class RooPlot ;
class RooRealVar ;
class RooWorkspace ;
class RooAbsGenContext ;

#include "RooArgSet.h"
#include "RooLinkedList.h"
#include "RooAbsPdf.h"
#include <string>

class RooGenFitStudy : public RooAbsStudy {
public:

  RooGenFitStudy(const char* name=0, const char* title=0) ;
  RooGenFitStudy(const RooGenFitStudy& other) ;
  virtual ~RooGenFitStudy() ;
  virtual RooAbsStudy* clone(const char* newname="") const { return new RooGenFitStudy(newname?newname:GetName(),GetTitle()) ; }

  void setGenConfig(const char* pdfName, const char* obsName, const RooCmdArg& arg1=RooCmdArg(),const RooCmdArg& arg2=RooCmdArg(),const RooCmdArg& arg3=RooCmdArg()) ;
  void setFitConfig(const char* pdfName, const char* obsName, const RooCmdArg& arg1=RooCmdArg(),const RooCmdArg& arg2=RooCmdArg(),const RooCmdArg& arg3=RooCmdArg()) ;
 
  virtual Bool_t attach(RooWorkspace& w) ;
  virtual Bool_t initialize() ;
  virtual Bool_t execute() ;
  virtual Bool_t finalize() ;

  void Print(Option_t *options= 0) const;

 protected:


  std::string _genPdfName ;
  std::string _genObsName ;
  std::string _fitPdfName ;
  std::string _fitObsName ;
  RooLinkedList _genOpts ; 
  RooLinkedList _fitOpts ;

  RooAbsPdf* _genPdf ; //!
  RooArgSet _genObs ;  //!
  RooAbsPdf* _fitPdf ; //!
  RooArgSet _fitObs ; //!

  RooAbsPdf::GenSpec* _genSpec ; //!
  RooRealVar* _nllVar ; //!
  RooRealVar* _ngenVar ; //!
  RooArgSet* _params ; //!
  RooArgSet* _initParams; //!
  
  ClassDef(RooGenFitStudy,1) // Generate-and-Fit study module
} ;


#endif

