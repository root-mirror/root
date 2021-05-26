/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooCurve.h,v 1.24 2007/05/11 09:11:30 verkerke Exp $
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
#ifndef ROO_CURVE
#define ROO_CURVE

#include "TGraph.h"
#include "RooPlotable.h"
#include <list>
#include <vector>
#include "TMatrixDfwd.h"

class RooAbsReal;
class RooRealVar;
class RooAbsFunc;
class RooArgSet;
class RooAbsRealLValue ;
class RooHist ;

class RooCurve : public TGraph, public RooPlotable {
public:
  RooCurve();
  enum WingMode { NoWings=0 ,Straight=1, Extended=2 } ;
  RooCurve(const RooAbsReal &func, RooAbsRealLValue &x, Double_t xlo, Double_t xhi, Int_t xbins,
      Double_t scaleFactor= 1, const RooArgSet *normVars= 0, Double_t prec= 1e-3, Double_t resolution= 1e-3,
      Bool_t shiftToZero=kFALSE, WingMode wmode=Extended, Int_t nEvalError=-1, Int_t doEEVal=kFALSE, Double_t eeVal=0,
      Bool_t showProgress=kFALSE);
  RooCurve(const char *name, const char *title, const RooAbsFunc &func, Double_t xlo,
      Double_t xhi, UInt_t minPoints, Double_t prec= 1e-3, Double_t resolution= 1e-3,
      Bool_t shiftToZero=kFALSE, WingMode wmode=Extended, Int_t nEvalError=-1, Int_t doEEVal=kFALSE, Double_t eeVal=0);
  virtual ~RooCurve();

  RooCurve(const char* name, const char* title, const RooCurve& c1, const RooCurve& c2, Double_t scale1=1., Double_t scale2=1.) ;

  void addPoint(Double_t x, Double_t y);

  Double_t getFitRangeBinW() const;
  Double_t getFitRangeNEvt(Double_t xlo, Double_t xhi) const ;
  Double_t getFitRangeNEvt() const;


  virtual void printName(std::ostream& os) const ;
  virtual void printTitle(std::ostream& os) const ;
  virtual void printClassName(std::ostream& os) const ;
  virtual void printMultiline(std::ostream& os, Int_t contents, Bool_t verbose=kFALSE, TString indent="") const;

  inline virtual void Print(Option_t *options= 0) const {
    // Printing interface
    printStream(defaultPrintStream(),defaultPrintContents(options),defaultPrintStyle(options));
  }

  Double_t chiSquare(const RooHist& hist, int nFitParam) const ;
  Int_t findPoint(Double_t value, Double_t tolerance=1e-10) const ;
  Double_t average(Double_t lo, Double_t hi) const ;
  Double_t interpolate(Double_t x, Double_t tolerance=1e-10) const ;

  Bool_t isIdentical(const RooCurve& other, Double_t tol=1e-6) const ;

  RooCurve* makeErrorBand(const std::vector<RooCurve*>& variations, Double_t Z=1) const ;
  RooCurve* makeErrorBand(const std::vector<RooCurve*>& plusVar, const std::vector<RooCurve*>& minusVar, const TMatrixD& V, Double_t Z=1) const ;

protected:

  void calcBandInterval(const std::vector<RooCurve*>& variations,Int_t i,Double_t Z,Double_t& lo, Double_t& hi, Bool_t approxGauss) const ;
  void calcBandInterval(const std::vector<RooCurve*>& plusVar, const std::vector<RooCurve*>& minusVar, Int_t i, const TMatrixD& V,
         Double_t Z,Double_t& lo, Double_t& hi) const ;

  void initialize();
  void addPoints(const RooAbsFunc &func, Double_t xlo, Double_t xhi,
       Int_t minPoints, Double_t prec, Double_t resolution, WingMode wmode,
       Int_t numee=0, Bool_t doEEVal=kFALSE, Double_t eeVal=0.,std::list<Double_t>* samplingHint=0) ;
  void addRange(const RooAbsFunc& func, Double_t x1, Double_t x2, Double_t y1,
      Double_t y2, Double_t minDy, Double_t minDx,
      Int_t numee=0, Bool_t doEEVal=kFALSE, Double_t eeVal=0.)  ;


  void shiftCurveToZero(Double_t prevYMax) ;

  Bool_t _showProgress ; //! Show progress indication when adding points

  ClassDef(RooCurve,1) // 1-dimensional smooth curve for use in RooPlots
};

#endif
