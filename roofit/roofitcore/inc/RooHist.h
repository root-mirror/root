/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooHist.h,v 1.22 2007/05/11 09:11:30 verkerke Exp $
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
#ifndef ROO_HIST
#define ROO_HIST

#include "TGraphAsymmErrors.h"
#include "RooPlotable.h"
#include "RooAbsData.h"
#include "RooAbsRealLValue.h"

class TH1;
class RooCurve ;

class RooHist : public TGraphAsymmErrors, public RooPlotable {
public:
  RooHist() ;
  RooHist(Double_t nominalBinWidth, Double_t nSigma= 1, Double_t xErrorFrac=1.0, Double_t scaleFactor=1.0);
  RooHist(const TH1 &data, Double_t nominalBinWidth= 0, Double_t nSigma= 1, RooAbsData::ErrorType=RooAbsData::Poisson,
     Double_t xErrorFrac=1.0, Bool_t correctForBinWidth=kTRUE, Double_t scaleFactor=1.);
  RooHist(const TH1 &data1, const TH1 &data2, Double_t nominalBinWidth= 0, Double_t nSigma= 1, RooAbsData::ErrorType=RooAbsData::Poisson,
     Double_t xErrorFrac=1.0, Bool_t efficiency=kFALSE, Double_t scaleFactor=1.0);
  RooHist(const RooHist& hist1, const RooHist& hist2, Double_t wgt1=1.0, Double_t wgt2=1.0,
     RooAbsData::ErrorType etype=RooAbsData::Poisson, Double_t xErrorFrac=1.0) ;
  RooHist(const RooAbsReal &f, RooAbsRealLValue &x, Double_t xErrorFrac=1.0, Double_t scaleFactor=1.0, const RooArgSet *normVars = 0, const RooFitResult* fr = 0);
  virtual ~RooHist();

  // add a datapoint for a bin with n entries, using a Poisson error
  void addBin(Axis_t binCenter, Double_t n, Double_t binWidth= 0, Double_t xErrorFrac=1.0, Double_t scaleFactor=1.0);
  // add a datapoint for a bin with n entries, using a given error
  void addBinWithError(Axis_t binCenter, Double_t n, Double_t elow, Double_t ehigh, Double_t binWidth= 0,
             Double_t xErrorFrac=1.0, Bool_t correctForBinWidth=kTRUE, Double_t scaleFactor=1.0);
  // add a datapoint for a bin with n entries, using a given x and y error
  void addBinWithXYError(Axis_t binCenter, Double_t n, Double_t exlow, Double_t exhigh, Double_t eylow, Double_t eyhigh,
                         Double_t scaleFactor=1.0);
  // add a datapoint for the asymmetry (n1-n2)/(n1+n2), using a binomial error
  void addAsymmetryBin(Axis_t binCenter, Int_t n1, Int_t n2, Double_t binWidth= 0, Double_t xErrorFrac=1.0, Double_t scaleFactor=1.0);
  // add a datapoint for the asymmetry (n1-n2)/(n1+n2), using sum-of-weights error
  void addAsymmetryBinWithError(Axis_t binCenter, Double_t n1, Double_t n2, Double_t en1, Double_t en2, Double_t binWidth= 0, Double_t xErrorFrac=1.0, Double_t scaleFactor=1.0);

  // add a datapoint for the efficiency (n1)/(n1+n2), using a binomial error
  void addEfficiencyBin(Axis_t binCenter, Int_t n1, Int_t n2, Double_t binWidth= 0, Double_t xErrorFrac=1.0, Double_t scaleFactor=1.0);
  // add a datapoint for the efficiency (n1)/(n1+n2), using a sum-of-weights error
  void addEfficiencyBinWithError(Axis_t binCenter, Double_t n1, Double_t n2, Double_t en1, Double_t en2, Double_t binWidth= 0, Double_t xErrorFrac=1.0, Double_t scaleFactor=1.0);

  virtual void printName(std::ostream& os) const ;
  virtual void printTitle(std::ostream& os) const ;
  virtual void printClassName(std::ostream& os) const ;
  virtual void printMultiline(std::ostream& os, Int_t content, Bool_t verbose=kFALSE, TString indent= "") const;

  inline virtual void Print(Option_t *options= 0) const {
    // Printing interface
    printStream(defaultPrintStream(),defaultPrintContents(options),defaultPrintStyle(options));
  }

  Double_t getFitRangeNEvt() const;
  Double_t getFitRangeNEvt(Double_t xlo, Double_t xhi) const ;
  Double_t getFitRangeBinW() const;
  inline Double_t getNominalBinWidth() const { return _nominalBinWidth; }
  inline void setRawEntries(Double_t n) { _rawEntries = n ; }

  Bool_t hasIdenticalBinning(const RooHist& other) const ;

  RooHist* makeResidHist(const RooCurve& curve,bool normalize=false, bool useAverage=false) const;
  RooHist* makePullHist(const RooCurve& curve, bool useAverage=false) const
    {return makeResidHist(curve,true,useAverage); }


  Bool_t isIdentical(const RooHist& other, Double_t tol=1e-6) const ;


protected:
  void initialize();
  Int_t roundBin(Double_t y);

private:
  Double_t _nominalBinWidth ; // Average bin width
  Double_t _nSigma ;          // Number of 'sigmas' error bars represent
  Double_t _entries ;         // Number of entries in histogram
  Double_t _rawEntries;        // Number of entries in source dataset

  ClassDef(RooHist,1) // 1-dimensional histogram with error bars
};

#endif
