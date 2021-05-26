/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooHistError.h,v 1.14 2007/05/11 09:11:30 verkerke Exp $
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
#ifndef ROO_HIST_ERROR
#define ROO_HIST_ERROR

#include "Rtypes.h"
#include "RooNumber.h"
#include "RooAbsFunc.h"
#include <cmath>
#include <iostream>

class RooHistError {
public:
  static const RooHistError &instance();
  virtual ~RooHistError() {} ;

  Bool_t getPoissonInterval(Int_t n, Double_t &mu1, Double_t &mu2, Double_t nSigma= 1) const;
  Bool_t getBinomialIntervalAsym(Int_t n, Int_t m, Double_t &a1, Double_t &a2, Double_t nSigma= 1) const;
  Bool_t getBinomialIntervalEff(Int_t n, Int_t m, Double_t &a1, Double_t &a2, Double_t nSigma= 1) const;
  Bool_t getInterval(const RooAbsFunc *Qu, const RooAbsFunc *Ql, Double_t pointEstimate, Double_t stepSize,
           Double_t &lo, Double_t &hi, Double_t nSigma) const;

  static RooAbsFunc *createPoissonSum(Int_t n) ;
  static RooAbsFunc *createBinomialSum(Int_t n, Int_t m, Bool_t eff) ;

private:


  Bool_t getPoissonIntervalCalc(Int_t n, Double_t &mu1, Double_t &mu2, Double_t nSigma= 1) const;
  Double_t _poissonLoLUT[1000] ;
  Double_t _poissonHiLUT[1000] ;

  RooHistError();
  Double_t seek(const RooAbsFunc &f, Double_t startAt, Double_t step, Double_t value) const;

  // -----------------------------------------------------------
  // Define a 1-dim RooAbsFunc of mu that evaluates the sum:
  //
  //  Q(n|mu) = Sum_{k=0}^{n} P(k|mu)
  //
  // where P(n|mu) = exp(-mu) mu**n / n! is the Poisson PDF.
  // -----------------------------------------------------------
  class PoissonSum : public RooAbsFunc {
  public:
    inline PoissonSum(Int_t n) : RooAbsFunc(1), _n(n) { }
    inline Double_t operator()(const Double_t xvec[]) const {
      Double_t mu(xvec[0]),result(1),factorial(1);
      for(Int_t k= 1; k <= _n; k++) {
   factorial*= k;
   result+= pow(mu,k)/factorial;
      }
      return exp(-mu)*result;
    };
    inline Double_t getMinLimit(UInt_t /*index*/) const { return 0; }
    inline Double_t getMaxLimit(UInt_t /*index*/) const { return RooNumber::infinity() ; }
  private:
    Int_t _n;
  };

  // -----------------------------------------------------------
  // Define a 1-dim RooAbsFunc of a that evaluates the sum:
  //
  //  Q(n|n+m,a) = Sum_{k=0}^{n} B(k|n+m,a)
  //
  // where B(n|n+m,a) = (n+m)!/(n!m!) ((1+a)/2)**n ((1-a)/2)**m
  // is the Binomial PDF.
  // -----------------------------------------------------------
  class BinomialSumAsym : public RooAbsFunc {
  public:
    BinomialSumAsym(Int_t n, Int_t m) : RooAbsFunc(1), _n1(n), _N1(n+m) {
    }
    inline Double_t operator()(const Double_t xvec[]) const
      {
   Double_t p1(0.5*(1+xvec[0])),p2(1-p1),result(0),fact1(1),fact2(1);
   for(Int_t k= 0; k <= _n1; k++) {
     if(k > 0) { fact2*= k; fact1*= _N1-k+1; }
     result+= fact1/fact2*pow(p1,k)*pow(p2,_N1-k);
   }
   return result;
      };

    inline Double_t getMinLimit(UInt_t /*index*/) const { return -1; }
    inline Double_t getMaxLimit(UInt_t /*index*/) const { return +1; }

  private:
    Int_t _n1 ; // WVE Solaris CC5 doesn't want _n or _N here (likely compiler bug)
    Int_t _N1 ;
  } ;


  // -----------------------------------------------------------
  // Define a 1-dim RooAbsFunc of a that evaluates the sum:
  //
  //  Q(n|n+m,a) = Sum_{k=0}^{n} B(k|n+m,a)
  //
  // where B(n|n+m,a) = (n+m)!/(n!m!) ((1+a)/2)**n ((1-a)/2)**m
  // is the Binomial PDF.
  // -----------------------------------------------------------
  class BinomialSumEff : public RooAbsFunc {
  public:
    BinomialSumEff(Int_t n, Int_t m) : RooAbsFunc(1), _n1(n), _N1(n+m) {
    }
    inline Double_t operator()(const Double_t xvec[]) const
      {
   Double_t p1(xvec[0]),p2(1-p1),result(0),fact1(1),fact2(1);
   for(Int_t k= 0; k <= _n1; k++) {
     if(k > 0) { fact2*= k; fact1*= _N1-k+1; }
     result+= fact1/fact2*pow(p1,k)*pow(p2,_N1-k);
   }
   return result;
      };

    inline Double_t getMinLimit(UInt_t /*index*/) const { return  0; }
    inline Double_t getMaxLimit(UInt_t /*index*/) const { return +1; }

  private:
    Int_t _n1 ; // WVE Solaris CC5 doesn't want _n or _N here (likely compiler bug)
    Int_t _N1 ;
  } ;

  ClassDef(RooHistError,1) // Utility class for calculating histogram errors
};

#endif
