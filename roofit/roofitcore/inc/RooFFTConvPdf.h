/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

#ifndef ROOFFTCONVPDF
#define ROOFFTCONVPDF

#include "RooAbsCachedPdf.h"
#include "RooRealProxy.h"
#include "RooSetProxy.h"
#include "RooAbsReal.h"
#include "RooHistPdf.h"
#include "TVirtualFFT.h"

class RooRealVar;

///PDF for the numerical (FFT) convolution of two PDFs.
class RooFFTConvPdf : public RooAbsCachedPdf {
public:

  RooFFTConvPdf() {
    // coverity[UNINIT_CTOR]
  } ;
  RooFFTConvPdf(const char *name, const char *title, RooRealVar& convVar, RooAbsPdf& pdf1, RooAbsPdf& pdf2, Int_t ipOrder=2);
  RooFFTConvPdf(const char *name, const char *title, RooAbsReal& pdfConvVar, RooRealVar& convVar, RooAbsPdf& pdf1, RooAbsPdf& pdf2, Int_t ipOrder=2);
  RooFFTConvPdf(const RooFFTConvPdf& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooFFTConvPdf(*this,newname); }
  virtual ~RooFFTConvPdf() ;

  void setShift(Double_t val1, Double_t val2) { _shift1 = val1 ; _shift2 = val2 ; }
  void setCacheObservables(const RooArgSet& obs) { _cacheObs.removeAll() ; _cacheObs.add(obs) ; }
  const RooArgSet& cacheObservables() const { return _cacheObs ; }

  Double_t bufferFraction() const {
    // Return value of buffer fraction applied in FFT calculation array beyond either
    // end of the observable domain to reduce cyclical effects
    return _bufFrac ;
  }

  enum BufStrat { Extend=0, Mirror=1, Flat=2 } ;
  BufStrat bufferStrategy() const {
    // Return the strategy currently used to fill the buffer:
    // 'Extend' means is that the input p.d.f convolution observable range is widened to include the buffer range
    // 'Flat' means that the buffer is filled with the p.d.f. value at the boundary of the observable range
    // 'Mirror' means that the buffer is filled with a mirror image of the p.d.f. around the convolution observable boundary
    return _bufStrat ;
  }
  void setBufferStrategy(BufStrat bs) ;
  void setBufferFraction(Double_t frac) ;

  void printMetaArgs(std::ostream& os) const ;

  // Propagate maximum value estimate of pdf1 as convolution can only result in lower max values
  virtual Int_t getMaxVal(const RooArgSet& vars) const { return _pdf1.arg().getMaxVal(vars) ; }
  virtual Double_t maxVal(Int_t code) const { return _pdf1.arg().maxVal(code) ; }


protected:

  RooRealProxy _x ;       // Convolution observable
  RooRealProxy _xprime ;  // Input function representing value of convolution observable
  RooRealProxy _pdf1 ; // First input p.d.f
  RooRealProxy _pdf2 ; // Second input p.d.f
  RooSetProxy _params ; // Effective parameters of this p.d.f.

  void calcParams() ;
  Bool_t redirectServersHook(const RooAbsCollection& newServerList, Bool_t mustReplaceAll, Bool_t nameChange, Bool_t isRecursive) ;

  Double_t*  scanPdf(RooRealVar& obs, RooAbsPdf& pdf, const RooDataHist& hist, const RooArgSet& slicePos, Int_t& N, Int_t& N2, Int_t& zeroBin, Double_t shift) const ;

  class FFTCacheElem : public PdfCacheElem {
  public:
    FFTCacheElem(const RooFFTConvPdf& self, const RooArgSet* nset) ;
    ~FFTCacheElem() ;

    virtual RooArgList containedArgs(Action) ;

    TVirtualFFT* fftr2c1 ;
    TVirtualFFT* fftr2c2 ;
    TVirtualFFT* fftc2r ;

    RooAbsPdf* pdf1Clone ;
    RooAbsPdf* pdf2Clone ;

    RooAbsBinning* histBinning ;
    RooAbsBinning* scanBinning ;

  };

  friend class FFTCacheElem ;

  virtual Double_t evaluate() const { RooArgSet dummy(_x.arg()) ; return getVal(&dummy) ; } ; // dummy
  virtual const char* inputBaseName() const ;
  virtual RooArgSet* actualObservables(const RooArgSet& nset) const ;
  virtual RooArgSet* actualParameters(const RooArgSet& nset) const ;
  virtual RooAbsArg& pdfObservable(RooAbsArg& histObservable) const ;
  virtual void fillCacheObject(PdfCacheElem& cache) const ;
  void fillCacheSlice(FFTCacheElem& cache, const RooArgSet& slicePosition) const ;

  virtual PdfCacheElem* createCache(const RooArgSet* nset) const ;
  virtual TString histNameSuffix() const ;

  // mutable std:: map<const RooHistPdf*,CacheAuxInfo*> _cacheAuxInfo ; //! Auxiliary Cache information (do not persist)
  Double_t _bufFrac ; // Sampling buffer size as fraction of domain size
  BufStrat _bufStrat ; // Strategy to fill the buffer

  Double_t  _shift1 ;
  Double_t  _shift2 ;

  virtual RooAbsGenContext* genContext(const RooArgSet &vars, const RooDataSet *prototype=0,
                                       const RooArgSet* auxProto=0, Bool_t verbose= kFALSE) const ;

  friend class RooConvGenContext ;
  RooSetProxy  _cacheObs ; // Non-convolution observables that are also cached

private:

  void prepareFFTBinning(RooRealVar& convVar) const;

  ClassDef(RooFFTConvPdf,1) // Convolution operator p.d.f based on numeric Fourier transforms
};

#endif
