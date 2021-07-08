/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels
 *
 * authors: Stefan A. Schmitz, Gregory Schott
 * implementation of the Gamma distribution (class structure derived
 * from the class RooGaussian by Wouter Verkerke and David Kirkby)
 *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

#ifndef ROO_GAMMA
#define ROO_GAMMA

#include "RooAbsPdf.h"
#include "RooRealProxy.h"

class RooGamma : public RooAbsPdf {
public:
  RooGamma() {} ;
  RooGamma(const char *name, const char *title,
         RooAbsReal& _x, RooAbsReal& _gamma, RooAbsReal& _beta, RooAbsReal& _mu);
  RooGamma(const RooGamma& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooGamma(*this,newname); }
  inline virtual ~RooGamma() { }

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;

  Int_t getGenerator(const RooArgSet& directVars, RooArgSet &generateVars, Bool_t staticInitOK=kTRUE) const;
  void generateEvent(Int_t code);

protected:

  RooRealProxy x ;
  RooRealProxy gamma ;
  RooRealProxy beta ;
  RooRealProxy mu ;

  Double_t evaluate() const ;
  void computeBatch(double* output, size_t nEvents, rbc::DataMap& dataMap) const;
  inline bool canComputeBatchWithCuda() const { return true; }

private:

  ClassDef(RooGamma,1) // Gaussian PDF
};

#endif
