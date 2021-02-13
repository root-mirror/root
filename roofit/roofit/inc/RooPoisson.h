 /***************************************************************************** 
  * Project: RooFit                                                           * 
  *                                                                           * 
  * Simple Poisson PDF
  * author: Kyle Cranmer <cranmer@cern.ch>
  *                                                                           * 
  *****************************************************************************/ 

#ifndef ROOPOISSON
#define ROOPOISSON

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
#include "RooTrace.h"
 
class RooPoisson : public RooAbsPdf {
public:
  RooPoisson() { _noRounding = kFALSE ;   } ;
  RooPoisson(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _mean, Bool_t noRounding=kFALSE);
  RooPoisson(const RooPoisson& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const override { return new RooPoisson(*this,newname); }
  inline virtual ~RooPoisson() {  }

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const override;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const override;

  Int_t getGenerator(const RooArgSet& directVars, RooArgSet &generateVars, Bool_t staticInitOK=kTRUE) const override;
  void generateEvent(Int_t code) override;
  
  /// Switch off/on rounding of `x` to the nearest integer.
  void setNoRounding(bool flag = kTRUE) {_noRounding = flag;}
  /// Switch on or off protection against negative means.
  void protectNegativeMean(bool flag = kTRUE) {_protectNegative = flag;}

protected:

  RooRealProxy x ;
  RooRealProxy mean ;
  Bool_t  _noRounding ;
  Bool_t  _protectNegative{true};
  
  Double_t evaluate() const override;
  RooSpan<double> evaluateSpan(RooBatchCompute::RunContext& evalData, const RooArgSet* normSet) const override;

  ClassDefOverride(RooPoisson,3) // A Poisson PDF
};
 
#endif
