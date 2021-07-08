 /*****************************************************************************
  * Project: RooFit                                                           *
  * @(#)root/roofit:$Id$ *
  *                                                                           *
  * RooFit Lognormal PDF                                                      *
  *                                                                           *
  * Author: Gregory Schott and Stefan Schmitz                                 *
  *                                                                           *
  *****************************************************************************/

#ifndef ROO_LOGNORMAL
#define ROO_LOGNORMAL

#include "RooAbsPdf.h"
#include "RooRealProxy.h"

class RooRealVar;

class RooLognormal : public RooAbsPdf {
public:
  RooLognormal() {} ;
  RooLognormal(const char *name, const char *title,
         RooAbsReal& _x, RooAbsReal& _m0, RooAbsReal& _k);
  RooLognormal(const RooLognormal& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooLognormal(*this,newname); }
  inline virtual ~RooLognormal() { }

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;

  Int_t getGenerator(const RooArgSet& directVars, RooArgSet &generateVars, Bool_t staticInitOK=kTRUE) const;
  void generateEvent(Int_t code);

protected:

  RooRealProxy x ;
  RooRealProxy m0 ;
  RooRealProxy k ;

  Double_t evaluate() const ;
  void computeBatch(double* output, size_t nEvents, rbc::DataMap& dataMap) const;
  inline bool canComputeBatchWithCuda() const { return true; }
  
private:

  ClassDef(RooLognormal,1) // log-normal PDF
};

#endif
