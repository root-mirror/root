/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 *    File: $Id: RooDSCBShape.h$                                             *
 * Authors:                                                                  *
 *    T. Skwarnicki modify RooCBShape to Asymmetrical Double-Sided CB        *
 *    Michael Wilkinson add to RooFit source                                 *
 *****************************************************************************/
#ifndef ROO_DSCB_SHAPE
#define ROO_DSCB_SHAPE

#include "RooAbsPdf.h"
#include "RooRealProxy.h"

class RooRealVar;

class RooDSCBShape : public RooAbsPdf {
public:
  RooDSCBShape() {} ;
  RooDSCBShape(const char *name, const char *title, RooAbsReal& _m,
	     RooAbsReal& _m0, RooAbsReal& _sigma,
	     RooAbsReal& _alphaL, RooAbsReal& _nL,
	     RooAbsReal& _alphaR, RooAbsReal& _nR);

  RooDSCBShape(const RooDSCBShape& other, const char* name = 0);
  virtual TObject* clone(const char* newname) const { return new RooDSCBShape(*this,newname); }

  inline virtual ~RooDSCBShape() { }

  virtual Int_t getAnalyticalIntegral( RooArgSet& allVars,  RooArgSet& analVars, const char* rangeName=0 ) const;
  virtual Double_t analyticalIntegral( Int_t code, const char* rangeName=0 ) const;

  // Optimized accept/reject generator support
  virtual Int_t getMaxVal(const RooArgSet& vars) const ;
  virtual Double_t maxVal(Int_t code) const ;

protected:

  Double_t ApproxErf(Double_t arg) const ;

  RooRealProxy m;
  RooRealProxy m0;
  RooRealProxy sigma;
  RooRealProxy alphaL;
  RooRealProxy nL;
  RooRealProxy alphaR;
  RooRealProxy nR;

  Double_t evaluate() const;

private:

  ClassDef(RooDSCBShape,1) // Asymmetrical Double-Sided Crystal Ball lineshape PDF
};

#endif
