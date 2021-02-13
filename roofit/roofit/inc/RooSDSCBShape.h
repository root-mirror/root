/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 *    File: $Id: RooSDSCBShape.h$                                            *
 * Authors:                                                                  *
 *    T. Skwarnicki modify RooCBShape to Symmetrical Double-Sided CB         *
 *    Michael Wilkinson add to RooFit source                                 *
 *****************************************************************************/
#ifndef ROO_SDSCB_SHAPE
#define ROO_SDSCB_SHAPE

#include "RooAbsPdf.h"
#include "RooRealProxy.h"

class RooRealVar;

class RooSDSCBShape : public RooAbsPdf {
public:
  RooSDSCBShape() {} ;
  RooSDSCBShape(const char *name, const char *title, RooAbsReal& _m,
	     RooAbsReal& _m0, RooAbsReal& _sigma,
	     RooAbsReal& _alpha, RooAbsReal& _n);

  RooSDSCBShape(const RooSDSCBShape& other, const char* name = 0);
  virtual TObject* clone(const char* newname) const { return new RooSDSCBShape(*this,newname); }

  inline virtual ~RooSDSCBShape() { }

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
  RooRealProxy alpha;
  RooRealProxy n;

  Double_t evaluate() const;

private:

  ClassDef(RooSDSCBShape,1) // Symmetrical Double-Sided Crystal Ball lineshape PDF
};

#endif
