/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 *    File: $Id$
 * Authors:                                                                  *
 *   Kyle Cranmer
 *                                                                           *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/
#ifndef ROO_BERNSTEIN
#define ROO_BERNSTEIN

#include "RooAbsPdf.h"
#include "RooTemplateProxy.h"
#include "RooRealVar.h"
#include "RooListProxy.h"
#include "RooAbsRealLValue.h"
#include <string>

class RooRealVar;
class RooArgList;

class RooBernstein : public RooAbsPdf {
public:

  RooBernstein() ;
  RooBernstein(const char *name, const char *title,
               RooAbsRealLValue& _x, const RooArgList& _coefList) ;

  RooBernstein(const RooBernstein& other, const char* name = 0);
  virtual TObject* clone(const char* newname) const { return new RooBernstein(*this, newname); }
  inline virtual ~RooBernstein() { }

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;
  void selectNormalizationRange(const char* rangeName=0, Bool_t force=kFALSE) ;

private:
  
  RooTemplateProxy<RooAbsRealLValue> _x ;
  RooListProxy _coefList ;
  std::string _refRangeName ;

  Double_t evaluate() const;
  void computeBatch(double* output, size_t nEvents, rbc::DataMap& dataMap) const;

  ClassDef(RooBernstein,2) // Bernstein polynomial PDF
};

#endif
