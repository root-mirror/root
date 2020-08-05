/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 *    File: $Id: RooPolynomial.h,v 1.8 2007/05/11 09:13:07 verkerke Exp $
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
#ifndef ROO_EXPPOLY
#define ROO_EXPPOLY

#include <vector>

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"

class RooRealVar;
class RooArgList ;

class RooExpPoly : public RooAbsPdf {
public:

  RooExpPoly() ;
  RooExpPoly(const char* name, const char* title, RooAbsReal& x) ;
  RooExpPoly(const char *name, const char *title,
      RooAbsReal& _x, const RooArgList& _coefList, Int_t lowestOrder=1) ;

  RooExpPoly(const RooExpPoly& other, const char* name = 0);
  virtual TObject* clone(const char* newname) const { return new RooExpPoly(*this, newname); }
  virtual ~RooExpPoly() ;

  Double_t getLogVal(const RooArgSet* nset) const override;

  std::string getFormulaExpression(bool expand) const;
  
  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const override;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const override; 
  
  void adjustLimits();
  
protected:

  RooRealProxy _x;
  RooListProxy _coefList ;
  Int_t _lowestOrder ;

  /// Evaluation
  Double_t evaluate() const override;
  Double_t evaluateLog() const;  

  ClassDefOverride(RooExpPoly,1) // ExpPoly PDF
};

#endif
