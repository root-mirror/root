/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooRealAnalytic.h,v 1.8 2007/05/11 09:11:30 verkerke Exp $
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
#ifndef ROO_REAL_ANALYTIC
#define ROO_REAL_ANALYTIC

#include "RooRealBinding.h"

class RooRealAnalytic : public RooRealBinding {
public:
  inline RooRealAnalytic(const RooAbsReal &func, const RooArgSet &vars, Int_t code, const RooArgSet* normSet=0, const TNamed* rangeName=0) :
    RooRealBinding(func,vars,normSet,rangeName), _code(code) { }
  inline virtual ~RooRealAnalytic() { }

  virtual Double_t operator()(const Double_t xvector[]) const override;
  RooSpan<const double> getValues(std::vector<RooSpan<const double>> coordinates) const override;

protected:
  Int_t _code;

private:
  mutable std::unique_ptr<std::vector<double>> _batchBuffer; //! Buffer for handing out spans.

  ClassDefOverride(RooRealAnalytic,0) // Function binding to an analytical integral of a RooAbsReal
};

#endif

