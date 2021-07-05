/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 *    File: $Id: RooBukinPdf.h,v 1.5 2007/07/12 20:30:49 wouter Exp $
 * Authors:                                                                  *
 *   RW, Ruddick William  UC Colorado        wor@slac.stanford.edu           *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/


// -- CLASS DESCRIPTION [PDF] --
// RooBukinPdf implements the NovosibirskA function

#ifndef ROO_BUKINPDF
#define ROO_BUKINPDF

#include "RooAbsPdf.h"
#include "RooRealProxy.h"

class RooRealVar;
class RooAbsReal;

class RooBukinPdf : public RooAbsPdf {
public:

  RooBukinPdf() {} ;
  RooBukinPdf(const char *name, const char *title,
         RooAbsReal& _x, RooAbsReal& _Xp,
         RooAbsReal& _sigp, RooAbsReal& _xi,
              RooAbsReal& _rho1, RooAbsReal& _rho2);

  RooBukinPdf(const RooBukinPdf& other,const char* name=0) ;

  virtual TObject* clone(const char* newname) const { return new RooBukinPdf(*this,newname);   }
  inline virtual ~RooBukinPdf() { }

protected:
  RooRealProxy x;
  RooRealProxy Xp;
  RooRealProxy sigp;
  RooRealProxy xi;
  RooRealProxy rho1;
  RooRealProxy rho2;
  Double_t evaluate() const;
  void computeBatch(double* output, size_t nEvents, rbc::DataMap& dataMap) const;

private:

  ClassDef(RooBukinPdf,2) // Variation of Novosibirsk PDF
};

#endif
