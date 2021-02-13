/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 * @(#)root/roofitcore:$Id$
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

/**
\file RooConvCoefVar.cxx
\class RooConvCoefVar
\ingroup Roofitcore

RooConvCoefVar is an auxilary class that represents the coefficient
of a RooAbsAnaConvPdf implementation as a separate RooAbsReal object
to be able to interface these coefficient terms with the generic
RooRealIntegral integration mechanism
**/

#include "RooFit.h"

#include "RooAbsAnaConvPdf.h"
#include "RooConvCoefVar.h"

using namespace std;

ClassImp(RooConvCoefVar);
;


////////////////////////////////////////////////////////////////////////////////
/// Constuctor given a RooAbsAnaConvPdf a coefficient index and a set with the
/// convoluted observable(s)

RooConvCoefVar::RooConvCoefVar(const char *name, const char *title, const RooAbsAnaConvPdf& input, 
			       Int_t coefIdx, const RooArgSet* varList) :
  RooAbsReal(name,title),
  _varSet("varSet","Set of coefficient variables",this),
  _convPdf("convPdf","Convoluted PDF",this,(RooAbsReal&)input,kFALSE,kFALSE),
  _coefIdx(coefIdx)
{
  if (varList) _varSet.add(*varList) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Copy constructor

RooConvCoefVar::RooConvCoefVar(const RooConvCoefVar& other, const char* name) :
  RooAbsReal(other,name),
  _varSet("varSet",this,other._varSet),
  _convPdf("convPdf",this,other._convPdf),
  _coefIdx(other._coefIdx)
{
}



////////////////////////////////////////////////////////////////////////////////
/// Return value of chosen coefficient

Double_t RooConvCoefVar::getValV(const RooArgSet*) const 
{ 
  return evaluate() ; 
}



////////////////////////////////////////////////////////////////////////////////
/// Return value of chosen coefficient

Double_t RooConvCoefVar::evaluate() const 
{
  return ((RooAbsAnaConvPdf&)_convPdf.arg()).coefficient(_coefIdx) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Return analytical integration capabilities of chosen coefficient

Int_t RooConvCoefVar::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const 
{
  Int_t code = ((RooAbsAnaConvPdf&)_convPdf.arg()).getCoefAnalyticalIntegral(_coefIdx,allVars,analVars,rangeName) ;
  return code ;
}



////////////////////////////////////////////////////////////////////////////////
/// Return analytical integral of chosen coefficient

Double_t RooConvCoefVar::analyticalIntegral(Int_t code, const char* rangeName) const 
{
  return ((RooAbsAnaConvPdf&)_convPdf.arg()).coefAnalyticalIntegral(_coefIdx,code,rangeName) ;
}

