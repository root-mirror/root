/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 * @(#)root/roofit:$Id$
 * Authors:                                                                  *
 *   UE, Ulrik Egede,     RAL,               U.Egede@rl.ac.uk                *
 *   MT, Max Turri,       UC Santa Cruz      turri@slac.stanford.edu         *
 *   CC, Chih-hsiang Cheng, Stanford         chcheng@slac.stanford.edu       *
 *   DK, David Kirkby,    UC Irvine,         dkirkby@uci.edu                 *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          RAL and Stanford University. All rights reserved.*
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

/** \class RooDstD0BG
    \ingroup Roofit

Special p.d.f shape that can be used to model the background of
D*-D0 mass difference distributions. It computes

\f[
  \mathrm{RooDSTD0}(m \, | \, m_0, A, B, C) =
    \left(1 - \exp\left(-\frac{m - m_0}{C}\right) \right)
    \cdot \left(\frac{m}{m_0}\right)^A + B
    \cdot \left(\frac{m}{m_0} - 1 \right)
\f]
**/

#include "RooDstD0BG.h"
#include "RooFit.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooIntegrator1D.h"
#include "RooAbsFunc.h"
#include "RooBatchCompute.h"

#include "TMath.h"

#include <cmath>
using namespace std;

ClassImp(RooDstD0BG);

////////////////////////////////////////////////////////////////////////////////

RooDstD0BG::RooDstD0BG(const char *name, const char *title,
             RooAbsReal& _dm, RooAbsReal& _dm0,
             RooAbsReal& _c, RooAbsReal& _a, RooAbsReal& _b) :
  RooAbsPdf(name,title),
  dm("dm","Dstar-D0 Mass Diff",this, _dm),
  dm0("dm0","Threshold",this, _dm0),
  C("C","Shape Parameter",this, _c),
  A("A","Shape Parameter 2",this, _a),
  B("B","Shape Parameter 3",this, _b)
{
}

////////////////////////////////////////////////////////////////////////////////

RooDstD0BG::RooDstD0BG(const RooDstD0BG& other, const char *name) :
  RooAbsPdf(other,name), dm("dm",this,other.dm), dm0("dm0",this,other.dm0),
  C("C",this,other.C), A("A",this,other.A), B("B",this,other.B)
{
}

////////////////////////////////////////////////////////////////////////////////

Double_t RooDstD0BG::evaluate() const
{
  Double_t arg= dm- dm0;
  if (arg <= 0 ) return 0;
  Double_t ratio= dm/dm0;
  Double_t val= (1- exp(-arg/C))* TMath::Power(ratio, A) + B*(ratio-1);

  return (val > 0 ? val : 0) ;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute multiple values of D*-D0 mass difference distribution.  
RooSpan<double> RooDstD0BG::evaluateSpan(RooBatchCompute::RunContext& evalData, const RooArgSet* normSet) const {
  return RooBatchCompute::dispatch->computeDstD0BG(this, evalData, dm->getValues(evalData, normSet), dm0->getValues(evalData, normSet), C->getValues(evalData, normSet), A->getValues(evalData, normSet), B->getValues(evalData, normSet));
}


////////////////////////////////////////////////////////////////////////////////
/// if (matchArgs(allVars,analVars,dm)) return 1 ;

Int_t RooDstD0BG::getAnalyticalIntegral(RooArgSet& /*allVars*/, RooArgSet& /*analVars*/, const char* /*rangeName*/) const
{
  return 0 ;
}

////////////////////////////////////////////////////////////////////////////////

Double_t RooDstD0BG::analyticalIntegral(Int_t code, const char* rangeName) const
{
  switch(code) {
  case 1:
    {
      Double_t min= dm.min(rangeName);
      Double_t max= dm.max(rangeName);
      if (max <= dm0 ) return 0;
      else if (min < dm0) min = dm0;

      Bool_t doNumerical= kFALSE;
      if ( A != 0 ) doNumerical= kTRUE;
      else if (B < 0) {
   // If b<0, pdf can be negative at large dm, the integral should
   // only up to where pdf hits zero. Better solution should be
   // solve the zero and use it as max.
   // Here we check this whether pdf(max) < 0. If true, let numerical
   // integral take care of. ( kind of ugly!)
   if ( 1- exp(-(max-dm0)/C) + B*(max/dm0 -1) < 0) doNumerical= kTRUE;
      }
      if ( ! doNumerical ) {
   return (max-min)+ C* exp(dm0/C)* (exp(-max/C)- exp(-min/C)) +
     B * (0.5* (max*max - min*min)/dm0 - (max- min));
      } else {
   // In principle the integral for a!=0  can be done analytically.
   // It involves incomplete Gamma function, TMath::Gamma(a+1,m/c),
   // which is not defined for a < -1. And the whole expression is
   // not stable for m/c >> 1.
   // Do numerical integral
   RooArgSet vset(dm.arg(),"vset");
   RooAbsFunc *func= bindVars(vset);
   RooIntegrator1D integrator(*func,min,max);
   return integrator.integral();
      }
    }
  }

  assert(0) ;
  return 0 ;
}
