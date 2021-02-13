/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 * @(#)root/roofit:$Id$
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

/** \class RooBCPEffDecay
    \ingroup Roofit

PDF describing decay time distribution of B meson including effects of standard model CP violation.
This function can be analytically convolved with any RooResolutionModel implementation.
*/


#include "RooFit.h"

#include "Riostream.h"
#include "RooRealVar.h"
#include "RooRandom.h"
#include "RooBCPEffDecay.h"
#include "RooRealIntegral.h"

using namespace std;

ClassImp(RooBCPEffDecay);

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

RooBCPEffDecay::RooBCPEffDecay(const char *name, const char *title,
                RooRealVar& t, RooAbsCategory& tag,
                RooAbsReal& tau, RooAbsReal& dm,
                RooAbsReal& avgMistag, RooAbsReal& CPeigenval,
                RooAbsReal& a, RooAbsReal& b,
                RooAbsReal& effRatio, RooAbsReal& delMistag,
                const RooResolutionModel& model, DecayType type) :
  RooAbsAnaConvPdf(name,title,model,t),
  _absLambda("absLambda","Absolute value of lambda",this,a),
  _argLambda("argLambda","Arg(Lambda)",this,b),
  _effRatio("effRatio","B0/B0bar efficiency ratio",this,effRatio),
  _CPeigenval("CPeigenval","CP eigen value",this,CPeigenval),
  _avgMistag("avgMistag","Average mistag rate",this,avgMistag),
  _delMistag("delMistag","Delta mistag rate",this,delMistag),
  _t("t","time",this,t),
  _tau("tau","decay time",this,tau),
  _dm("dm","mixing frequency",this,dm),
  _tag("tag","CP state",this,tag),
  _genB0Frac(0),
  _type(type)
{
  switch(type) {
  case SingleSided:
    _basisExp = declareBasis("exp(-@0/@1)",RooArgList(tau,dm)) ;
    _basisSin = declareBasis("exp(-@0/@1)*sin(@0*@2)",RooArgList(tau,dm)) ;
    _basisCos = declareBasis("exp(-@0/@1)*cos(@0*@2)",RooArgList(tau,dm)) ;
    break ;
  case Flipped:
    _basisExp = declareBasis("exp(@0)/@1)",RooArgList(tau,dm)) ;
    _basisSin = declareBasis("exp(@0/@1)*sin(@0*@2)",RooArgList(tau,dm)) ;
    _basisCos = declareBasis("exp(@0/@1)*cos(@0*@2)",RooArgList(tau,dm)) ;
    break ;
  case DoubleSided:
    _basisExp = declareBasis("exp(-abs(@0)/@1)",RooArgList(tau,dm)) ;
    _basisSin = declareBasis("exp(-abs(@0)/@1)*sin(@0*@2)",RooArgList(tau,dm)) ;
    _basisCos = declareBasis("exp(-abs(@0)/@1)*cos(@0*@2)",RooArgList(tau,dm)) ;
    break ;
  }
}

////////////////////////////////////////////////////////////////////////////////
/// Copy constructor.

RooBCPEffDecay::RooBCPEffDecay(const RooBCPEffDecay& other, const char* name) :
  RooAbsAnaConvPdf(other,name),
  _absLambda("absLambda",this,other._absLambda),
  _argLambda("argLambda",this,other._argLambda),
  _effRatio("effRatio",this,other._effRatio),
  _CPeigenval("CPeigenval",this,other._CPeigenval),
  _avgMistag("avgMistag",this,other._avgMistag),
  _delMistag("delMistag",this,other._delMistag),
  _t("t",this,other._t),
  _tau("tau",this,other._tau),
  _dm("dm",this,other._dm),
  _tag("tag",this,other._tag),
  _genB0Frac(other._genB0Frac),
  _type(other._type),
  _basisExp(other._basisExp),
  _basisSin(other._basisSin),
  _basisCos(other._basisCos)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

RooBCPEffDecay::~RooBCPEffDecay()
{
}

////////////////////////////////////////////////////////////////////////////////
/// B0    : _tag = +1
///
/// B0bar : _tag = -1
/// \param[in] basisIndex

Double_t RooBCPEffDecay::coefficient(Int_t basisIndex) const
{
  if (basisIndex==_basisExp) {
    //exp term: (1 -/+ dw)(1+a^2)/2
    return (1 - _tag*_delMistag)*(1+_absLambda*_absLambda)/2 ;
    // =    1 + _tag*deltaDil/2
  }

  if (basisIndex==_basisSin) {
    //sin term: +/- (1-2w)*ImLambda(= -etaCP * |l| * sin2b)
    return -1*_tag*(1-2*_avgMistag)*_CPeigenval*_absLambda*_argLambda ;
    // =   _tag*avgDil * ...
  }

  if (basisIndex==_basisCos) {
    //cos term: +/- (1-2w)*(1-a^2)/2
    return -1*_tag*(1-2*_avgMistag)*(1-_absLambda*_absLambda)/2 ;
    // =   -_tag*avgDil * ...
  }

  return 0 ;
}

////////////////////////////////////////////////////////////////////////////////

Int_t RooBCPEffDecay::getCoefAnalyticalIntegral(Int_t /*code*/, RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const
{
  if (rangeName) return 0 ;

  if (matchArgs(allVars,analVars,_tag)) return 1 ;
  return 0 ;
}

////////////////////////////////////////////////////////////////////////////////

Double_t RooBCPEffDecay::coefAnalyticalIntegral(Int_t basisIndex, Int_t code, const char* /*rangeName*/) const
{
  switch(code) {
    // No integration
  case 0: return coefficient(basisIndex) ;

    // Integration over 'tag'
  case 1:
    if (basisIndex==_basisExp) {
      return (1+_absLambda*_absLambda) ;
    }

    if (basisIndex==_basisSin || basisIndex==_basisCos) {
      return 0 ;
    }
    break ;

  default:
    assert(0) ;
  }

  return 0 ;
}

////////////////////////////////////////////////////////////////////////////////

Int_t RooBCPEffDecay::getGenerator(const RooArgSet& directVars, RooArgSet &generateVars, Bool_t staticInitOK) const
{
  if (staticInitOK) {
    if (matchArgs(directVars,generateVars,_t,_tag)) return 2 ;
  }
  if (matchArgs(directVars,generateVars,_t)) return 1 ;
  return 0 ;
}

////////////////////////////////////////////////////////////////////////////////

void RooBCPEffDecay::initGenerator(Int_t code)
{
  if (code==2) {
    // Calculate the fraction of mixed events to generate
    Double_t sumInt = RooRealIntegral("sumInt","sum integral",*this,RooArgSet(_t.arg(),_tag.arg())).getVal() ;
    _tag = 1 ;
    Double_t b0Int = RooRealIntegral("mixInt","mix integral",*this,RooArgSet(_t.arg())).getVal() ;
    _genB0Frac = b0Int/sumInt ;
  }
}

////////////////////////////////////////////////////////////////////////////////
/// Generates mix-state dependent.
/// \param[in] code

void RooBCPEffDecay::generateEvent(Int_t code)
{
  if (code==2) {
    Double_t rand = RooRandom::uniform() ;
    _tag = (rand<=_genB0Frac) ? 1 : -1 ;
  }

  // Generate delta-t dependent
  while(1) {
    Double_t rand = RooRandom::uniform() ;
    Double_t tval(0) ;

    switch(_type) {
    case SingleSided:
      tval = -_tau*log(rand);
      break ;
    case Flipped:
      tval= +_tau*log(rand);
      break ;
    case DoubleSided:
      tval = (rand<=0.5) ? -_tau*log(2*rand) : +_tau*log(2*(rand-0.5)) ;
      break ;
    }

    // Accept event if T is in generated range
    Double_t maxDil = 1.0 ;
    Double_t al2 = _absLambda*_absLambda ;
    Double_t maxAcceptProb = (1+al2) + fabs(maxDil*_CPeigenval*_absLambda*_argLambda) + fabs(maxDil*(1-al2)/2);
    Double_t acceptProb    = (1+al2)/2*(1-_tag*_delMistag)
                           - (_tag*(1-2*_avgMistag))*(_CPeigenval*_absLambda*_argLambda)*sin(_dm*tval)
                           - (_tag*(1-2*_avgMistag))*(1-al2)/2*cos(_dm*tval);

    Bool_t accept = maxAcceptProb*RooRandom::uniform() < acceptProb ? kTRUE : kFALSE ;

    if (tval<_t.max() && tval>_t.min() && accept) {
      _t = tval ;
      break ;
    }
  }

}
