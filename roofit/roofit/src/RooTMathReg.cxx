/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id$
 * Authors:                                                                  *
 *   WV, Wouter Verkerke, NIKHEF, verkerke@nikhef.nl                         *
 *                                                                           *
 * Copyright (c) 2000-2008, NIKHEF, Regents of the University of California  *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 *****************************************************************************/

/** \class RooTMathReg
    \ingroup Roofit

**/

#include "Riostream.h"
#include "RooTMathReg.h"
#include "RooCFunction1Binding.h"
#include "RooCFunction2Binding.h"
#include "RooCFunction3Binding.h"
#include "RooCFunction4Binding.h"
#include "TMath.h"

namespace {

RooTMathReg dummy ;

}

RooTMathReg::RooTMathReg()
{

  // Import function from TMath namespace
  RooCFunction1Ref<double,double>::fmap().add("TMath::Abs",TMath::Abs,"d") ;
  RooCFunction1Ref<double,double>::fmap().add("TMath::ACos",TMath::ACos) ;
  RooCFunction1Ref<double,double>::fmap().add("TMath::ACosH",TMath::ACosH,"t") ;
  RooCFunction1Ref<double,double>::fmap().add("TMath::ASin",TMath::ASin) ;
  RooCFunction1Ref<double,double>::fmap().add("TMath::ASinH",TMath::ASinH,"t") ;
  RooCFunction1Ref<double,double>::fmap().add("TMath::ATan",TMath::ATan) ;
  RooCFunction1Ref<double,double>::fmap().add("TMath::ATanH",TMath::ATanH,"t") ;
  RooCFunction1Ref<double,double>::fmap().add("TMath::BesselI0",TMath::BesselI0) ;
  RooCFunction1Ref<double,double>::fmap().add("TMath::BesselI1",TMath::BesselI1) ;
  RooCFunction1Ref<double,double>::fmap().add("TMath::BesselJ0",TMath::BesselJ0) ;
  RooCFunction1Ref<double,double>::fmap().add("TMath::BesselJ1",TMath::BesselJ1) ;
  RooCFunction1Ref<double,double>::fmap().add("TMath::BesselK0",TMath::BesselK0) ;
  RooCFunction1Ref<double,double>::fmap().add("TMath::BesselK1",TMath::BesselK1) ;
  RooCFunction1Ref<double,double>::fmap().add("TMath::BesselY0",TMath::BesselY0) ;
  RooCFunction1Ref<double,double>::fmap().add("TMath::BesselY1",TMath::BesselY1) ;
  RooCFunction1Ref<double,double>::fmap().add("TMath::Cos",TMath::Cos) ;
  RooCFunction1Ref<double,double>::fmap().add("TMath::CosH",TMath::CosH) ;
  RooCFunction1Ref<double,double>::fmap().add("TMath::DiLog",TMath::DiLog) ;
  RooCFunction1Ref<double,double>::fmap().add("TMath::Erf",TMath::Erf) ;
  RooCFunction1Ref<double,double>::fmap().add("TMath::Erfc",TMath::Erfc) ;
  RooCFunction1Ref<double,double>::fmap().add("TMath::ErfcInverse",TMath::ErfcInverse) ;
  RooCFunction1Ref<double,double>::fmap().add("TMath::ErfInverse",TMath::ErfInverse) ;
  RooCFunction1Ref<double,double>::fmap().add("TMath::Exp",TMath::Exp) ;
  RooCFunction1Ref<double,double>::fmap().add("TMath::Freq",TMath::Freq) ;
  RooCFunction1Ref<double,double>::fmap().add("TMath::Gamma",TMath::Gamma,"z") ;
  RooCFunction1Ref<double,double>::fmap().add("TMath::KolmogorovProb",TMath::KolmogorovProb,"z") ;
  RooCFunction1Ref<double,double>::fmap().add("TMath::LandauI",TMath::LandauI) ;
  RooCFunction1Ref<double,double>::fmap().add("TMath::LnGamma",TMath::LnGamma,"z") ;
  RooCFunction1Ref<double,double>::fmap().add("TMath::Log",TMath::Log) ;
  RooCFunction1Ref<double,double>::fmap().add("TMath::Log10",TMath::Log10) ;
  RooCFunction1Ref<double,double>::fmap().add("TMath::Log2",TMath::Log2) ;
  RooCFunction1Ref<double,double>::fmap().add("TMath::NormQuantile",TMath::NormQuantile,"p") ;
  RooCFunction1Ref<double,double>::fmap().add("TMath::Sin",TMath::Sin) ;
  RooCFunction1Ref<double,double>::fmap().add("TMath::SinH",TMath::SinH) ;
  RooCFunction1Ref<double,double>::fmap().add("TMath::Sqrt",TMath::Sqrt) ;
  RooCFunction1Ref<double,double>::fmap().add("TMath::StruveH0",TMath::StruveH0) ;
  RooCFunction1Ref<double,double>::fmap().add("TMath::StruveH1",TMath::StruveH1) ;
  RooCFunction1Ref<double,double>::fmap().add("TMath::StruveL0",TMath::StruveL0) ;
  RooCFunction1Ref<double,double>::fmap().add("TMath::StruveL1",TMath::StruveL1) ;
  RooCFunction1Ref<double,double>::fmap().add("TMath::Tan",TMath::Tan) ;
  RooCFunction1Ref<double,double>::fmap().add("TMath::TanH",TMath::TanH) ;
  RooCFunction1Ref<double,int>::fmap().add("TMath::Factorial",TMath::Factorial,"i") ;

  RooCFunction2Ref<double,double,double>::fmap().add("TMath::ATan2",TMath::ATan2, "y", "x") ;
  RooCFunction2Ref<double,double,double>::fmap().add("TMath::Beta",TMath::Beta, "p", "q") ;
  RooCFunction2Ref<double,double,double>::fmap().add("TMath::ChisquareQuantile",TMath::ChisquareQuantile, "p", "ndf") ;
  RooCFunction2Ref<double,double,double>::fmap().add("TMath::Gamma",TMath::Gamma, "a", "x") ;
  RooCFunction2Ref<double,double,double>::fmap().add("TMath::Hypot",TMath::Hypot, "x", "y") ;
  RooCFunction2Ref<double,double,double>::fmap().add("TMath::Poisson",TMath::Poisson, "x", "par") ;
  RooCFunction2Ref<double,double,double>::fmap().add("TMath::PoissonI",TMath::PoissonI, "x", "par") ;
  RooCFunction2Ref<double,double,double>::fmap().add("TMath::Power",TMath::Power, "x", "y") ;
  RooCFunction2Ref<double,double,double>::fmap().add("TMath::Sign",TMath::Sign, "a", "b") ;
  RooCFunction2Ref<double,double,double>::fmap().add("TMath::Student",TMath::Student, "T", "ndf") ;
  RooCFunction2Ref<double,double,double>::fmap().add("TMath::StudentI",TMath::StudentI, "T", "ndf") ;
  RooCFunction2Ref<double,int,double>::fmap().add("TMath::BesselI",TMath::BesselI, "n","x") ;
  RooCFunction2Ref<double,int,double>::fmap().add("TMath::BesselK",TMath::BesselK, "n", "x") ;
  RooCFunction2Ref<double,double,int>::fmap().add("TMath::Prob",TMath::Prob,"chi2","ndf") ;
  RooCFunction2Ref<double,double,int>::fmap().add("TMath::Ldexp",TMath::Ldexp,"x","exp") ;
  RooCFunction2Ref<double,int,int>::fmap().add("TMath::Binomial",TMath::Binomial,"n","k") ;

  RooCFunction3Ref<double,double,double,double>::fmap().add("TMath::BetaCf",TMath::BetaCf,"x","a","b") ;
  RooCFunction3Ref<double,double,double,double>::fmap().add("TMath::BetaDist",TMath::BetaDist,"x","p","q") ;
  RooCFunction3Ref<double,double,double,double>::fmap().add("TMath::BetaDistI",TMath::BetaDistI,"x","p","q") ;
  RooCFunction3Ref<double,double,double,double>::fmap().add("TMath::BetaIncomplete",TMath::BetaIncomplete,"x","a","b") ;
  RooCFunction3Ref<double,double,double,double>::fmap().add("TMath::BreitWigner",TMath::BreitWigner,"x","mean","gamma") ;
  RooCFunction3Ref<double,double,double,double>::fmap().add("TMath::CauchyDist",TMath::CauchyDist,"x","t","s") ;
  RooCFunction3Ref<double,double,double,double>::fmap().add("TMath::FDist",TMath::FDist,"F","N","M") ;
  RooCFunction3Ref<double,double,double,double>::fmap().add("TMath::FDistI",TMath::FDistI,"F","N","M") ;
  RooCFunction3Ref<double,double,double,double>::fmap().add("TMath::LaplaceDist",TMath::LaplaceDist,"x","alpha","beta") ;
  RooCFunction3Ref<double,double,double,double>::fmap().add("TMath::LaplaceDistI",TMath::LaplaceDistI,"x","alpha","beta") ;
  RooCFunction3Ref<double,double,double,double>::fmap().add("TMath::Vavilov",TMath::Vavilov,"x","kappa","beta2") ;
  RooCFunction3Ref<double,double,double,double>::fmap().add("TMath::VavilovI",TMath::VavilovI,"x","kappa","beta2") ;
  RooCFunction3Ref<double,double,double,bool>::fmap().add("TMath::StudentQuantile",TMath::StudentQuantile,"p","ndf","lower_tail") ;
  RooCFunction3Ref<double,double,int,int>::fmap().add("TMath::BinomialI",TMath::BinomialI,"p","n","k") ;

  RooCFunction4Ref<double,double,double,double,double>::fmap().add("TMath::GammaDist",TMath::GammaDist,"x","gamma","mu","beta") ;
  RooCFunction4Ref<double,double,double,double,double>::fmap().add("TMath::LogNormal",TMath::LogNormal,"x","sigma","theta","m") ;
  RooCFunction4Ref<double,double,double,double,int>::fmap().add("TMath::Voigt",TMath::Voigt,"x","sigma","lg","R") ;
  RooCFunction4Ref<double,double,double,double,bool>::fmap().add("TMath::Gaus",TMath::Gaus,"x","mean","sigma","norm") ;
  RooCFunction4Ref<double,double,double,double,bool>::fmap().add("TMath::Landau",TMath::Landau,"x","mpv","sigma","norm") ;

}
