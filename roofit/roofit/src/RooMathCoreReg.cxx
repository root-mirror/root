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

/** \class RooMathCoreReg
    \ingroup Roofit

**/

#include "Riostream.h"
#include "RooMathCoreReg.h"
#include "RooCFunction1Binding.h"
#include "RooCFunction2Binding.h"
#include "RooCFunction3Binding.h"
#include "RooCFunction4Binding.h"
#include "Math/SpecFuncMathCore.h"
#include "Math/DistFuncMathCore.h"

namespace {

RooMathCoreReg dummy ;

}

RooMathCoreReg::RooMathCoreReg()
{
  // Import MathCore 'special' functions from ROOT::Math namespace
  RooCFunction1Ref<double,double>::fmap().add("ROOT::Math::erf",ROOT::Math::erf,"x") ;
  RooCFunction1Ref<double,double>::fmap().add("ROOT::Math::erfc",ROOT::Math::erfc,"x") ;
  RooCFunction1Ref<double,double>::fmap().add("ROOT::Math::tgamma",ROOT::Math::tgamma,"x") ;
  RooCFunction1Ref<double,double>::fmap().add("ROOT::Math::lgamma",ROOT::Math::lgamma,"x") ;
  RooCFunction2Ref<double,double,double>::fmap().add("ROOT::Math::inc_gamma",ROOT::Math::inc_gamma,"a","x") ;
  RooCFunction2Ref<double,double,double>::fmap().add("ROOT::Math::inc_gamma_c",ROOT::Math::inc_gamma_c,"a","x") ;
  RooCFunction2Ref<double,double,double>::fmap().add("ROOT::Math::beta",ROOT::Math::beta,"x","y") ;
  RooCFunction3Ref<double,double,double,double>::fmap().add("ROOT::Math::inc_beta",ROOT::Math::inc_beta,"x","a","b") ;

  // MathCore pdf functions from ROOT::Math namespace
  RooCFunction2Ref<double,unsigned int,double>::fmap().add("ROOT::Math::poisson_pdf",ROOT::Math::poisson_pdf, "n","mu") ;
  RooCFunction3Ref<double,double,double,double>::fmap().add("ROOT::Math::beta_pdf",ROOT::Math::beta_pdf,"x","a","b") ;
  RooCFunction3Ref<double,double,double,double>::fmap().add("ROOT::Math::breitwigner_pdf",ROOT::Math::breitwigner_pdf,"x","gamma","x0") ;
  RooCFunction3Ref<double,double,double,double>::fmap().add("ROOT::Math::cauchy_pdf",ROOT::Math::cauchy_pdf,"x","b","x0") ;
  RooCFunction3Ref<double,double,double,double>::fmap().add("ROOT::Math::chisquared_pdf",ROOT::Math::chisquared_pdf,"x","r","x0") ;
  RooCFunction3Ref<double,double,double,double>::fmap().add("ROOT::Math::exponential_pdf",ROOT::Math::exponential_pdf,"x","lambda","x0") ;
  RooCFunction3Ref<double,double,double,double>::fmap().add("ROOT::Math::gaussian_pdf",ROOT::Math::gaussian_pdf,"x","sigma","x0") ;
  RooCFunction3Ref<double,double,double,double>::fmap().add("ROOT::Math::landau_pdf",ROOT::Math::landau_pdf,"x","sigma","x0.") ;
  RooCFunction3Ref<double,double,double,double>::fmap().add("ROOT::Math::normal_pdf",ROOT::Math::normal_pdf,"x","sigma","x0") ;
  RooCFunction3Ref<double,double,double,double>::fmap().add("ROOT::Math::tdistribution_pdf",ROOT::Math::tdistribution_pdf,"x","r","x0") ;
  RooCFunction3Ref<double,unsigned int,double,unsigned int>::fmap().add("ROOT::Math::binomial_pdf",ROOT::Math::binomial_pdf,"int","double","unsigned") ;
  RooCFunction4Ref<double,double,double,double,double>::fmap().add("ROOT::Math::fdistribution_pdf",ROOT::Math::fdistribution_pdf,"x","n","m","x0") ;
  RooCFunction4Ref<double,double,double,double,double>::fmap().add("ROOT::Math::gamma_pdf",ROOT::Math::gamma_pdf,"x","alpha","theta","x0") ;
  RooCFunction4Ref<double,double,double,double,double>::fmap().add("ROOT::Math::lognormal_pdf",ROOT::Math::lognormal_pdf,"x","m","s","x0") ;
  RooCFunction4Ref<double,double,double,double,double>::fmap().add("ROOT::Math::uniform_pdf",ROOT::Math::uniform_pdf,"x","a","b","x0") ;

  // MathCore cdf functions from ROOT::Math namespace uint
  RooCFunction2Ref<double,unsigned int,double>::fmap().add("ROOT::Math::poisson_cdf_c",ROOT::Math::poisson_cdf_c, "n","mu") ;
  RooCFunction2Ref<double,unsigned int,double>::fmap().add("ROOT::Math::poisson_cdf",ROOT::Math::poisson_cdf, "n","mu") ;
  RooCFunction3Ref<double,double,double,double>::fmap().add("ROOT::Math::beta_cdf_c",ROOT::Math::beta_cdf_c,"x","a","b") ;
  RooCFunction3Ref<double,double,double,double>::fmap().add("ROOT::Math::beta_cdf",ROOT::Math::beta_cdf,"x","a","b") ;
  RooCFunction3Ref<double,double,double,double>::fmap().add("ROOT::Math::breitwigner_cdf_c",ROOT::Math::breitwigner_cdf_c,"x","gamma","x0") ;
  RooCFunction3Ref<double,double,double,double>::fmap().add("ROOT::Math::breitwigner_cdf",ROOT::Math::breitwigner_cdf,"x","gamma","x0") ;
  RooCFunction3Ref<double,double,double,double>::fmap().add("ROOT::Math::cauchy_cdf_c",ROOT::Math::cauchy_cdf_c,"x","b","x0") ;
  RooCFunction3Ref<double,double,double,double>::fmap().add("ROOT::Math::cauchy_cdf",ROOT::Math::cauchy_cdf,"x","b","x0") ;
  RooCFunction3Ref<double,double,double,double>::fmap().add("ROOT::Math::chisquared_cdf_c",ROOT::Math::chisquared_cdf_c,"x","r","x0") ;
  RooCFunction3Ref<double,double,double,double>::fmap().add("ROOT::Math::chisquared_cdf",ROOT::Math::chisquared_cdf,"x","r","x0") ;
  RooCFunction3Ref<double,double,double,double>::fmap().add("ROOT::Math::exponential_cdf_c",ROOT::Math::exponential_cdf_c,"x","lambda","x0") ;
  RooCFunction3Ref<double,double,double,double>::fmap().add("ROOT::Math::exponential_cdf",ROOT::Math::exponential_cdf,"x","lambda","x0") ;
  RooCFunction3Ref<double,double,double,double>::fmap().add("ROOT::Math::landau_cdf",ROOT::Math::landau_cdf,"x","sigma","x0") ;
  RooCFunction3Ref<double,double,double,double>::fmap().add("ROOT::Math::normal_cdf_c",ROOT::Math::normal_cdf_c,"x","sigma","x0") ;
  RooCFunction3Ref<double,double,double,double>::fmap().add("ROOT::Math::normal_cdf",ROOT::Math::normal_cdf,"x","sigma","x0") ;
  RooCFunction3Ref<double,double,double,double>::fmap().add("ROOT::Math::tdistribution_cdf_c",ROOT::Math::tdistribution_cdf_c,"x","r","x0") ;
  RooCFunction3Ref<double,double,double,double>::fmap().add("ROOT::Math::tdistribution_cdf",ROOT::Math::tdistribution_cdf,"x","r","x0") ;
  RooCFunction3Ref<double,unsigned int,double,unsigned int>::fmap().add("ROOT::Math::binomial_cdf_c",ROOT::Math::binomial_cdf_c,"int","double","unsigned") ;
  RooCFunction3Ref<double,unsigned int,double,unsigned int>::fmap().add("ROOT::Math::binomial_cdf",ROOT::Math::binomial_cdf,"int","double","unsigned") ;
  RooCFunction4Ref<double,double,double,double,double>::fmap().add("ROOT::Math::fdistribution_cdf_c",ROOT::Math::fdistribution_cdf_c,"x","n","m","x0") ;
  RooCFunction4Ref<double,double,double,double,double>::fmap().add("ROOT::Math::fdistribution_cdf",ROOT::Math::fdistribution_cdf,"x","n","m","x0") ;
  RooCFunction4Ref<double,double,double,double,double>::fmap().add("ROOT::Math::gamma_cdf_c",ROOT::Math::gamma_cdf_c,"x","alpha","theta","x0") ;
  RooCFunction4Ref<double,double,double,double,double>::fmap().add("ROOT::Math::gamma_cdf",ROOT::Math::gamma_cdf,"x","alpha","theta","x0") ;
  RooCFunction4Ref<double,double,double,double,double>::fmap().add("ROOT::Math::lognormal_cdf_c",ROOT::Math::lognormal_cdf_c,"x","m","s","x0") ;
  RooCFunction4Ref<double,double,double,double,double>::fmap().add("ROOT::Math::lognormal_cdf",ROOT::Math::lognormal_cdf,"x","m","s","x0") ;
  RooCFunction4Ref<double,double,double,double,double>::fmap().add("ROOT::Math::uniform_cdf_c",ROOT::Math::uniform_cdf_c,"x","a","b","x0") ;
  RooCFunction4Ref<double,double,double,double,double>::fmap().add("ROOT::Math::uniform_cdf",ROOT::Math::uniform_cdf,"x","a","b","x0") ;

  // MathCore quantile functions from ROOT::Math namespace
  RooCFunction2Ref<double,double,double>::fmap().add("ROOT::Math::cauchy_quantile_c",ROOT::Math::cauchy_quantile_c, "z", "b") ;
  RooCFunction2Ref<double,double,double>::fmap().add("ROOT::Math::cauchy_quantile",ROOT::Math::cauchy_quantile, "z", "b") ;
  RooCFunction2Ref<double,double,double>::fmap().add("ROOT::Math::breitwigner_quantile_c",ROOT::Math::breitwigner_quantile_c, "z", "gamma") ;
  RooCFunction2Ref<double,double,double>::fmap().add("ROOT::Math::breitwigner_quantile",ROOT::Math::breitwigner_quantile, "z", "gamma") ;
  RooCFunction2Ref<double,double,double>::fmap().add("ROOT::Math::chisquared_quantile_c",ROOT::Math::chisquared_quantile_c, "z", "r") ;
  RooCFunction2Ref<double,double,double>::fmap().add("ROOT::Math::exponential_quantile_c",ROOT::Math::exponential_quantile_c, "z", "lambda") ;
  RooCFunction2Ref<double,double,double>::fmap().add("ROOT::Math::exponential_quantile",ROOT::Math::exponential_quantile, "z", "lambda") ;
  RooCFunction2Ref<double,double,double>::fmap().add("ROOT::Math::gaussian_quantile_c",ROOT::Math::gaussian_quantile_c, "z", "sigma") ;
  RooCFunction2Ref<double,double,double>::fmap().add("ROOT::Math::gaussian_quantile",ROOT::Math::gaussian_quantile, "z", "sigma") ;
  RooCFunction2Ref<double,double,double>::fmap().add("ROOT::Math::normal_quantile_c",ROOT::Math::normal_quantile_c, "z", "sigma") ;
  RooCFunction2Ref<double,double,double>::fmap().add("ROOT::Math::normal_quantile",ROOT::Math::normal_quantile, "z", "sigma") ;
  //RooCFunction2Ref<double,double,double>::fmap().add("ROOT::Math::chisquared_quantile",ROOT::Math::chisquared_quantile, "z", "r") ;
  RooCFunction3Ref<double,double,double,double>::fmap().add("ROOT::Math::beta_quantile",ROOT::Math::beta_quantile,"x","a","b") ;
  RooCFunction3Ref<double,double,double,double>::fmap().add("ROOT::Math::beta_quantile_c",ROOT::Math::beta_quantile_c,"x","a","b") ;
  RooCFunction3Ref<double,double,double,double>::fmap().add("ROOT::Math::fdistribution_quantile",ROOT::Math::fdistribution_quantile,"z","n","m") ;
  RooCFunction3Ref<double,double,double,double>::fmap().add("ROOT::Math::fdistribution_quantile_c",ROOT::Math::fdistribution_quantile_c,"z","n","m") ;
  RooCFunction3Ref<double,double,double,double>::fmap().add("ROOT::Math::gamma_quantile_c",ROOT::Math::gamma_quantile_c,"z","alpha","theta") ;
  RooCFunction3Ref<double,double,double,double>::fmap().add("ROOT::Math::lognormal_quantile_c",ROOT::Math::lognormal_quantile_c,"x","m","s") ;
  RooCFunction3Ref<double,double,double,double>::fmap().add("ROOT::Math::lognormal_quantile",ROOT::Math::lognormal_quantile,"x","m","s") ;
  RooCFunction3Ref<double,double,double,double>::fmap().add("ROOT::Math::uniform_quantile_c",ROOT::Math::uniform_quantile_c,"z","a","b") ;
  RooCFunction3Ref<double,double,double,double>::fmap().add("ROOT::Math::uniform_quantile",ROOT::Math::uniform_quantile,"z","a","b") ;
  //RooCFunction3Ref<double,double,double,double>::fmap().add("ROOT::Math::gamma_quantile",ROOT::Math::gamma_quantile,"z","alpha","theta") ;
  RooCFunction3Ref<double,double,double,double>::fmap().add("ROOT::Math::gamma_quantile_c",ROOT::Math::gamma_quantile_c,"z","alpha","theta") ;

}
