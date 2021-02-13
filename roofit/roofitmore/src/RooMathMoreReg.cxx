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

/** \class RooMathMoreReg
    \ingroup Roofit

**/

#include "Riostream.h"
#include "RooMathMoreReg.h"
#include "RooCFunction1Binding.h"
#include "RooCFunction2Binding.h"
#include "RooCFunction3Binding.h"
#include "RooCFunction4Binding.h"
#include "Math/SpecFuncMathMore.h"
#include "Math/DistFuncMathMore.h"

namespace {

RooMathMoreReg dummy ;

}

RooMathMoreReg::RooMathMoreReg()
{
#ifdef MATHMORE

  // Import MathMore 'special' functions from ROOT::Math namespace
  RooCFunction1Ref<double,double>::fmap().add("ROOT::Math::comp_ellint_1",ROOT::Math::comp_ellint_1,"k") ;
  RooCFunction1Ref<double,double>::fmap().add("ROOT::Math::comp_ellint_2",ROOT::Math::comp_ellint_2,"k") ;
  RooCFunction1Ref<double,double>::fmap().add("ROOT::Math::expint",ROOT::Math::expint) ;
  RooCFunction1Ref<double,double>::fmap().add("ROOT::Math::riemann_zeta",ROOT::Math::riemann_zeta) ;
  RooCFunction2Ref<double,double,double>::fmap().add("ROOT::Math::cyl_bessel_i",ROOT::Math::cyl_bessel_i, "nu", "x") ;
  RooCFunction2Ref<double,double,double>::fmap().add("ROOT::Math::cyl_bessel_j",ROOT::Math::cyl_bessel_j, "nu", "x") ;
  RooCFunction2Ref<double,double,double>::fmap().add("ROOT::Math::cyl_bessel_k",ROOT::Math::cyl_bessel_k, "nu", "x") ;
  RooCFunction2Ref<double,double,double>::fmap().add("ROOT::Math::cyl_neumann",ROOT::Math::cyl_neumann, "nu", "x") ;
  RooCFunction2Ref<double,double,double>::fmap().add("ROOT::Math::ellint_1",ROOT::Math::ellint_1, "k", "phi") ;
  RooCFunction2Ref<double,double,double>::fmap().add("ROOT::Math::ellint_2",ROOT::Math::ellint_2, "k", "phi") ;
  RooCFunction2Ref<double,unsigned int,double>::fmap().add("ROOT::Math::laguerre",ROOT::Math::laguerre, "n", "x") ;
  RooCFunction2Ref<double,unsigned int,double>::fmap().add("ROOT::Math::legendre",ROOT::Math::legendre, "l", "x") ;
  RooCFunction2Ref<double,unsigned int,double>::fmap().add("ROOT::Math::sph_bessel",ROOT::Math::sph_bessel, "n", "x") ;
  RooCFunction2Ref<double,unsigned int,double>::fmap().add("ROOT::Math::sph_neumann",ROOT::Math::sph_neumann, "n", "x") ;
  RooCFunction3Ref<double,double,double,double>::fmap().add("ROOT::Math::conf_hyperg",ROOT::Math::conf_hyperg,"a","b","z") ;
  RooCFunction3Ref<double,double,double,double>::fmap().add("ROOT::Math::conf_hypergU",ROOT::Math::conf_hypergU,"a","b","z") ;
  RooCFunction3Ref<double,double,double,double>::fmap().add("ROOT::Math::ellint_3",ROOT::Math::ellint_3,"n","k","phi") ;
  RooCFunction3Ref<double,unsigned int,double,double>::fmap().add("ROOT::Math::assoc_laguerre",ROOT::Math::assoc_laguerre,"n","m","x") ;
  RooCFunction3Ref<double,unsigned int,unsigned int,double>::fmap().add("ROOT::Math::assoc_legendre",ROOT::Math::assoc_legendre,"l","m","x") ;
  RooCFunction3Ref<double,unsigned int,unsigned int,double>::fmap().add("ROOT::Math::sph_legendre",ROOT::Math::sph_legendre,"l","m","theta") ;

  // MathMore quantile functions from ROOT::Math namespace
  RooCFunction2Ref<double,double,double>::fmap().add("ROOT::Math::tdistribution_quantile_c",ROOT::Math::tdistribution_quantile_c,"z","r") ;
  RooCFunction2Ref<double,double,double>::fmap().add("ROOT::Math::tdistribution_quantile",ROOT::Math::tdistribution_quantile,"z","r") ;

#endif
}
