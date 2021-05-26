/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 * @(#)root/roofit:$Id$
 * Authors:                                                                  *
 *   DB, Dieter Best,     UC Irvine,         best@slac.stanford.edu          *
 *   HT, Hirohisa Tanaka  SLAC               tanaka@slac.stanford.edu        *
 *                                                                           *
 *   Updated version with analytical integral                                *
 *   MP, Marko Petric,   J. Stefan Institute,  marko.petric@ijs.si           *
 *                                                                           *
 * Copyright (c) 2000-2013, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

/** \class RooNovosibirsk
    \ingroup Roofit

RooNovosibirsk implements the Novosibirsk function

Function taken from H. Ikeda et al. NIM A441 (2000), p. 401 (Belle Collaboration)

**/
#include "RooNovosibirsk.h"
#include "RooFit.h"
#include "RooRealVar.h"
#include "RooBatchCompute.h"

#include "TMath.h"

#include <cmath>
using namespace std;

ClassImp(RooNovosibirsk);

////////////////////////////////////////////////////////////////////////////////

RooNovosibirsk::RooNovosibirsk(const char *name, const char *title,
              RooAbsReal& _x,     RooAbsReal& _peak,
              RooAbsReal& _width, RooAbsReal& _tail) :
  // The two addresses refer to our first dependent variable and
  // parameter, respectively, as declared in the rdl file
  RooAbsPdf(name, title),
  x("x","x",this,_x),
  width("width","width",this,_width),
  peak("peak","peak",this,_peak),
  tail("tail","tail",this,_tail)
{
}

////////////////////////////////////////////////////////////////////////////////

RooNovosibirsk::RooNovosibirsk(const RooNovosibirsk& other, const char *name):
  RooAbsPdf(other,name),
  x("x",this,other.x),
  width("width",this,other.width),
  peak("peak",this,other.peak),
  tail("tail",this,other.tail)
{
}

////////////////////////////////////////////////////////////////////////////////
///If tail=eta=0 the Belle distribution becomes gaussian

Double_t RooNovosibirsk::evaluate() const
{
  if (TMath::Abs(tail) < 1.e-7) {
    return TMath::Exp( -0.5 * TMath::Power( ( (x - peak) / width ), 2 ));
  }

  Double_t arg = 1.0 - ( x - peak ) * tail / width;

  if (arg < 1.e-7) {
    //Argument of logarithm negative. Real continuation -> function equals zero
    return 0.0;
  }

  Double_t log = TMath::Log(arg);
  static const Double_t xi = 2.3548200450309494; // 2 Sqrt( Ln(4) )

  Double_t width_zero = ( 2.0 / xi ) * TMath::ASinH( tail * xi * 0.5 );
  Double_t width_zero2 = width_zero * width_zero;
  Double_t exponent = ( -0.5 / (width_zero2) * log * log ) - ( width_zero2 * 0.5 );

  return TMath::Exp(exponent) ;
}
////////////////////////////////////////////////////////////////////////////////
/// Compute multiple values of Novosibirsk distribution.  
void RooNovosibirsk::computeBatch(double* output, size_t nEvents, rbc::DataMap& dataMap) const
{
  rbc::dispatch->compute(rbc::Novosibirsk, output, nEvents, dataMap, {&*x,&*peak,&*width,&*tail,&*_norm});
}

////////////////////////////////////////////////////////////////////////////////

Int_t RooNovosibirsk::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* ) const
{
  if (matchArgs(allVars,analVars,x)) return 1 ;
  if (matchArgs(allVars,analVars,peak)) return 2 ;

  //The other two integrals over tali and width are not integrable

  return 0 ;
}

////////////////////////////////////////////////////////////////////////////////

Double_t RooNovosibirsk::analyticalIntegral(Int_t code, const char* rangeName) const
{
  assert(code==1 || code==2) ;

  //The range is defined as [A,B]

  //Numerical values need for the evaluation of the integral
  static const Double_t sqrt2 = 1.4142135623730950; // Sqrt(2)
  static const Double_t sqlog2 = 0.832554611157697756; //Sqrt( Log(2) )
  static const Double_t sqlog4 = 1.17741002251547469; //Sqrt( Log(4) )
  static const Double_t log4 = 1.38629436111989062; //Log(2)
  static const Double_t rootpiby2 = 1.2533141373155003; // Sqrt(pi/2)
  static const Double_t sqpibylog2 = 2.12893403886245236; //Sqrt( pi/Log(2) )

  if (code==1) {
    Double_t A = x.min(rangeName);
    Double_t B = x.max(rangeName);

    Double_t result = 0;


    //If tail==0 the function becomes gaussian, thus we return a Gaussian integral
    if (TMath::Abs(tail) < 1.e-7) {

      Double_t xscale = sqrt2*width;

      result = rootpiby2*width*(TMath::Erf((B-peak)/xscale)-TMath::Erf((A-peak)/xscale));

      return result;

    }

    // If the range is not defined correctly the function becomes complex
    Double_t log_argument_A = ( (peak - A)*tail + width ) / width ;
    Double_t log_argument_B = ( (peak - B)*tail + width ) / width ;

    //lower limit
    if ( log_argument_A < 1.e-7) {
      log_argument_A = 1.e-7;
    }

    //upper limit
    if ( log_argument_B < 1.e-7) {
      log_argument_B = 1.e-7;
    }

    Double_t term1 =  TMath::ASinH( tail * sqlog4 );
    Double_t term1_2 =  term1 * term1;

    //Calculate the error function arguments
    Double_t erf_termA = ( term1_2 - log4 * TMath::Log( log_argument_A ) ) / ( 2 * term1 * sqlog2 );
    Double_t erf_termB = ( term1_2 - log4 * TMath::Log( log_argument_B ) ) / ( 2 * term1 * sqlog2 );

    result = 0.5 / tail * width * term1 * ( TMath::Erf(erf_termB) - TMath::Erf(erf_termA)) * sqpibylog2;

    return result;

  } else if (code==2) {
    Double_t A = x.min(rangeName);
    Double_t B = x.max(rangeName);

    Double_t result = 0;


    //If tail==0 the function becomes gaussian, thus we return a Gaussian integral
    if (TMath::Abs(tail) < 1.e-7) {

      Double_t xscale = sqrt2*width;

      result = rootpiby2*width*(TMath::Erf((B-x)/xscale)-TMath::Erf((A-x)/xscale));

      return result;

    }

    // If the range is not defined correctly the function becomes complex
    Double_t log_argument_A = ( (A - x)*tail + width ) / width;
    Double_t log_argument_B = ( (B - x)*tail + width ) / width;

    //lower limit
    if ( log_argument_A < 1.e-7) {
      log_argument_A = 1.e-7;
    }

    //upper limit
    if ( log_argument_B < 1.e-7) {
      log_argument_B = 1.e-7;
    }

    Double_t term1 =  TMath::ASinH( tail * sqlog4 );
    Double_t term1_2 =  term1 * term1;

    //Calculate the error function arguments
    Double_t erf_termA = ( term1_2 - log4 * TMath::Log( log_argument_A ) ) / ( 2 * term1 * sqlog2 );
    Double_t erf_termB = ( term1_2 - log4 * TMath::Log( log_argument_B ) ) / ( 2 * term1 * sqlog2 );

    result = 0.5 / tail * width * term1 * ( TMath::Erf(erf_termB) - TMath::Erf(erf_termA)) * sqpibylog2;

    return result;

  }

  // Emit fatal error
  coutF(Eval) << "Error in RooNovosibirsk::analyticalIntegral" << std::endl;

  // Put dummy return here to avoid compiler warnings
  return 1.0 ;
}
