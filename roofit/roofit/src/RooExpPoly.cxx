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

/** \class RooExpPoly
    \ingroup Roofit

RooExpPoly implements a polynomial p.d.f of the form \f[ f(x) =
\mathcal{N} \cdot \exp( \sum_{i} a_{i} * x^{i} ) \f] \f$ \mathcal{N}
\f$ is a normalisation constant that is automatically calculated when
the function is used in computations.

The sum can be truncated at the low end. See the main constructor
RooExpPoly::RooExpPoly(const char*, const char*, RooAbsReal&, const RooArgList&, Int_t)

\image html RooExpPoly.png

**/

#include <cmath>
#include <sstream>
#include <cassert>

#include "TMath.h"

#include "RooExpPoly.h"
#include "RooAbsReal.h"
#include "RooArgList.h"
#include "RooMsgService.h"

#include "TError.h"

using namespace std;

ClassImp(RooExpPoly);

////////////////////////////////////////////////////////////////////////////////
/// coverity[UNINIT_CTOR]

RooExpPoly::RooExpPoly()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Create a polynomial in the variable `x`.
/// \param[in] name Name of the PDF
/// \param[in] title Title for plotting the PDF
/// \param[in] x The variable of the polynomial
/// \param[in] coefList The coefficients \f$ a_i \f$
/// \param[in] lowestOrder [optional] Truncate the sum such that it skips the lower orders:
/// \f[
///     1. + \sum_{i=0}^{\mathrm{coefList.size()}} a_{i} * x^{(i + \mathrm{lowestOrder})}
/// \f]
///
/// This means that
/// \code{.cpp}
/// RooExpPoly pol("pol", "pol", x, RooArgList(a, b), lowestOrder = 2)
/// \endcode
/// computes
/// \f[
///   \mathrm{pol}(x) = 1 * x^0 + (0 * x^{\ldots}) + a * x^2 + b * x^3.
/// \f]


RooExpPoly::RooExpPoly(const char* name, const char* title,
                       RooAbsReal& x, const RooArgList& coefList, Int_t lowestOrder) :
  RooAbsPdf(name, title),
  _x("x", "Dependent", this, x),
  _coefList("coefList","List of coefficients",this),
  _lowestOrder(lowestOrder)
{
  // Check lowest order
  if (_lowestOrder<0) {
    coutE(InputArguments) << "RooExpPoly::ctor(" << GetName()
           << ") WARNING: lowestOrder must be >=0, setting value to 0" << endl ;
    _lowestOrder=0 ;
  }

  for(auto coef:coefList){
    if (!dynamic_cast<RooAbsReal*>(coef)) {
      coutE(InputArguments) << "RooExpPoly::ctor(" << GetName() << ") ERROR: coefficient " << coef->GetName()
                            << " is not of type RooAbsReal" << endl ;
      R__ASSERT(0) ;
    }
    _coefList.add(*coef) ;
  }
}

////////////////////////////////////////////////////////////////////////////////

RooExpPoly::RooExpPoly(const char* name, const char* title,
                           RooAbsReal& x) :
  RooAbsPdf(name, title),
  _x("x", "Dependent", this, x),
  _coefList("coefList","List of coefficients",this),
  _lowestOrder(1)
{ }

////////////////////////////////////////////////////////////////////////////////
/// Copy constructor

RooExpPoly::RooExpPoly(const RooExpPoly& other, const char* name) :
  RooAbsPdf(other, name),
  _x("x", this, other._x),
  _coefList("coefList",this,other._coefList),
  _lowestOrder(other._lowestOrder)
{ }

////////////////////////////////////////////////////////////////////////////////
/// Destructor

RooExpPoly::~RooExpPoly()
{ }

////////////////////////////////////////////////////////////////////////////////

Double_t RooExpPoly::evaluateLog() const
{
  // Calculate and return value of polynomial
  
  const unsigned sz = _coefList.getSize();
  const int lowestOrder = _lowestOrder;
  if (!sz) return lowestOrder ? 1. : 0.;
  std::vector<double> coefs;
  coefs.reserve(sz);

  const RooArgSet* nset = _coefList.nset();
  for(auto coef:_coefList) { coefs.push_back(static_cast<RooAbsReal*>(coef)->getVal(nset)); };
  const Double_t x = _x;
  double xpow = std::pow(x, lowestOrder);
  double retval = 0;
  for(size_t i=0; i<sz; ++i){
    retval += coefs[i]*xpow;
    xpow *= x;
  }
  return retval;
}

////////////////////////////////////////////////////////////////////////////////

Double_t RooExpPoly::evaluate() const {
  return exp(this->evaluateLog());
}

////////////////////////////////////////////////////////////////////////////////

Double_t RooExpPoly::getLogVal(const RooArgSet* nset) const 
{
  return RooAbsPdf::getLogVal(nset);
  //  return this->evaluateLog();
}

////////////////////////////////////////////////////////////////////////////////

Int_t RooExpPoly::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const
{
  
  if ((_coefList.getSize() + _lowestOrder < 4) && (( _coefList.getSize() + _lowestOrder < 3) || (static_cast<RooAbsReal*>(_coefList.at(2-_lowestOrder))->getVal() <= 0)) && matchArgs(allVars,analVars,_x)){
    return 1;
  }
  return 0 ;
}

////////////////////////////////////////////////////////////////////////////////

#define PI TMath::Pi()

namespace {
  double deltaerf(long double x1, long double x2){
    // several things happening here
    // 1. use "erf(x) = -erf(-x)" to evaluate the function only ever for the positive side (higher precision)
    // 2. use "erf(x) = 1.-erfc(x)" and, instead of "erf(x1) - erf(x2)", do "(1.-erfc(x1)) - (1.-erfc(x2)) = erfc(x2) - erfc(x1)"
    double y2 = x1 > 0 ? erfc(x1) : -erfc(-x1);
    double y1 = x2 > 0 ? erfc(x2) : -erfc(-x2);
    if(y1 == y2){
      std::cout << "WARNING in calculation of analytical integral limited by numerical precision" << std::endl;
      std::cout << "x: " << x1 << " , " << x2 << std::endl;
      std::cout << "y: " << y1 << " , " << y2 << std::endl;
    } 
    return y1 - y2;
  }
}

Double_t RooExpPoly::analyticalIntegral(Int_t code, const char* rangeName) const
{
  assert(code==1);

  const Double_t xmin = _x.min(rangeName), xmax = _x.max(rangeName);
  const unsigned sz = _coefList.getSize();
  if (!sz) return xmax - xmin;
  
  std::vector<double> coefs;
  coefs.reserve(sz);
  const RooArgSet* nset = _coefList.nset();
  for(auto c:_coefList){ coefs.push_back(static_cast<RooAbsReal*>(c)->getVal(nset)); }
  
  switch(_coefList.getSize()+_lowestOrder){
  case 1:
    return xmax - xmin;    
  case 2: {
    const double a = coefs[1-_lowestOrder];
    if(a!=0){
      return 1./a * (exp(a*xmax)-exp(a*xmin)) * (_lowestOrder == 0 ? exp(coefs[0]) : 1 );
    } else {
      return xmax - xmin;
    }
  }
  case 3: {
    const double a = coefs[2-_lowestOrder];
    const double b = _lowestOrder == 2 ? 0. : coefs[1-_lowestOrder];
    const double c = _lowestOrder == 0 ? coefs[0] : 0.;
    const double absa = std::abs(a);
    const double sqrta = sqrt(absa);
    if(a < 0){
      double d = ::deltaerf((-b + 2 * absa * xmax)/(2 * sqrta),(-b + 2 * absa * xmin)/(2 * sqrta));
      double retval = exp(b*b/(4 * absa) + c) * sqrt(PI) * d /(2 * sqrta);
      return retval;  
    } else if(a>0){
      throw std::runtime_error("dawson integral not yet implemented!");
    } else if(b!=0){
      return 1./b * (exp(b*xmax)-exp(b*xmin)) * exp(c);
    } else {      
      return xmax - xmin;
    }
  }
  }
  return 0.;
}

std::string RooExpPoly::getFormulaExpression(bool expand) const {
  std::stringstream ss;
  ss << "exp(" ;
  int order = _lowestOrder;
  for(auto coef:_coefList) {
    if(order != _lowestOrder) ss << "+";
    if(expand) ss << ((RooAbsReal*)coef)->getVal();
    else       ss << coef->GetName();
    ss << "*pow(" << _x.GetName() << "," << order << ")";
    ++order;
  }
  ss << ")";
  return ss.str();
}
