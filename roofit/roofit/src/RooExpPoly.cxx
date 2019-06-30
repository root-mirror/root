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
**/

#include <cmath>
#include <cassert>

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

  RooFIter coefIter = coefList.fwdIterator() ;
  RooAbsArg* coef ;
  while((coef = (RooAbsArg*)coefIter.next())) {
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

Double_t RooExpPoly::evaluate() const
{
  // Calculate and return value of polynomial

  const unsigned sz = _coefList.getSize();
  const int lowestOrder = _lowestOrder;
  if (!sz) return lowestOrder ? 1. : 0.;
  std::vector<double> coefs;
  coefs.reserve(sz);
  {
    const RooArgSet* nset = _coefList.nset();
    RooFIter it = _coefList.fwdIterator();
    RooAbsReal* c;
    while ((c = (RooAbsReal*) it.next())) coefs.push_back(c->getVal(nset));
  }
  const Double_t x = _x;
  double xpow = std::pow(x, lowestOrder);
  double retval = 0;
  for(size_t i=0; i<sz; ++i){
    retval += coefs[i]*xpow;
    xpow *= x;
  }
  return exp(retval);
}
