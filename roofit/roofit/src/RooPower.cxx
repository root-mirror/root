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

/** \class RooPower
    \ingroup Roofit

RooPower implements a polynomial p.d.f of the form
\f[ f(x) = \mathcal{N} \cdot \sum_{i} a_{i} * x^i \f]
By default, the coefficient \f$ a_0 \f$ is chosen to be 1, as polynomial
probability density functions have one degree of freedom
less than polynomial functions due to the normalisation condition. \f$ \mathcal{N} \f$
is a normalisation constant that is automatically calculated when the polynomial is used
in computations.

The sum can be truncated at the low end. See the main constructor
RooPower::RooPower(const char*, const char*, RooAbsReal&, const RooArgList&, Int_t)
**/

#include <cmath>
#include <cassert>

#include "RooPower.h"
#include "RooAbsReal.h"
#include "RooArgList.h"
#include "RooMsgService.h"

#include "TError.h"

using namespace std;

ClassImp(RooPower);

////////////////////////////////////////////////////////////////////////////////
/// coverity[UNINIT_CTOR]

RooPower::RooPower()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Create a power law in the variable `x`.
/// \param[in] name Name of the PDF
/// \param[in] title Title for plotting the PDF
/// \param[in] x The variable of the polynomial
/// \param[in] coefList The coefficients \f$ a_i \f$
/// \param[in] expList The exponentials \f$ b_i \f$
/// \f[
///     1. + \sum_{i=0}^{n} a_{i} * x^{b_{i}}
/// \f]
///
/// This means that
/// \code{.cpp}
/// RooPower pol("pow", "pow", x, RooArgList(a1, a2), RooArgList(b1,b2))
/// \endcode
/// computes
/// \f[
///   \mathrm{pol}(x) = a1 * x^b1 + a2 * x^b2
/// \f]


RooPower::RooPower(const char* name, const char* title,
              RooAbsReal& x, const RooArgList& coefList, const RooArgList& expList) :
  RooAbsPdf(name, title),
  _x("x", "Dependent", this, x),
  _coefList("coefList","List of coefficients",this),
  _expList("expList","List of exponents",this)  
{
  RooFIter coefIter = coefList.fwdIterator() ;
  RooFIter expIter = expList.fwdIterator() ;
  RooAbsArg* coef, * exp ;
  while((coef = (RooAbsArg*)coefIter.next()) && (exp = (RooAbsArg*)expIter.next())) {
    if (!dynamic_cast<RooAbsReal*>(coef)) {
      coutE(InputArguments) << "RooPower::ctor(" << GetName() << ") ERROR: coefficient " << coef->GetName()
             << " is not of type RooAbsReal" << endl ;
      R__ASSERT(0) ;
    }
    if (!dynamic_cast<RooAbsReal*>(exp)) {
      coutE(InputArguments) << "RooPower::ctor(" << GetName() << ") ERROR: coefficient " << exp->GetName()
             << " is not of type RooAbsReal" << endl ;
      R__ASSERT(0) ;
    }
    _coefList.add(*coef) ;    
    _expList.add(*exp) ;
  }
}

////////////////////////////////////////////////////////////////////////////////

RooPower::RooPower(const char* name, const char* title,
                           RooAbsReal& x) :
  RooAbsPdf(name, title),
  _x("x", "Dependent", this, x),
  _coefList("coefList","List of coefficients",this),
  _expList("expList","List of exponents",this)  
{ }

////////////////////////////////////////////////////////////////////////////////
/// Copy constructor

RooPower::RooPower(const RooPower& other, const char* name) :
  RooAbsPdf(other, name),
  _x("x", this, other._x),
  _coefList("coefList",this,other._coefList),
  _expList("expList",this,other._expList)  
{ }

////////////////////////////////////////////////////////////////////////////////
/// Destructor

RooPower::~RooPower()
{ }

////////////////////////////////////////////////////////////////////////////////

Double_t RooPower::evaluate() const
{
  // Calculate and return value of polynomial

  const unsigned sz = _coefList.getSize();
  if (!sz) return 0.;

  std::vector<double> coefs;
  std::vector<double> exps;  
  coefs.reserve(sz);
  exps.reserve(sz);  
  {
    const RooArgSet* nset = _coefList.nset();
    RooAbsReal* c;
    RooFIter coef_it = _coefList.fwdIterator();
    while ((c = (RooAbsReal*) coef_it.next())) coefs.push_back(c->getVal(nset));
    RooFIter exp_it = _expList.fwdIterator();
    while ((c = (RooAbsReal*) exp_it.next())) exps.push_back(c->getVal(nset));
    
  }
  double x = this->_x;
  Double_t retval = 0;
  for(unsigned int i=0; i<sz; ++i){
    retval += coefs[i] * pow(x,exps[i]);
  }
  return retval;
}
