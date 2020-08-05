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

RooPower implements a power law p.d.f of the form
\f[ f(x) = \mathcal{N} \cdot \sum_{i} a_{i} * x^{b_i} \f]

\image html RooPower.png
**/

#include <cmath>
#include <cassert>
#include <sstream>

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
///     \sum_{i=0}^{n} a_{i} * x^{b_{i}}
/// \f]
///
/// This means that
/// \code{.cpp}
/// RooPower powl("pow", "pow", x, RooArgList(a1, a2), RooArgList(b1,b2))
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
  if(coefList.getSize() != expList.getSize()){
    coutE(InputArguments) << "RooPower::ctor(" << GetName() << ") ERROR: coefficient list and exponent list must be of same length" << endl;
    return;
  }
  for(auto coef : coefList){
    if (!dynamic_cast<RooAbsReal*>(coef)) {
      coutE(InputArguments) << "RooPower::ctor(" << GetName() << ") ERROR: coefficient " << coef->GetName()
                            << " is not of type RooAbsReal" << endl ;
      R__ASSERT(0) ;
    }
    _coefList.add(*coef) ;    
  }
  for(auto exp : expList){
    if (!dynamic_cast<RooAbsReal*>(exp)) {
      coutE(InputArguments) << "RooPower::ctor(" << GetName() << ") ERROR: coefficient " << exp->GetName()
             << " is not of type RooAbsReal" << endl ;
      R__ASSERT(0) ;
    }
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
  const RooArgSet* nset = _coefList.nset();
  for(auto c:_coefList){ coefs.push_back(static_cast<RooAbsReal*>(c)->getVal(nset)); }
  for(auto c:_expList ) { exps.push_back(static_cast<RooAbsReal*>(c)->getVal(nset)); }    
  double x = this->_x;
  Double_t retval = 0;
  for(unsigned int i=0; i<sz; ++i){
    retval += coefs[i] * pow(x,exps[i]);
  }
  return retval;
}

////////////////////////////////////////////////////////////////////////////////
/// Advertise to RooFit that this function can be analytically integrated.
Int_t RooPower::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const
{
  if (matchArgs(allVars, analVars, _x)) return 1;
  return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Do the analytical integral according to the code that was returned by getAnalyticalIntegral().
Double_t RooPower::analyticalIntegral(Int_t code, const char* rangeName) const
{
  R__ASSERT(code==1) ;
  
  const Double_t xmin = _x.min(rangeName), xmax = _x.max(rangeName);
  const unsigned sz = _coefList.getSize();
  if (!sz) return xmax - xmin;
  
  std::vector<double> coefs;
  std::vector<double> exps;  
  coefs.reserve(sz);
  exps.reserve(sz);  
  const RooArgSet* nset = _coefList.nset();
  for(auto c:_coefList){ coefs.push_back(static_cast<RooAbsReal*>(c)->getVal(nset)); }
  for(auto c:_expList ) { exps.push_back(static_cast<RooAbsReal*>(c)->getVal(nset)); }    
  
  Double_t retval = 0;
  for(unsigned int i=0; i<sz; ++i){
    if(exps[i] == -1){
      retval += coefs[i] * (log(xmax)-log(xmin));
    } else {
      retval += coefs[i]/(exps[i]+1) * (pow(xmax,(exps[i]+1))-pow(xmin,(exps[i]+1)));
    }
  }
  return retval;
}

std::string RooPower::getFormulaExpression(bool expand) const {
  std::stringstream ss;
  for(int i=0; i<_coefList.getSize(); ++i){
    if(i!=0) ss << "+";
    if(expand) ss << ((RooAbsReal*)_coefList.at(i))->getVal();
    else       ss << _coefList.at(i)->GetName();
    ss << "*pow(" << _x.GetName() << "," << _expList.at(i)->GetName() << ")";
  }
  return ss.str().c_str();
}
