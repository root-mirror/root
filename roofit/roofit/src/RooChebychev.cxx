/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 * @(#)root/roofit:$Id$
 * Authors:                                                                  *
 *   GR, Gerhard Raven,   UC San Diego, Gerhard.Raven@slac.stanford.edu
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

/** \class RooChebychev
   \ingroup Roofit

Chebychev polynomial p.d.f. of the first kind.

The coefficient that goes with \f$ T_0(x)=1 \f$ (i.e. the constant polynomial) is
implicitly assumed to be 1, and the list of coefficients supplied by callers
starts with the coefficient that goes with \f$ T_1(x)=x \f$ (i.e. the linear term).
**/

#include "RooChebychev.h"
#include "RooFit.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooArgList.h"
#include "RooNameReg.h"
#include "RooBatchCompute.h"

#include <cmath>

ClassImp(RooChebychev);

namespace { // anonymous namespace to hide implementation details
/// use fast FMA if available, fall back to normal arithmetic if not
static inline double fast_fma(
	const double x, const double y, const double z) noexcept
{
#if defined(FP_FAST_FMA) // check if std::fma has fast hardware implementation
   return std::fma(x, y, z);
#else // defined(FP_FAST_FMA)
   // std::fma might be slow, so use a more pedestrian implementation
#if defined(__clang__)
#pragma STDC FP_CONTRACT ON // hint clang that using an FMA is okay here
#endif // defined(__clang__)
   return (x * y) + z;
#endif // defined(FP_FAST_FMA)
}

/// Chebychev polynomials of first or second kind
enum class Kind : int { First = 1, Second = 2 };

/** @brief ChebychevIterator evaluates increasing orders at given x
 *
 * @author Manuel Schiller <Manuel.Schiller@glasgow.ac.uk>
 * @date 2019-03-24
 */
template <typename T, Kind KIND>
class ChebychevIterator {
private:
   T _last = 1;
   T _curr = 0;
   T _twox = 0;

public:
   /// default constructor
   constexpr ChebychevIterator() = default;
   /// copy constructor
   ChebychevIterator(const ChebychevIterator &) = default;
   /// move constructor
   ChebychevIterator(ChebychevIterator &&) = default;
   /// construct from given x in [-1, 1]
   constexpr ChebychevIterator(const T &x)
       : _curr(static_cast<int>(KIND) * x), _twox(2 * x)
   {}

   /// (copy) assignment
   ChebychevIterator &operator=(const ChebychevIterator &) = default;
   /// move assignment
   ChebychevIterator &operator=(ChebychevIterator &&) = default;

   /// get value of Chebychev polynomial at current order
   constexpr inline T operator*() const noexcept { return _last; }
   // get value of Chebychev polynomial at (current + 1) order
   constexpr inline T lookahead() const noexcept { return _curr; }
   /// move on to next order, return reference to new value
   inline ChebychevIterator &operator++() noexcept
   {
      //T newval = fast_fma(_twox, _curr, -_last);
      T newval = _twox*_curr -_last;
      _last = _curr;
      _curr = newval;
      return *this;
   }
   /// move on to next order, return copy of new value
   inline ChebychevIterator operator++(int) noexcept
   {
      ChebychevIterator retVal(*this);
      operator++();
      return retVal;
   }
};
} // anonymous namespace

////////////////////////////////////////////////////////////////////////////////

RooChebychev::RooChebychev() : _refRangeName(0)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor

RooChebychev::RooChebychev(const char* name, const char* title,
                           RooAbsReal& x, const RooArgList& coefList):
  RooAbsPdf(name, title),
  _x("x", "Dependent", this, x),
  _coefList("coefficients","List of coefficients",this),
  _refRangeName(0)
{
  for (const auto coef : coefList) {
    if (!dynamic_cast<RooAbsReal*>(coef)) {
      coutE(InputArguments) << "RooChebychev::ctor(" << GetName() <<
       ") ERROR: coefficient " << coef->GetName() <<
       " is not of type RooAbsReal" << std::endl ;
      throw std::invalid_argument("Wrong input arguments for RooChebychev");
    }
    _coefList.add(*coef) ;
  }
}

////////////////////////////////////////////////////////////////////////////////

RooChebychev::RooChebychev(const RooChebychev& other, const char* name) :
  RooAbsPdf(other, name),
  _x("x", this, other._x),
  _coefList("coefList",this,other._coefList),
  _refRangeName(other._refRangeName)
{
}

////////////////////////////////////////////////////////////////////////////////

void RooChebychev::selectNormalizationRange(const char* rangeName, Bool_t force)
{
  if (rangeName && (force || !_refRangeName)) {
    _refRangeName = (TNamed*) RooNameReg::instance().constPtr(rangeName) ;
  }
  if (!rangeName) {
    _refRangeName = 0 ;
  }
}

////////////////////////////////////////////////////////////////////////////////

Double_t RooChebychev::evaluate() const
{
  // first bring the range of the variable _x to the normalised range [-1, 1]
  // calculate sum_k c_k T_k(x) where x is given in the normalised range,
  // c_0 = 1, and the higher coefficients are given in _coefList
  const Double_t xmax = _x.max(_refRangeName?_refRangeName->GetName():0);
  const Double_t xmin = _x.min(_refRangeName?_refRangeName->GetName():0);
  // transform to range [-1, +1]
  const Double_t x = (_x - 0.5 * (xmax + xmin)) / (0.5 * (xmax - xmin));
  // extract current values of coefficients
  using size_type = typename RooListProxy::Storage_t::size_type;
  const size_type iend = _coefList.size();
  double sum = 1.;
  if (iend > 0) {
     ChebychevIterator<double, Kind::First> chit(x);
     ++chit;
     for (size_type i = 0; iend != i; ++i, ++chit) {
        auto c = static_cast<const RooAbsReal &>(_coefList[i]).getVal();
        //sum = fast_fma(*chit, c, sum);
        sum += *chit*c;
     }
  }
  return sum;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute multiple values of Chebychev.  
RooSpan<double> RooChebychev::evaluateSpan(RooBatchCompute::RunContext& evalData, const RooArgSet* normSet) const {

  RooSpan<const double> xData = _x->getValues(evalData, normSet);  
  size_t batchSize = xData.size();
  RooSpan<double> output = evalData.makeBatch(this, batchSize);
  const Double_t xmin = _x.min(_refRangeName?_refRangeName->GetName() : nullptr);
  const Double_t xmax = _x.max(_refRangeName?_refRangeName->GetName() : nullptr);

  const size_t nCoef = _coefList.size();
  std::vector<double> coef(nCoef);
  for (size_t i=0; i<nCoef; i++) {
    coef[i] = static_cast<const RooAbsReal &>(_coefList[i]).getVal();
  }
  RooBatchCompute::dispatch->computeChebychev(batchSize, output.data(), xData.data(), xmin, xmax, coef);
  return output;
}

////////////////////////////////////////////////////////////////////////////////


Int_t RooChebychev::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /* rangeName */) const
{
  if (matchArgs(allVars, analVars, _x)) return 1;
  return 0;
}

////////////////////////////////////////////////////////////////////////////////

Double_t RooChebychev::analyticalIntegral(Int_t code, const char* rangeName) const
{
  assert(1 == code); (void)code;

  const Double_t xmax = _x.max(_refRangeName?_refRangeName->GetName():0);
  const Double_t xmin = _x.min(_refRangeName?_refRangeName->GetName():0);
  const Double_t halfrange = .5 * (xmax - xmin), mid = .5 * (xmax + xmin);
  // the full range of the function is mapped to the normalised [-1, 1] range
  const Double_t b = (_x.max(rangeName) - mid) / halfrange;
  const Double_t a = (_x.min(rangeName) - mid) / halfrange;

  // take care to multiply with the right factor to account for the mapping to
  // normalised range [-1, 1]
  return halfrange * evalAnaInt(a, b);
}

////////////////////////////////////////////////////////////////////////////////

Double_t RooChebychev::evalAnaInt(const Double_t a, const Double_t b) const
{
   // coefficient for integral(T_0(x)) is 1 (implicit), integrate by hand
   // T_0(x) and T_1(x), and use for n > 1: integral(T_n(x) dx) =
   // (T_n+1(x) / (n + 1) - T_n-1(x) / (n - 1)) / 2
   double sum = b - a; // integrate T_0(x) by hand

   using size_type = typename RooListProxy::Storage_t::size_type;
   const size_type iend = _coefList.size();
   if (iend > 0) {
      {
         // integrate T_1(x) by hand...
         const double c = static_cast<const RooAbsReal &>(_coefList[0]).getVal();
         sum = fast_fma(0.5 * (b + a) * (b - a), c, sum);
      }
      if (1 < iend) {
         ChebychevIterator<double, Kind::First> bit(b), ait(a);
         ++bit, ++ait;
         double nminus1 = 1.;
         for (size_type i = 1; iend != i; ++i) {
            // integrate using recursion relation
            const double c = static_cast<const RooAbsReal &>(_coefList[i]).getVal();
            const double term2 = (*bit - *ait) / nminus1;
            ++bit, ++ait, ++nminus1;
            const double term1 = (bit.lookahead() - ait.lookahead()) / (nminus1 + 1.);
            const double intTn = 0.5 * (term1 - term2);
            sum = fast_fma(intTn, c, sum);
         }
      }
  }
  return sum;
}
