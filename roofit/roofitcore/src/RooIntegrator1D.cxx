/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 * @(#)root/roofitcore:$Id$
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

/**
\file RooIntegrator1D.cxx
\class RooIntegrator1D
\ingroup Roofitcore

RooIntegrator1D implements an adaptive one-dimensional 
numerical integration algorithm. 

It uses Romberg's method with trapezoids or midpoints.
The integrand is approximated by \f$ 1, 2, 4, 8, \ldots, 2^n \f$ trapezoids, and
Richardson series acceleration is applied to the series. For refinement step \f$ n \f$, that is
\f[
  R(n,m) = \frac{1}{4^m - 1} \left( 4^m R(n, m-1) - R(n-1, m-1) \right)
\f]

The integrator will evaluate the first six refinements (i.e. 64 points) in one go,
apply a five-point series acceleration, and successively add more steps until the
desired precision is reached.

\since In ROOT 6.24, the implementation of the integrator was revised, since it often
stopped early, not reaching the desired relative precision. The old (less accurate) integrator
is available under the name OldIntegrator1D. If less precision is actually desired (to speed up the
integration), a relative epsilon 5, 10 or more times higher than for the old integrator can be used.
**/


#include "RooIntegrator1D.h"

#include "RooArgSet.h"
#include "RooRealVar.h"
#include "RooIntegratorBinding.h"
#include "RooNumIntConfig.h"
#include "RooNumIntFactory.h"
#include "RooMsgService.h"

#include "Math/Util.h"

#include <numeric>
#include <array>
#include <assert.h>
#include <iomanip>

using namespace std;

ClassImp(RooIntegrator1D);


////////////////////////////////////////////////////////////////////////////////
/// Register RooIntegrator1D, its parameters and capabilities with RooNumIntFactory

void RooIntegrator1D::registerIntegrator(RooNumIntFactory& fact)
{
  RooCategory sumRule("sumRule","Summation Rule") ;
  sumRule.defineType("Trapezoid",RooIntegrator1D::Trapezoid) ;
  sumRule.defineType("Midpoint",RooIntegrator1D::Midpoint) ;
  sumRule.setLabel("Trapezoid") ;
  RooCategory extrap("extrapolation","Extrapolation procedure") ;
  extrap.defineType("None",0) ;
  extrap.defineType("Wynn-Epsilon",1) ;
  extrap.setLabel("Wynn-Epsilon") ;
  RooRealVar maxSteps("maxSteps","Maximum number of steps",20) ;
  RooRealVar minSteps("minSteps","Minimum number of steps",999) ;
  RooRealVar fixSteps("fixSteps","Fixed number of steps",0) ;

  RooIntegrator1D* proto = new RooIntegrator1D() ;
  fact.storeProtoIntegrator(proto,RooArgSet(sumRule,extrap,maxSteps,minSteps,fixSteps)) ;
  RooNumIntConfig::defaultConfig().method1D().setLabel(proto->IsA()->GetName()) ;
}


////////////////////////////////////////////////////////////////////////////////
/// Construct integrator on given function binding, using specified summation
/// rule, maximum number of steps and conversion tolerance. The integration
/// limits are taken from the function binding

RooIntegrator1D::RooIntegrator1D(const RooAbsFunc& function, SummationRule rule,
				 Int_t maxSteps, Double_t eps) : 
  RooAbsIntegrator(function), _rule(rule), _maxSteps(maxSteps),  _minStepsZero(999), _fixSteps(0), _epsAbs(eps), _epsRel(eps), _doExtrap(kTRUE)
{
  _useIntegrandLimits= kTRUE;
  _valid= initialize();
} 


////////////////////////////////////////////////////////////////////////////////
/// Construct integrator on given function binding for given range,
/// using specified summation rule, maximum number of steps and
/// conversion tolerance. The integration limits are taken from the
/// function binding

RooIntegrator1D::RooIntegrator1D(const RooAbsFunc& function, Double_t xmin, Double_t xmax,
				 SummationRule rule, Int_t maxSteps, Double_t eps) : 
  RooAbsIntegrator(function), 
  _rule(rule), 
  _maxSteps(maxSteps), 
  _minStepsZero(999), 
  _fixSteps(0),
  _epsAbs(eps), 
  _epsRel(eps),
  _doExtrap(kTRUE)
{
  _useIntegrandLimits= kFALSE;
  _xmin= xmin;
  _xmax= xmax;
  _valid= initialize();
} 


////////////////////////////////////////////////////////////////////////////////
/// Construct integrator on given function binding, using specified
/// configuration object. The integration limits are taken from the
/// function binding

RooIntegrator1D::RooIntegrator1D(const RooAbsFunc& function, const RooNumIntConfig& config) :
  RooAbsIntegrator(function,config.printEvalCounter()), 
  _epsAbs(config.epsAbs()),
  _epsRel(config.epsRel())
{
  // Extract parameters from config object
  const RooArgSet& configSet = config.getConfigSection(IsA()->GetName()) ;  
  _rule = (SummationRule) configSet.getCatIndex("sumRule",Trapezoid) ;
  _maxSteps = (Int_t) configSet.getRealValue("maxSteps",20) ;
  _minStepsZero = (Int_t) configSet.getRealValue("minSteps",999) ;
  _fixSteps = (Int_t) configSet.getRealValue("fixSteps",0) ;
  _doExtrap = (Bool_t) configSet.getCatIndex("extrapolation",1) ;

  if (_fixSteps>_maxSteps) {
    oocoutE((TObject*)0,Integration) << "RooIntegrator1D::ctor() ERROR: fixSteps>maxSteps, fixSteps set to maxSteps" << endl ;
    _fixSteps = _maxSteps ;
  }

  _useIntegrandLimits= kTRUE;
  _valid= initialize();
} 



////////////////////////////////////////////////////////////////////////////////
/// Construct integrator on given function binding, using specified
/// configuration object and integration range

RooIntegrator1D::RooIntegrator1D(const RooAbsFunc& function, Double_t xmin, Double_t xmax,
				const RooNumIntConfig& config) :
  RooAbsIntegrator(function,config.printEvalCounter()), 
  _epsAbs(config.epsAbs()),
  _epsRel(config.epsRel())
{
  // Extract parameters from config object
  const RooArgSet& configSet = config.getConfigSection(IsA()->GetName()) ;  
  _rule = (SummationRule) configSet.getCatIndex("sumRule",Trapezoid) ;
  _maxSteps = (Int_t) configSet.getRealValue("maxSteps",20) ;
  _minStepsZero = (Int_t) configSet.getRealValue("minSteps",999) ;
  _fixSteps = (Int_t) configSet.getRealValue("fixSteps",0) ;  
  _doExtrap = (Bool_t) configSet.getCatIndex("extrapolation",1) ;

  _useIntegrandLimits= kFALSE;
  _xmin= xmin;
  _xmax= xmax;
  _valid= initialize();
} 



////////////////////////////////////////////////////////////////////////////////
/// Clone integrator with new function binding and configuration. Needed by RooNumIntFactory

RooAbsIntegrator* RooIntegrator1D::clone(const RooAbsFunc& function, const RooNumIntConfig& config) const
{
  return new RooIntegrator1D(function,config) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Initialize the integrator

Bool_t RooIntegrator1D::initialize()
{
  // apply defaults if necessary
  if(_maxSteps <= 0) {
    _maxSteps= (_rule == Trapezoid) ? 20 : 14;
  }

  if(_epsRel <= 0) _epsRel= 1e-6;
  if(_epsAbs <= 0) _epsAbs= 1e-6;

  // check that the integrand is a valid function
  if(!isValid()) {
    oocoutE((TObject*)0,Integration) << "RooIntegrator1D::initialize: cannot integrate invalid function" << endl;
    return kFALSE;
  }


  // Allocate workspace for numerical integration engine
  _s.resize(_maxSteps);

  return checkLimits();
}


////////////////////////////////////////////////////////////////////////////////
/// Change our integration limits. Return kTRUE if the new limits are
/// ok, or otherwise kFALSE. Always returns kFALSE and does nothing
/// if this object was constructed to always use our integrand's limits.

Bool_t RooIntegrator1D::setLimits(Double_t *xmin, Double_t *xmax) 
{
  if(_useIntegrandLimits) {
    oocoutE((TObject*)0,Integration) << "RooIntegrator1D::setLimits: cannot override integrand's limits" << endl;
    return kFALSE;
  }
  _xmin= *xmin;
  _xmax= *xmax;

  return checkLimits();
}


////////////////////////////////////////////////////////////////////////////////
/// Check that our integration range is finite and otherwise return kFALSE.
/// Update the limits from the integrand if requested.

Bool_t RooIntegrator1D::checkLimits() const 
{
  if(_useIntegrandLimits) {
    assert(0 != integrand() && integrand()->isValid());
    const_cast<double&>(_xmin) = integrand()->getMinLimit(0);
    const_cast<double&>(_xmax) = integrand()->getMaxLimit(0);
  }

  if (_xmax - _xmin < 0.) {
    oocoutE((TObject*)0,Integration) << "RooIntegrator1D::checkLimits: bad range with min > max (_xmin = " << _xmin << " _xmax = " << _xmax << ")" << endl;
    return kFALSE;
  }
  return !(RooNumber::isInfinite(_xmin) || RooNumber::isInfinite(_xmax));
}


namespace {
////////////////////////////////////////////////////////////////////////////////
/// Apply Richardson series acceleration to trapezoid / midpoint sums.
/// \tparam[in] numAccel Number of series acceleration steps to apply.
/// \param[in] nSeries Current refinement step. This number or more elements of `series` have been computed.
/// \param[in] series Series to accelerate.
/// \param[in] printStr Print series acceleration computation to this stream.
template<int numAccel>
std::pair<double,double> richardsonExtrapolation(int nSeries, const std::vector<double>& series, std::ostream* printStr = nullptr) {
  static_assert(numAccel >= 1, "Series acceleration only makes sense when applied at least once.");
  assert(nSeries <= static_cast<int>(series.size()));
  assert(nSeries > 0);
  if (nSeries == 1)
    return {series[0], series[0]};

  std::vector<std::vector<double>> R(nSeries);

  for (int n = 0; n < static_cast<int>(R.size()); ++n) {
    R[n].push_back(series[n]);
    unsigned int fourToTheMth = 1;

    for (int m = 1; m < numAccel && m <= n && m <= numAccel - (nSeries - n); ++m) {
      fourToTheMth *= 4;

      // compute R[n][m]
      R[n].push_back( 1./(fourToTheMth - 1.) * ( fourToTheMth * R[n][m-1] - R[n-1][m-1] ) );
    }

    if (printStr) {
      *printStr << "Integration step n=" << std::setw(2) << std::right << n;
      for (unsigned int m=0; m < R[n].size(); ++m)
        *printStr << "  m=" << m << ": " << std::setw(14) << std::setprecision(9) << R[n][m];
      *printStr << std::endl;
    }
  }

  const double ret = R.back().back();
  const double prev = R[R.size()-2].back();

  if (printStr)
    *printStr << "Integral = " << ret << "\t+/-\t" << fabs( ret - prev) << std::endl;

  return {ret, fabs(ret - prev)};
}

}


////////////////////////////////////////////////////////////////////////////////
/// Calculate numeric integral at given set of function parameters.
/// \param xValues Values at which to evaluate the function.
/// \param parameters If the function is more than 1-dimensional, the parameters of the function
/// can be set using this array.
/// \param nPar Size of the parameter array.
/// Parameters that are not passed in the `parameters` array will be taken at the current values
/// of the parameter objects.
RooSpan<const double> RooIntegrator1D::evalIntegrand(const std::vector<double>& xValues, const double* parameters, std::size_t nPar) const {
  std::vector<RooSpan<const double>> arguments;
  arguments.emplace_back(xValues);

  for (unsigned int i=0; i < nPar; ++i) {
    arguments.emplace_back(&parameters[i], 1);
  }

  auto results = _function->getValues(arguments);

  return results;
}


////////////////////////////////////////////////////////////////////////////////
/// Calculate numeric integral at given set of function binding parameters.
/// If the function is more than 1-dimensional, the parameters for dimensions
/// 1 to N can be passed in the `yvec`.
/// \warning N-1 parameters will be copied from `yvec`, irrespective of its true size.
Double_t RooIntegrator1D::integral(const Double_t *yvec) {
  assert(isValid());

  const std::size_t nPar = yvec != nullptr ? _function->getDimension() - 1 : 0;

  if (_xmax - _xmin == 0.)
    return 0.;

  const unsigned int initialSteps = std::min(6, _maxSteps - 1); // 65 points
  constexpr unsigned int seriesAccelerationSteps = 5;
  _s = _rule == Trapezoid ? computeTrapezoids(0, initialSteps, 0., yvec, nPar) : computeMidpoints(0, initialSteps, 0., yvec, nPar);
  _s.resize(_maxSteps, 0.);

  const double zeroThresh = _epsAbs/(_xmax - _xmin);

  for (unsigned int j = initialSteps; j < static_cast<unsigned int>(_maxSteps); ++j) {
    // refine our estimate using the appropriate summation rule
    _s[j] = (_rule == Trapezoid) ?
        computeTrapezoids(j, j+1, j > 0 ? _s[j-1] : 0., yvec, nPar).front() :
        computeMidpoints( j, j+1, j > 0 ? _s[j-1] : 0., yvec, nPar).front();

    if (j > static_cast<unsigned int>(_minStepsZero)
        && std::all_of(_s.begin(), _s.begin()+j+1, [=](double val){ return val < zeroThresh; })) {
      // All sums are so close to zero that we can stop here
      return 0.;
    }

    if (_fixSteps > 0 && j+1 == static_cast<unsigned int>(_fixSteps)) {
      return _s[j];
    }

    if(j > seriesAccelerationSteps) {
      double extrapValue, extrapError;

      // extrapolate the results of recent refinements and check for a stable result
      if (_doExtrap) {
        std::tie(extrapValue, extrapError) = richardsonExtrapolation<seriesAccelerationSteps>(j, _s);
      } else {
        extrapValue = _s[j];
        extrapError = _s[j] - _s[j-1];
      }

      // The error estimated from comparing the last two points of the series
      // is about 20 times larger than the precision that the integrator actually reaches.
      // (See unit tests in roofit/roofitcore/test)
      // That's why we add a fudge factor here that makes the integrator stop faster:
      if(fabs(extrapError) <= _epsRel*fabs(extrapValue) * 15.) {
        return extrapValue;
      }
      if(fabs(extrapError) <= _epsAbs) {
        return extrapValue ;
      }

    }
  }

  ostream& str = oocoutW((TObject*)0,Integration) << "RooIntegrator1D::integral: integral of " << _function->getName()
      << " over range (" << _xmin << "," << _xmax << ") did not converge after "
      << _maxSteps << " steps.\nSeries is:" << endl;
  if (_doExtrap) {
    auto result = richardsonExtrapolation<seriesAccelerationSteps>(_maxSteps, _s, &str);
    str << std::endl;
    return result.first;
  } else {
    for (unsigned int i=0; i < _s.size(); ++i) {
      str << "\n\ti=" << i << "\t" << _s[i];
    }
    str << std::endl;
    return _s[_maxSteps-1];
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Calculate the first n stages of refinement of the extended trapezoidal
/// summation rule. Start by approximating the integral with a single trapezoid,
/// and then split existing trapezoids in the middle n times.
///
/// This function takes a shortcut by only invoking the evaluation of the integrand
/// once for all points that are needed for the first n steps of the trapezoid series.
/// \return vector with summation results when splitting trapezoids in the middle n times.
std::vector<double> RooIntegrator1D::computeTrapezoids(unsigned int start, unsigned int end, double previousSum, const double* parameters, std::size_t nPar) const {
  std::vector<double> xValues;
  xValues.reserve((1 << end) - (1 << start) + 1);

  // Prepare all x values where we want to evaluate
  for (unsigned int step = start; step < end; ++step) {
    if (step == 0) {
      xValues.push_back(_xmin);
      xValues.push_back(_xmax);
      continue;
    }

    const unsigned int nTrapezoids = 1 << step; //== 2^step
    const double width = (_xmax - _xmin) / nTrapezoids;

    // Skip every second point, since they were already in the previous iteration
    for (unsigned int i=0; i < nTrapezoids/2; ++i) {
      xValues.push_back(_xmin + (2*i + 1) * width);
    }
  }

  if (xValues.empty())
    return {};



  // Now evaluate function, and compute trapezoids
  auto functionValues = evalIntegrand(xValues, parameters, nPar);
  assert(functionValues.size() == xValues.size());


  std::vector<double> results;
  results.reserve(_maxSteps); // Outside, we may be swapping into _s, so can reserve _maxSteps directly
  unsigned int nProcessed = 0;
  for (unsigned int step = start; step < end; ++step) {
    const unsigned int nTrapezoids = 1 << step; //== 2^step
    const double width = (_xmax - _xmin) / nTrapezoids;

    if (step == 0) {
      previousSum = 0.5 * width * (functionValues[0] + functionValues[1]);
      results.push_back(previousSum);
      nProcessed = 2;
      continue;
    }

    // Previous trapezoids were split in two, so start with their surface:
    ROOT::Math::KahanSum<double> accumulator(0.5 * previousSum);

    // Add the missing points.
    // Half of the points have been taken into account, so select only new ones
    for (unsigned int k = 0; k < nTrapezoids/2; ++k) {
      accumulator += width * functionValues[nProcessed++];
    }

    previousSum = accumulator;
    results.push_back(accumulator);
  }

  assert(nProcessed == functionValues.size());

  return results;
}


////////////////////////////////////////////////////////////////////////////////
/// Calculate the first n stages of refinement of the Second Euler-Maclaurin
/// summation rule which has the useful property of not evaluating the
/// integrand at either of its endpoints. However, it requires more function
/// evaluations than the trapezoidal rule. This rule can be used with
/// a suitable change of variables to estimate improper integrals.
std::vector<double> RooIntegrator1D::computeMidpoints(unsigned int start, unsigned int end, double previousSum, const double* parameters, std::size_t nPar) const {
  std::vector<double> xValues;
  xValues.reserve(static_cast<unsigned int>(std::pow(3, end-start)));

  // Prepare all x values where we want to evaluate
  for (unsigned int step = start; step < end; ++step) {
    unsigned int nBins = 1;
    for (unsigned int i=0; i < step; ++i)
      nBins *= 3; //== 3^step

    const double width = (_xmax - _xmin) / nBins;

    // Skip every point with index divisible by 3, since they were already in previous iterations
    for (unsigned int i=0; i < nBins; ++i) {
      const unsigned int elm = 2*i+1;
      if (elm % 3 == 0)
        continue;

      xValues.push_back(_xmin + elm * 0.5*width);
    }
  }

  if (xValues.empty())
    return {};



  // Now evaluate function, and compute bin integrals using midpoints
  auto functionValues = evalIntegrand(xValues, parameters, nPar);

  std::vector<double> results;
  results.reserve(_maxSteps - 1);
  unsigned int nProcessed = 0;
  for (unsigned int step = start; step < end; ++step) {
    unsigned int nBins = 1;
    for (unsigned int i=0; i < step; ++i)
      nBins *= 3; //== 3^step

    const double width = (_xmax - _xmin) / nBins;

    if (step == 0) {
      results.push_back(width * functionValues[nProcessed++]);
      previousSum = results.front();
      continue;
    }

    // Divide previous bins in three
    ROOT::Math::KahanSum<double> accumulator(previousSum/3.);

    const unsigned int previousNBins = nBins/3;
    for (unsigned int i=0; i < nBins - previousNBins; ++i) {
      accumulator += width * functionValues[nProcessed++];
    }

    results.push_back(accumulator);
    previousSum = results.back();
  }

  assert(nProcessed == functionValues.size());

  return results;
}
