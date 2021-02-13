// Author: Stephan Hageboeck, CERN, Oct 2019
// Based on RooIpatia2 by Diego Martinez Santos, Nikhef, Diego.Martinez.Santos@cern.ch
/*****************************************************************************
 * Project: RooFit                                                           *
 * Authors:                                                                  *
 *   WV, Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu       *
 *   DK, David Kirkby,    UC Irvine,         dkirkby@uci.edu                 *
 *                                                                           *
 * Copyright (c) 2000-2019, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

/**
 * \class RooHypatia2
 *
 * RooHypatia2 is the two-sided version of the Hypatia distribution described in https://arxiv.org/abs/1312.5000.
 * \image html RooHypatia2.png
 *
 * It has a hyperbolic core of a crystal-ball-like function \f$ G \f$ and two tails:
 * \f[
 *  \mathrm{Hypatia2}(x, \mu, \sigma, \lambda, \zeta, \beta, a_l, n_l, a_r, n_r) =
 *  \begin{cases}
 *  \frac{ G(\mu - a_l \sigma, \mu, \sigma, \lambda, \zeta, \beta)                             }
 *       { \left( 1 - \frac{x}{n_l G(\ldots)/G'(\ldots) - a_l\sigma } \right)^{n_l} }
 *       & \text{if } \frac{x-\mu}{\sigma} < -a_l \\
 *  \left( (x-\mu)^2 + A^2_\lambda(\zeta)\sigma^2 \right)^{\frac{1}{2}\lambda-\frac{1}{4}} e^{\beta(x-\mu)} K_{\lambda-\frac{1}{2}}
 *      \left( \zeta \sqrt{1+\left( \frac{x-\mu}{A_\lambda(\zeta)\sigma} \right)^2 } \right) \equiv G(x, \mu, \ldots)
 *      & \text{otherwise} \\
 *  \frac{ G(\mu + a_r \sigma, \mu, \sigma, \lambda, \zeta, \beta)                               }
 *       { \left( 1 - \frac{x}{-n_r G(\ldots)/G'(\ldots) - a_r\sigma } \right)^{n_r} }
 *       & \text{if } \frac{x-\mu}{\sigma} > a_r \\
 *  \end{cases}
 * \f]
 * Here, \f$ K_\lambda \f$ are the modified Bessel functions of the second kind
 * ("irregular modified cylindrical Bessel functions" from the gsl,
 * "special Bessel functions of the third kind"),
 * and \f$ A^2_\lambda(\zeta) \f$ is a ratio of these:
 * \f[
 *   A_\lambda^{2}(\zeta) = \frac{\zeta K_\lambda(\zeta)}{K_{\lambda+1}(\zeta)}
 * \f]
 *
 * \if false
 * TODO Enable once analytic integrals work.
 * ### Analytical Integration
 * The Hypatia distribution can be integrated analytically if \f$ \beta = \zeta = 0 \f$ and
 * \f$ \lambda < 0 \f$. An analytic integral will only be used, though, if the parameters are **constant**
 * at zero, and if \f$ \lambda < 0 \f$. This can be ensured as follows:
 * ```
 * RooRealVar beta("beta", "beta", 0.); // NOT beta("beta", "beta", 0., -1., 1.) This would allow it to float.
 * RooRealVar zeta("zeta", "zeta", 0.);
 * RooRealVar lambda("lambda", "lambda", -1., -10., -0.00001);
 * ```
 * In all other cases, slower / less accurate numeric integration will be used.
 * Note that including `0.` in the value range of lambda excludes using analytic integrals.
 * \endif
 *
 * ### Concavity
 * Note that unless the parameters \f$ a_l,\ a_r \f$ are very large, the function has non-hyperbolic tails. This requires
 * \f$ G \f$ to be strictly concave, *i.e.*, peaked, as otherwise the tails would yield imaginary numbers. Choosing \f$ \lambda,
 * \beta, \zeta \f$ inappropriately will therefore lead to evaluation errors.
 *
 * Further, the original paper establishes that to keep the tails from rising,
 * \f[
 * \begin{split}
 * \beta^2 &< \alpha^2 \\
 * \Leftrightarrow \beta^2 &< \frac{\zeta^2}{\delta^2} = \frac{\zeta^2}{\sigma^2 A_{\lambda}^2(\zeta)}
 * \end{split}
 * \f]
 * needs to be satisfied, unless the fit range is very restricted, because otherwise, the function rises in the tails.
 *
 *
 * In case of evaluation errors, it is advisable to choose very large values for \f$ a_l,\ a_r \f$, tweak the parameters of the core region to
 * make it concave, and re-enable the tails. Especially \f$ \beta \f$ needs to be close to zero.
 *
 * ## Relation to RooIpatia2
 * This implementation is largely based on RooIpatia2, https://gitlab.cern.ch/lhcb/Urania/blob/master/PhysFit/P2VV/src/RooIpatia2.cxx,
 * but there are differences:
 * - At construction time, the Hypatia implementation checks if the range of parameters extends into regions where
 *   the function might be undefined.
 * - Hypatia supports I/O to ROOT files.
 * - Hypatia will support faster batched function evaluations.
 * - Hypatia might support analytical integration for the case \f$ \zeta = \beta = 0, \lambda < 1 \f$.
 *
 * Because of these differences, and to avoid name clashes for setups where RooFit is used in an environment that also has
 * RooIpatia2, class names need to be different.
 */

#include "RooHypatia2.h"
#include "RooBatchCompute.h"
#include "RooAbsReal.h"
#include "RooHelpers.h"

#include "TMath.h"
#include "Math/SpecFunc.h"
#include "ROOT/RConfig.hxx"

#include <cmath>
#include <algorithm>

///////////////////////////////////////////////////////////////////////////////////////////
/// Construct a new Hypatia2 PDF.
/// \param[in] name Name of this instance.
/// \param[in] title Title (for plotting)
/// \param[in] x The variable of this distribution
/// \param[in] lambda Shape parameter. Note that \f$ \lambda < 0 \f$ is required if \f$ \zeta = 0 \f$.
/// \param[in] zeta Shape parameter (\f$ \zeta >= 0 \f$).
/// \param[in] beta Asymmetry parameter \f$ \beta \f$. Symmetric case is \f$ \beta = 0 \f$,
/// choose values close to zero.
/// \param[in] sigma Width parameter. If \f$ \beta = 0, \ \sigma \f$ is the RMS width.
/// \param[in] mu Location parameter. Shifts the distribution left/right.
/// \param[in] a Start of the left tail (\f$ a \geq 0 \f$, to the left of the peak). Note that when setting \f$ a = \sigma = 1 \f$,
/// the tail region is to the left of \f$ x = \mu - 1 \f$, so a should be positive.
/// \param[in] n Shape parameter of left tail (\f$ n \ge 0 \f$). With \f$ n = 0 \f$, the function is constant.
/// \param[in] a2 Start of right tail.
/// \param[in] n2 Shape parameter of right tail (\f$ n2 \ge 0 \f$). With \f$ n2 = 0 \f$, the function is constant.
RooHypatia2::RooHypatia2(const char *name, const char *title, RooAbsReal& x, RooAbsReal& lambda,
    RooAbsReal& zeta, RooAbsReal& beta, RooAbsReal& argSigma, RooAbsReal& mu, RooAbsReal& a,
    RooAbsReal& n, RooAbsReal& a2, RooAbsReal& n2) :
                RooAbsPdf(name, title),
                _x("x", "x", this, x),
                _lambda("lambda", "Lambda", this, lambda),
                _zeta("zeta", "zeta", this, zeta),
                _beta("beta", "Asymmetry parameter beta", this, beta),
                _sigma("sigma", "Width parameter sigma", this, argSigma),
                _mu("mu", "Location parameter mu", this, mu),
                _a("a", "Left tail location a", this, a),
                _n("n", "Left tail parameter n", this, n),
                _a2("a2", "Right tail location a2", this, a2),
                _n2("n2", "Right tail parameter n2", this, n2)
{
  RooHelpers::checkRangeOfParameters(this, {&argSigma}, 0.);
  RooHelpers::checkRangeOfParameters(this, {&zeta, &n, &n2, &a, &a2}, 0., std::numeric_limits<double>::max(), true);
  if (zeta.getVal() == 0. && zeta.isConstant()) {
    RooHelpers::checkRangeOfParameters(this, {&lambda}, -std::numeric_limits<double>::max(), 0., false,
        std::string("Lambda needs to be negative when ") + _zeta.GetName() + " is zero.");
  }

#ifndef R__HAS_MATHMORE
  throw std::logic_error("RooHypatia2 needs ROOT with mathmore enabled to access gsl functions.");
#endif
}


///////////////////////////////////////////////////////////////////////////////////////////
/// Copy a new Hypatia2 PDF.
/// \param[in] other Original to copy from.
/// \param[in] name Optional new name.
RooHypatia2::RooHypatia2(const RooHypatia2& other, const char* name) :
                   RooAbsPdf(other, name),
                   _x("x", this, other._x),
                   _lambda("lambda", this, other._lambda),
                   _zeta("zeta", this, other._zeta),
                   _beta("beta", this, other._beta),
                   _sigma("sigma", this, other._sigma),
                   _mu("mu", this, other._mu),
                   _a("a", this, other._a),
                   _n("n", this, other._n),
                   _a2("a2", this, other._a2),
                   _n2("n2", this, other._n2)
{
#ifndef R__HAS_MATHMORE
  throw std::logic_error("RooHypatia2 needs ROOT with mathmore enabled to access gsl functions.");
#endif
}

namespace {
const double sq2pi_inv = 1./std::sqrt(TMath::TwoPi());
const double logsq2pi = std::log(std::sqrt(TMath::TwoPi()));
const double ln2 = std::log(2.);

double low_x_BK(double nu, double x){
  return TMath::Gamma(nu)*std::pow(2., nu-1.)*std::pow(x, -nu);
}

double low_x_LnBK(double nu, double x){
  return std::log(TMath::Gamma(nu)) + ln2*(nu-1.) - std::log(x) * nu;
}

double besselK(double ni, double x) {
  const double nu = std::fabs(ni);
  if ((x < 1.e-06 && nu > 0.) ||
      (x < 1.e-04 && nu > 0. && nu < 55.) ||
      (x < 0.1 && nu >= 55.) )
    return low_x_BK(nu, x);

#ifdef R__HAS_MATHMORE
  return ROOT::Math::cyl_bessel_k(nu, x);
#else
  return std::numeric_limits<double>::signaling_NaN();
#endif

}

double LnBesselK(double ni, double x) {
  const double nu = std::fabs(ni);
  if ((x < 1.e-06 && nu > 0.) ||
      (x < 1.e-04 && nu > 0. && nu < 55.) ||
      (x < 0.1 && nu >= 55.) )
    return low_x_LnBK(nu, x);

#ifdef R__HAS_MATHMORE
  return std::log(ROOT::Math::cyl_bessel_k(nu, x));
#else
  return std::numeric_limits<double>::signaling_NaN();
#endif
}


double LogEval(double d, double l, double alpha, double beta, double delta) {
  const double gamma = alpha;//std::sqrt(alpha*alpha-beta*beta);
  const double dg = delta*gamma;
  const double thing = delta*delta + d*d;
  const double logno = l*std::log(gamma/delta) - logsq2pi - LnBesselK(l, dg);

  return std::exp(logno + beta*d
      + (0.5-l)*(std::log(alpha)-0.5*std::log(thing))
      + LnBesselK(l-0.5, alpha*std::sqrt(thing)) );// + std::log(std::fabs(beta)+0.0001) );

}


double diff_eval(double d, double l, double alpha, double beta, double delta){
  const double gamma = alpha;
  const double dg = delta*gamma;

  const double thing = delta*delta + d*d;
  const double sqrthing = std::sqrt(thing);
  const double alphasq = alpha*sqrthing;
  const double no = std::pow(gamma/delta, l)/besselK(l, dg)*sq2pi_inv;
  const double ns1 = 0.5-l;

  return no * std::pow(alpha, ns1) * std::pow(thing, l/2. - 1.25)
  * ( -d * alphasq * (besselK(l - 1.5, alphasq)
      + besselK(l + 0.5, alphasq))
      + (2.*(beta*thing + d*l) - d) * besselK(ns1, alphasq) )
      * std::exp(beta*d) * 0.5;
}

/*
double Gauss2F1(double a, double b, double c, double x){
  if (fabs(x) <= 1.) {
    return ROOT::Math::hyperg(a, b, c, x);
  } else {
    return ROOT::Math::hyperg(c-a, b, c, 1-1/(1-x))/std::pow(1-x, b);
  }
}

double stIntegral(double d1, double delta, double l){
  return d1 * Gauss2F1(0.5, 0.5-l, 3./2, -d1*d1/(delta*delta));
}
*/
}

double RooHypatia2::evaluate() const
{
  const double d = _x-_mu;
  const double cons0 = std::sqrt(_zeta);
  const double asigma = _a*_sigma;
  const double a2sigma = _a2*_sigma;
  const double beta = _beta;
  double out = 0.;

  if (_zeta > 0.) {
    // careful if zeta -> 0. You can implement a function for the ratio,
    // but careful again that |nu + 1 | != |nu| + 1 so you have to deal with the signs
    const double phi = besselK(_lambda+1., _zeta) / besselK(_lambda, _zeta);
    const double cons1 = _sigma/std::sqrt(phi);
    const double alpha  = cons0/cons1;
    const double delta = cons0*cons1;

    if (d < -asigma){
      const double k1 = LogEval(-asigma, _lambda, alpha, beta, delta);
      const double k2 = diff_eval(-asigma, _lambda, alpha, beta, delta);
      const double B = -asigma + _n*k1/k2;
      const double A = k1 * std::pow(B+asigma, _n);

      out = A * std::pow(B-d, -_n);
    }
    else if (d>a2sigma) {
      const double k1 = LogEval(a2sigma, _lambda, alpha, beta, delta);
      const double k2 = diff_eval(a2sigma, _lambda, alpha, beta, delta);
      const double B = -a2sigma - _n2*k1/k2;
      const double A = k1 * std::pow(B+a2sigma, _n2);

      out = A * std::pow(B+d, -_n2);
    }
    else {
      out = LogEval(d, _lambda, alpha, beta, delta);
    }
  }
  else if (_zeta < 0.) {
    coutE(Eval) << "The parameter " << _zeta.GetName() << " of the RooHypatia2 " << GetName() << " cannot be < 0." << std::endl;
  }
  else if (_lambda < 0.) {
    const double delta = _sigma;

    if (d < -asigma ) {
      const double cons1 = std::exp(-beta*asigma);
      const double phi = 1. + _a * _a;
      const double k1 = cons1 * std::pow(phi, _lambda-0.5);
      const double k2 = beta*k1 - cons1*(_lambda-0.5) * std::pow(phi, _lambda-1.5) * 2.*_a/delta;
      const double B = -asigma + _n*k1/k2;
      const double A = k1*std::pow(B+asigma, _n);

      out = A*std::pow(B-d, -_n);
    }
    else if (d > a2sigma) {
      const double cons1 = std::exp(beta*a2sigma);
      const double phi = 1. + _a2*_a2;
      const double k1 = cons1 * std::pow(phi, _lambda-0.5);
      const double k2 = beta*k1 + cons1*(_lambda-0.5) * std::pow(phi, _lambda-1.5) * 2.*_a2/delta;
      const double B = -a2sigma - _n2*k1/k2;
      const double A = k1*std::pow(B+a2sigma, _n2);

      out =  A*std::pow(B+d, -_n2);
    }
    else {
      out = std::exp(beta*d) * std::pow(1. + d*d/(delta*delta), _lambda - 0.5);
    }
  }
  else {
    coutE(Eval) << "zeta = 0 only supported for lambda < 0. lambda = " << double(_lambda) << std::endl;
  }

  return out;
}

namespace {
//////////////////////////////////////////////////////////////////////////////////////////
/// The generic compute function that recalculates everything for every loop iteration.
/// This is only needed in the rare case where a parameter is used as an observable.
template<typename Tx, typename Tl, typename Tz, typename Tb, typename Ts, typename Tm, typename Ta, typename Tn,
typename Ta2, typename Tn2>
void compute(RooSpan<double> output, Tx x, Tl lambda, Tz zeta, Tb beta, Ts sigma, Tm mu, Ta a, Tn n, Ta2 a2, Tn2 n2) {
  const auto N = output.size();
  const bool zetaIsAlwaysLargerZero = !zeta.isBatch() && zeta[0] > 0.;
  const bool zetaIsAlwaysZero = !zeta.isBatch() && zeta[0] == 0.;

  for (std::size_t i = 0; i < N; ++i) {

    const double d = x[i] - mu[i];
    const double cons0 = std::sqrt(zeta[i]);
    const double asigma = a[i]*sigma[i];
    const double a2sigma = a2[i]*sigma[i];
//    const double beta = beta[i];

    if (zetaIsAlwaysLargerZero || zeta[i] > 0.) {
      const double phi = besselK(lambda[i]+1., zeta[i]) / besselK(lambda[i], zeta[i]);
      const double cons1 = sigma[i]/std::sqrt(phi);
      const double alpha  = cons0/cons1;
      const double delta = cons0*cons1;

      if (d < -asigma){
        const double k1 = LogEval(-asigma, lambda[i], alpha, beta[i], delta);
        const double k2 = diff_eval(-asigma, lambda[i], alpha, beta[i], delta);
        const double B = -asigma + n[i]*k1/k2;
        const double A = k1 * std::pow(B+asigma, n[i]);

        output[i] = A * std::pow(B-d, -n[i]);
      }
      else if (d>a2sigma) {
        const double k1 = LogEval(a2sigma, lambda[i], alpha, beta[i], delta);
        const double k2 = diff_eval(a2sigma, lambda[i], alpha, beta[i], delta);
        const double B = -a2sigma - n2[i]*k1/k2;
        const double A = k1 * std::pow(B+a2sigma, n2[i]);

        output[i] = A * std::pow(B+d, -n2[i]);
      }
      else {
        output[i] = LogEval(d, lambda[i], alpha, beta[i], delta);
      }
    }

    if ((!zetaIsAlwaysLargerZero && zetaIsAlwaysZero) || (zeta[i] == 0. && lambda[i] < 0.)) {
      const double delta = sigma[i];

      if (d < -asigma ) {
        const double cons1 = std::exp(-beta[i]*asigma);
        const double phi = 1. + a[i] * a[i];
        const double k1 = cons1 * std::pow(phi, lambda[i]-0.5);
        const double k2 = beta[i]*k1 - cons1*(lambda[i]-0.5) * std::pow(phi, lambda[i]-1.5) * 2.*a[i]/delta;
        const double B = -asigma + n[i]*k1/k2;
        const double A = k1*std::pow(B+asigma, n[i]);

        output[i] = A*std::pow(B-d, -n[i]);
      }
      else if (d > a2sigma) {
        const double cons1 = std::exp(beta[i]*a2sigma);
        const double phi = 1. + a2[i]*a2[i];
        const double k1 = cons1 * std::pow(phi, lambda[i]-0.5);
        const double k2 = beta[i]*k1 + cons1*(lambda[i]-0.5) * std::pow(phi, lambda[i]-1.5) * 2.*a2[i]/delta;
        const double B = -a2sigma - n2[i]*k1/k2;
        const double A = k1*std::pow(B+a2sigma, n2[i]);

        output[i] =  A*std::pow(B+d, -n2[i]);
      }
      else {
        output[i] = std::exp(beta[i]*d) * std::pow(1. + d*d/(delta*delta), lambda[i] - 0.5);
      }
    }
  }
}

template<bool right>
std::pair<double, double> computeAB_zetaNZero(double asigma, double lambda, double alpha, double beta, double delta, double n) {
  const double k1 = LogEval(right ? asigma : -asigma, lambda, alpha, beta, delta);
  const double k2 = diff_eval(right ? asigma : -asigma, lambda, alpha, beta, delta);
  const double B = -asigma + (right ? -1 : 1.) * n*k1/k2;
  const double A = k1 * std::pow(B+asigma, n);

  return {A, B};
}

template<bool right>
std::pair<double, double> computeAB_zetaZero(double beta, double asigma, double a, double lambda, double delta, double n) {
  const double cons1 = std::exp((right ? 1. : -1.) * beta*asigma);
  const double phi = 1. + a * a;
  const double k1 = cons1 * std::pow(phi, lambda-0.5);
  const double k2 = beta*k1 + (right ? 1. : -1) * cons1*(lambda-0.5) * std::pow(phi, lambda-1.5) * 2.*a/delta;
  const double B = -asigma + (right ? -1. : 1.) * n*k1/k2;
  const double A = k1*std::pow(B+asigma, n);

  return {A, B};
}

using RooBatchCompute::BracketAdapter;
//////////////////////////////////////////////////////////////////////////////////////////
/// A specialised compute function where x is an observable, and all parameters are used as
/// parameters. Since many things can be calculated outside of the loop, it is faster.
void compute(RooSpan<double> output, RooSpan<const double> x,
    BracketAdapter<double> lambda, BracketAdapter<double> zeta, BracketAdapter<double> beta,
    BracketAdapter<double> sigma, BracketAdapter<double> mu,
    BracketAdapter<double> a, BracketAdapter<double> n, BracketAdapter<double> a2, BracketAdapter<double> n2) {
  const auto N = output.size();

  const double cons0 = std::sqrt(zeta);
  const double asigma = a*sigma;
  const double a2sigma = a2*sigma;

  if (zeta > 0.) {
    const double phi = besselK(lambda+1., zeta) / besselK(lambda, zeta);
    const double cons1 = sigma/std::sqrt(phi);
    const double alpha  = cons0/cons1;
    const double delta = cons0*cons1;

    const auto AB_l = computeAB_zetaNZero<false>(asigma, lambda, alpha, beta, delta, n);
    const auto AB_r = computeAB_zetaNZero<true>(a2sigma, lambda, alpha, beta, delta, n2);

    for (std::size_t i = 0; i < N; ++i) {
      const double d = x[i] - mu[i];

      if (d < -asigma){
        output[i] = AB_l.first * std::pow(AB_l.second - d, -n);
      }
      else if (d>a2sigma) {
        output[i] = AB_r.first * std::pow(AB_r.second + d, -n2);
      }
      else {
        output[i] = LogEval(d, lambda, alpha, beta, delta);
      }
    }
  } else if (zeta == 0. && lambda < 0.) {
    const double delta = sigma;

    const auto AB_l = computeAB_zetaZero<false>(beta, asigma,  a,  lambda, delta, n);
    const auto AB_r = computeAB_zetaZero<true>(beta, a2sigma, a2, lambda, delta, n2);

    for (std::size_t i = 0; i < N; ++i) {
      const double d = x[i] - mu[i];

      if (d < -asigma ) {
        output[i] = AB_l.first*std::pow(AB_l.second - d, -n);
      }
      else if (d > a2sigma) {
        output[i] = AB_r.first * std::pow(AB_r.second + d, -n2);
      }
      else {
        output[i] = std::exp(beta*d) * std::pow(1. + d*d/(delta*delta), lambda - 0.5);
      }
    }
  }
}

}

RooSpan<double> RooHypatia2::evaluateSpan(RooBatchCompute::RunContext& evalData, const RooArgSet* normSet) const {
  using namespace RooBatchCompute;

  auto x = _x->getValues(evalData, normSet);
  auto lambda = _lambda->getValues(evalData, normSet);
  auto zeta = _zeta->getValues(evalData, normSet);
  auto beta = _beta->getValues(evalData, normSet);
  auto sig = _sigma->getValues(evalData, normSet);
  auto mu = _mu->getValues(evalData, normSet);
  auto a = _a->getValues(evalData, normSet);
  auto n = _n->getValues(evalData, normSet);
  auto a2 = _a2->getValues(evalData, normSet);
  auto n2 = _n2->getValues(evalData, normSet);

  size_t paramSizeSum=0, batchSize = x.size();
  for (const auto& i:{lambda, zeta, beta, sig, mu, a, n, a2, n2}) {
    paramSizeSum += i.size();
    batchSize = std::max(batchSize, i.size());
  }
  RooSpan<double> output = evalData.makeBatch(this, batchSize);

  // Run high performance compute if only x has multiple values
  if (x.size()>1 && paramSizeSum==9) {
    compute(output, x,
        BracketAdapter<double>(_lambda), BracketAdapter<double>(_zeta),
        BracketAdapter<double>(_beta), BracketAdapter<double>(_sigma), BracketAdapter<double>(_mu),
        BracketAdapter<double>(_a), BracketAdapter<double>(_n),
        BracketAdapter<double>(_a2), BracketAdapter<double>(_n2));
  } else {
    compute(output, BracketAdapterWithMask(_x, x),
        BracketAdapterWithMask(_lambda, lambda), BracketAdapterWithMask(_zeta, zeta),
        BracketAdapterWithMask(_beta, beta), BracketAdapterWithMask(_sigma, sig),
        BracketAdapterWithMask(_mu, mu),
        BracketAdapterWithMask(_a, a), BracketAdapterWithMask(_n, n),
        BracketAdapterWithMask(_a2, a2), BracketAdapterWithMask(_n2, n2));
  }
  return output;
}


/* Analytical integrals need testing.

Int_t RooHypatia2::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char*) const
{
  if (matchArgs(allVars, analVars, _x)
      && _beta == 0. && _beta.arg().isConstant()
      && _zeta == 0. && _zeta.arg().isConstant()
      && _lambda.max() < 0.) return 1;
  return 0 ;
}



double RooHypatia2::analyticalIntegral(Int_t code, const char* rangeName) const
{
  if (_beta != 0. || _zeta != 0. || _lambda >= 0) {
    auto& logstream = coutF(Integration) << "When the PDF " << GetName()
        << " was constructed, beta,zeta were 0, lambda<0 and all three were constant.\n"
        << "This allowed for analytic integration, but now, numeric integration would be required.\n"
        << "These parameters must either be kept constant, or be floating from the beginning. "
        << " Terminating fit ..."
        << std::endl;
    RooArgSet vars;
    vars.add(_beta.arg());
    vars.add(_zeta.arg());
    vars.add(_lambda.arg());
    vars.printStream(logstream, vars.defaultPrintContents(nullptr), RooPrintable::kVerbose);
    throw std::runtime_error("RooHypatia2 cannot be integrated analytically since parameters changed.");
  }

  // The formulae to follow still use beta and zeta to facilitate comparisons with the
  // evaluate function. Since beta and zeta are zero, all relevant terms will be disabled
  // by defining these two constexpr:
  constexpr double beta = 0.;
  constexpr double cons1 = 1.;

  if (code != 1) {
    throw std::logic_error("Trying to compute analytic integral with unknown configuration.");
  }

  const double asigma = _a * _sigma;
  const double a2sigma = _a2 * _sigma;
  const double d0 = _x.min(rangeName) - _mu;
  const double d1 = _x.max(rangeName) - _mu;


  double delta;
  if (_lambda <= -1.0) {
    delta = _sigma * sqrt(-2. + -2.*_lambda);
  }
  else {
    delta = _sigma;
  }
  const double deltaSq = delta*delta;

  if ((d0 > -asigma) && (d1 < a2sigma)) {
    return stIntegral(d1, delta, _lambda) - stIntegral(d0, delta, _lambda);
  }

  if (d0 > a2sigma) {
    const double phi = 1. + a2sigma*a2sigma/deltaSq;
    const double k1 = cons1*std::pow(phi,_lambda-0.5);
    const double k2 = beta*k1+ cons1*(_lambda-0.5)*std::pow(phi,_lambda-1.5)*2.*a2sigma/deltaSq;
    const double B = -a2sigma - _n2*k1/k2;
    const double A = k1*std::pow(B+a2sigma,_n2);
    return A*(std::pow(B+d1,1-_n2)/(1-_n2) -std::pow(B+d0,1-_n2)/(1-_n2) ) ;

  }

  if (d1 < -asigma) {
    const double phi = 1. + asigma*asigma/deltaSq;
    const double k1 = cons1*std::pow(phi,_lambda-0.5);
    const double k2 = beta*k1- cons1*(_lambda-0.5)*std::pow(phi,_lambda-1.5)*2*asigma/deltaSq;
    const double B = -asigma + _n*k1/k2;
    const double A = k1*std::pow(B+asigma,_n);
    const double I0 = A*std::pow(B-d0,1-_n)/(_n-1);
    const double I1 = A*std::pow(B-d1,1-_n)/(_n-1);
    return I1 - I0;
  }


  double I0;
  double I1a = 0;
  double I1b = 0;
  if (d0 <-asigma) {
    const double phi = 1. + asigma*asigma/deltaSq;
    const double k1 = cons1*std::pow(phi,_lambda-0.5);
    const double k2 = beta*k1- cons1*(_lambda-0.5)*std::pow(phi,_lambda-1.5)*2*asigma/deltaSq;
    const double B = -asigma + _n*k1/k2;
    const double A = k1*std::pow(B+asigma,_n);
    I0 = A*std::pow(B-d0,1-_n)/(_n-1);
    I1a = A*std::pow(B+asigma,1-_n)/(_n-1) - stIntegral(-asigma, delta, _lambda);
  }
  else {
    I0 = stIntegral(d0, delta, _lambda);
  }

  if (d1 > a2sigma) {
    const double phi = 1. + a2sigma*a2sigma/deltaSq;
    const double k1 = cons1*std::pow(phi,_lambda-0.5);
    const double k2 = beta*k1+ cons1*(_lambda-0.5)*std::pow(phi,_lambda-1.5)*2.*a2sigma/deltaSq;
    const double B = -a2sigma - _n2*k1/k2;
    const double A = k1*std::pow(B+a2sigma,_n2);
    I1b = A*(std::pow(B+d1,1-_n2)/(1-_n2) -std::pow(B+a2sigma,1-_n2)/(1-_n2) ) - stIntegral(d1, delta, _lambda) +  stIntegral(a2sigma,delta, _lambda) ;
  }

  const double I1 = stIntegral(d1, delta, _lambda) + I1a + I1b;
  return I1 - I0;
}

*/
