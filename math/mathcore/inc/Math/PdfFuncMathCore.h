// @(#)root/mathcore:$Id$
// Authors: Andras Zsenei & Lorenzo Moneta   06/2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 , LCG ROOT MathLib Team                         *
 *                                                                    *
 *                                                                    *
 **********************************************************************/



/**

Probability density functions, cumulative distribution functions
and their inverses (quantiles) for various statistical distributions (continuous and discrete).
Whenever possible the conventions followed are those of the
CRC Concise Encyclopedia of Mathematics, Second Edition
(or <A HREF="http://mathworld.wolfram.com/">Mathworld</A>).
By convention the distributions are centered around 0, so for
example in the case of a Gaussian there is no parameter mu. The
user must calculate the shift themselves if they wish.

MathCore provides the majority of the probability density functions, of the
cumulative distributions and of the quantiles (inverses of the cumulatives).
Additional distributions are also provided by the
<A HREF="../../MathMore/html/group__StatFunc.html">MathMore</A> library.


@defgroup StatFunc Statistical functions

@ingroup  MathCore
@ingroup  MathMore

*/

#ifndef ROOT_Math_PdfFuncMathCore
#define ROOT_Math_PdfFuncMathCore

#include "Math/Math.h"
#include "Math/SpecFuncMathCore.h"

#include <limits>

namespace ROOT {
namespace Math {



  /** @defgroup PdfFunc Probability Density Functions (PDF)
   *   @ingroup StatFunc
   *  Probability density functions of various statistical distributions
   *  (continuous and discrete).
   *  The probability density function returns the probability that
   *  the variate has the value x.
   *  In statistics the PDF is also called the frequency function.
   *
   *
   */

   /** @name Probability Density Functions from MathCore
   *   Additional PDF's are provided in the MathMore library
   *   (see PDF functions from MathMore)
   */

  //@{

  /**

  Probability density function of the beta distribution.

  \f[ p(x) = \frac{\Gamma (a + b) } {\Gamma(a)\Gamma(b) } x ^{a-1} (1 - x)^{b-1} \f]

  for \f$0 \leq x \leq 1 \f$. For detailed description see
  <A HREF="http://mathworld.wolfram.com/BetaDistribution.html">
  Mathworld</A>.

  @ingroup PdfFunc

  */

  inline double beta_pdf(double x, double a, double b) {
    // Inlined to enable clad-auto-derivation for this function.

    if (x < 0 || x > 1.0) return 0;
    if (x == 0 ) {
      // need this work Windows
      if (a < 1) return  std::numeric_limits<double>::infinity();
      else if (a > 1) return  0;
      else if ( a == 1) return b; // to avoid a nan from log(0)*0
    }
    if (x == 1 ) {
      // need this work Windows
      if (b < 1) return  std::numeric_limits<double>::infinity();
      else if (b > 1) return  0;
      else if ( b == 1) return a; // to avoid a nan from log(0)*0
    }
    return std::exp( ROOT::Math::lgamma(a + b) - ROOT::Math::lgamma(a) - ROOT::Math::lgamma(b) +
                     std::log(x) * (a -1.) + ROOT::Math::log1p(-x ) * (b - 1.) );
  }



  /**

  Probability density function of the binomial distribution.

  \f[ p(k) = \frac{n!}{k! (n-k)!} p^k (1-p)^{n-k} \f]

  for \f$ 0 \leq k \leq n \f$. For detailed description see
  <A HREF="http://mathworld.wolfram.com/BinomialDistribution.html">
  Mathworld</A>.

  @ingroup PdfFunc

  */

  inline double binomial_pdf(unsigned int k, double p, unsigned int n) {
    // Inlined to enable clad-auto-derivation for this function.
    if (k > n)
      return 0.0;

    double coeff = ROOT::Math::lgamma(n+1) - ROOT::Math::lgamma(k+1) - ROOT::Math::lgamma(n-k+1);
    return std::exp(coeff + k * std::log(p) + (n - k) * ROOT::Math::log1p(-p));
  }



  /**

  Probability density function of the negative binomial distribution.

  \f[ p(k) = \frac{(k+n-1)!}{k! (n-1)!} p^{n} (1-p)^{k} \f]

  For detailed description see
  <A HREF="http://mathworld.wolfram.com/NegativeBinomialDistribution.html">
  Mathworld</A> (where \f$k \to x\f$ and \f$n \to r\f$).
  The distribution in <A HREF="http://en.wikipedia.org/wiki/Negative_binomial_distribution">
  Wikipedia</A> is defined with a \f$p\f$ corresponding to \f$1-p\f$ in this case.


  @ingroup PdfFunc

  */

  inline double negative_binomial_pdf(unsigned int k, double p, double n) {
    // Inlined to enable clad-auto-derivation for this function.

    // implement in term of gamma function

    if (n < 0)  return 0.0;
    if (p < 0 || p > 1.0) return 0.0;

    double coeff = ROOT::Math::lgamma(k+n) - ROOT::Math::lgamma(k+1.0) - ROOT::Math::lgamma(n);
    return std::exp(coeff + n * std::log(p) + double(k) * ROOT::Math::log1p(-p));

  }




  /**

  Probability density function of Breit-Wigner distribution, which is similar, just
  a different definition of the parameters, to the Cauchy distribution
  (see  #cauchy_pdf )

  \f[ p(x) = \frac{1}{\pi} \frac{\frac{1}{2} \Gamma}{x^2 + (\frac{1}{2} \Gamma)^2} \f]


  @ingroup PdfFunc

  */

  inline double breitwigner_pdf(double x, double gamma, double x0 = 0) {
    // Inlined to enable clad-auto-derivation for this function.
    double gammahalf = gamma/2.0;
    return gammahalf/(M_PI * ((x-x0)*(x-x0) + gammahalf*gammahalf));
  }




  /**

  Probability density function of the Cauchy distribution which is also
  called Lorentzian distribution.


  \f[ p(x) = \frac{1}{\pi} \frac{ b }{ (x-m)^2 + b^2} \f]

  For detailed description see
  <A HREF="http://mathworld.wolfram.com/CauchyDistribution.html">
  Mathworld</A>. It is also related to the #breitwigner_pdf which
  will call the same implementation.

  @ingroup PdfFunc

  */

  inline double cauchy_pdf(double x, double b = 1, double x0 = 0) {

    return b/(M_PI * ((x-x0)*(x-x0) + b*b));

  }




  /**

  Probability density function of the \f$\chi^2\f$ distribution with \f$r\f$
  degrees of freedom.

  \f[ p_r(x) = \frac{1}{\Gamma(r/2) 2^{r/2}} x^{r/2-1} e^{-x/2} \f]

  for \f$x \geq 0\f$. For detailed description see
  <A HREF="http://mathworld.wolfram.com/Chi-SquaredDistribution.html">
  Mathworld</A>.

  @ingroup PdfFunc

  */

  inline double chisquared_pdf(double x, double r, double x0 = 0) {
    // Inlined to enable clad-auto-derivation for this function.

    if ((x-x0) <  0) {
      return 0.0;
    }
    double a = r/2 -1.;
    // let return inf for case x  = x0 and treat special case of r = 2 otherwise will return nan
    if (x == x0 && a == 0) return 0.5;

    return std::exp ((r/2 - 1) * std::log((x-x0)/2) - (x-x0)/2 - ROOT::Math::lgamma(r/2))/2;

  }


  /**

  Crystal ball function

  See the definition at
  <A HREF="http://en.wikipedia.org/wiki/Crystal_Ball_function">
  Wikipedia</A>.

  It is not really a pdf since it is not normalized

  @ingroup PdfFunc

  */

  inline double crystalball_function(double x, double alpha, double n, double sigma, double mean = 0) {
     // Inlined to enable clad-auto-derivation for this function.

     // evaluate the crystal ball function
     if (sigma < 0.)     return 0.;
     double z = (x - mean)/sigma;
     if (alpha < 0) z = -z;
     double abs_alpha = std::abs(alpha);
     // double C = n/abs_alpha * 1./(n-1.) * std::exp(-alpha*alpha/2.);
     // double D = std::sqrt(M_PI/2.)*(1.+ROOT::Math::erf(abs_alpha/std::sqrt(2.)));
     // double N = 1./(sigma*(C+D));
     if (z  > - abs_alpha)
        return std::exp(- 0.5 * z * z);
     //double A = std::pow(n/abs_alpha,n) * std::exp(-0.5*abs_alpha*abs_alpha);
     double nDivAlpha = n/abs_alpha;
     double AA =  std::exp(-0.5*abs_alpha*abs_alpha);
     double B = nDivAlpha -abs_alpha;
     double arg = nDivAlpha/(B-z);
     return AA * std::pow(arg,n);
   }

   /**
       pdf definition of the crystal_ball which is defined only for n > 1 otherwise integral is diverging
    */
  inline double crystalball_pdf(double x, double alpha, double n, double sigma, double mean = 0) {
    // Inlined to enable clad-auto-derivation for this function.

    // evaluation of the PDF ( is defined only for n >1)
    if (sigma < 0.)     return 0.;
    if ( n <= 1) return std::numeric_limits<double>::quiet_NaN();  // pdf is not normalized for n <=1
    double abs_alpha = std::abs(alpha);
    double C = n/abs_alpha * 1./(n-1.) * std::exp(-alpha*alpha/2.);
    double D = std::sqrt(M_PI/2.)*(1.+ROOT::Math::erf(abs_alpha/std::sqrt(2.)));
    double N = 1./(sigma*(C+D));
    return N * crystalball_function(x,alpha,n,sigma,mean);
   }

  /**

  Probability density function of the exponential distribution.

  \f[ p(x) = \lambda e^{-\lambda x} \f]

  for x>0. For detailed description see
  <A HREF="http://mathworld.wolfram.com/ExponentialDistribution.html">
  Mathworld</A>.


  @ingroup PdfFunc

  */

  inline double exponential_pdf(double x, double lambda, double x0 = 0) {
    // Inlined to enable clad-auto-derivation for this function.

    if ((x-x0) < 0)
      return 0.0;
    return lambda * std::exp (-lambda * (x-x0));
  }




  /**

  Probability density function of the F-distribution.

  \f[ p_{n,m}(x) = \frac{\Gamma(\frac{n+m}{2})}{\Gamma(\frac{n}{2}) \Gamma(\frac{m}{2})} n^{n/2} m^{m/2} x^{n/2 -1} (m+nx)^{-(n+m)/2} \f]

  for x>=0. For detailed description see
  <A HREF="http://mathworld.wolfram.com/F-Distribution.html">
  Mathworld</A>.

  @ingroup PdfFunc

  */


  inline double fdistribution_pdf(double x, double n, double m, double x0 = 0) {
    // Inlined to enable clad-auto-derivation for this function.

    // function is defined only for both n and m > 0
    if (n < 0 || m < 0)
      return std::numeric_limits<double>::quiet_NaN();
    if ((x-x0) < 0)
      return 0.0;

    return std::exp((n/2) * std::log(n) + (m/2) * std::log(m) + ROOT::Math::lgamma((n+m)/2) - ROOT::Math::lgamma(n/2) - ROOT::Math::lgamma(m/2)
                    + (n/2 -1) * std::log(x-x0) - ((n+m)/2) * std::log(m +  n*(x-x0)) );

  }




  /**

  Probability density function of the gamma distribution.

  \f[ p(x) = {1 \over \Gamma(\alpha) \theta^{\alpha}} x^{\alpha-1} e^{-x/\theta} \f]

  for x>0. For detailed description see
  <A HREF="http://mathworld.wolfram.com/GammaDistribution.html">
  Mathworld</A>.

  @ingroup PdfFunc

  */

  inline double gamma_pdf(double x, double alpha, double theta, double x0 = 0) {
    // Inlined to enable clad-auto-derivation for this function.

    if ((x-x0) < 0) {
      return 0.0;
    } else if ((x-x0) == 0) {

      if (alpha == 1) {
        return 1.0/theta;
      } else {
        return 0.0;
      }

    } else if (alpha == 1) {
      return std::exp(-(x-x0)/theta)/theta;
    } else {
      return std::exp((alpha - 1) * std::log((x-x0)/theta) - (x-x0)/theta - ROOT::Math::lgamma(alpha))/theta;
    }

  }




  /**

  Probability density function of the normal (Gaussian) distribution.

  \f[ p(x) = {1 \over \sqrt{2 \pi \sigma^2}} e^{-x^2 / 2\sigma^2} \f]

  For detailed description see
  <A HREF="http://mathworld.wolfram.com/NormalDistribution.html">
  Mathworld</A>. It can also be evaluated using #normal_pdf which will
  call the same implementation.

  @ingroup PdfFunc

  */

  inline double gaussian_pdf(double x, double sigma = 1, double x0 = 0) {

    double tmp = (x-x0)/sigma;
    return (1.0/(std::sqrt(2 * M_PI) * std::fabs(sigma))) * std::exp(-tmp*tmp/2);
  }

   /**

  Probability density function of the bi-dimensional (Gaussian) distribution.

  \f[ p(x) = {1 \over 2 \pi \sigma_x \sigma_y \sqrt{1-\rho^2}} \exp (-(x^2/\sigma_x^2 + y^2/\sigma_y^2 - 2 \rho x y/(\sigma_x\sigma_y))/2(1-\rho^2)) \f]

  For detailed description see
  <A HREF="http://mathworld.wolfram.com/BivariateNormalDistribution.html">
  Mathworld</A>. It can also be evaluated using #normal_pdf which will
  call the same implementation.

 @param rho correlation , must be between -1,1

  @ingroup PdfFunc

  */

  inline double bigaussian_pdf(double x, double y, double sigmax = 1, double sigmay = 1, double rho = 0, double x0 = 0, double y0 = 0) {
    double u = (x-x0)/sigmax;
    double v = (y-y0)/sigmay;
    double c = 1. - rho*rho;
    double z = u*u - 2.*rho*u*v + v*v;
    return  1./(2 * M_PI * sigmax * sigmay * std::sqrt(c) ) * std::exp(- z / (2. * c) );
  }

   /**

   Probability density function of the Landau distribution:
  \f[ p(x) = \frac{1}{\xi} \phi (\lambda) \f]
   with
   \f[  \phi(\lambda) = \frac{1}{2 \pi i}\int_{c-i\infty}^{c+i\infty} e^{\lambda s + s \log{s}} ds\f]
   where \f$\lambda = (x-x_0)/\xi\f$. For a detailed description see
   K.S. K&ouml;lbig and B. Schorr, A program package for the Landau distribution,
   <A HREF="http://dx.doi.org/10.1016/0010-4655(84)90085-7">Computer Phys. Comm. 31 (1984) 97-111</A>
   <A HREF="http://dx.doi.org/10.1016/j.cpc.2008.03.002">[Erratum-ibid. 178 (2008) 972]</A>.
   The same algorithms as in
   <A HREF="https://cern-tex.web.cern.ch/cern-tex/shortwrupsdir/g110/top.html">
   CERNLIB</A> (DENLAN)  is used

   @param x The argument \f$x\f$
   @param xi The width parameter \f$\xi\f$
   @param x0 The location parameter \f$x_0\f$

   @ingroup PdfFunc

   */

   double landau_pdf(double x, double xi = 1, double x0 = 0);



  /**

  Probability density function of the lognormal distribution.

  \f[ p(x) = {1 \over x \sqrt{2 \pi s^2} } e^{-(\ln{x} - m)^2/2 s^2} \f]

  for x>0. For detailed description see
  <A HREF="http://mathworld.wolfram.com/LogNormalDistribution.html">
  Mathworld</A>.
  @param s scale parameter (not the sigma of the distribution which is not even defined)
  @param x0  location parameter, corresponds approximately to the most probable value. For x0 = 0, sigma = 1, the x_mpv = -0.22278

  @ingroup PdfFunc

  */

  inline double lognormal_pdf(double x, double m, double s, double x0 = 0) {
    // Inlined to enable clad-auto-derivation for this function.
    if ((x-x0) <= 0)
      return 0.0;
    double tmp = (std::log((x-x0)) - m)/s;
    return 1.0 / ((x-x0) * std::fabs(s) * std::sqrt(2 * M_PI)) * std::exp(-(tmp * tmp) /2);
  }




  /**

  Probability density function of the normal (Gaussian) distribution.

  \f[ p(x) = {1 \over \sqrt{2 \pi \sigma^2}} e^{-x^2 / 2\sigma^2} \f]

  For detailed description see
  <A HREF="http://mathworld.wolfram.com/NormalDistribution.html">
  Mathworld</A>. It can also be evaluated using #gaussian_pdf which will call the same
  implementation.

  @ingroup PdfFunc

  */

  inline double normal_pdf(double x, double sigma =1, double x0 = 0) {
    // Inlined to enable clad-auto-derivation for this function.

    double tmp = (x-x0)/sigma;
    return (1.0/(std::sqrt(2 * M_PI) * std::fabs(sigma))) * std::exp(-tmp*tmp/2);

  }


  /**

  Probability density function of the Poisson distribution.

  \f[ p(n) = \frac{\mu^n}{n!} e^{- \mu} \f]

  For detailed description see
  <A HREF="http://mathworld.wolfram.com/PoissonDistribution.html">
  Mathworld</A>.

  @ingroup PdfFunc

  */

  inline double poisson_pdf(unsigned int n, double mu) {
    // Inlined to enable clad-auto-derivation for this function.

    if (n > 0 && mu >= 0)
      return std::exp (n*std::log(mu) - ROOT::Math::lgamma(n+1) - mu);

    //  when  n = 0 and mu = 0,  1 is returned
    if (mu >= 0)
      return std::exp(-mu);

    // return a nan for mu < 0 since it does not make sense
    return std::numeric_limits<double>::quiet_NaN();
  }




  /**

  Probability density function of Student's t-distribution.

  \f[ p_{r}(x) = \frac{\Gamma(\frac{r+1}{2})}{\sqrt{r \pi}\Gamma(\frac{r}{2})} \left( 1+\frac{x^2}{r}\right)^{-(r+1)/2}  \f]

  for \f$k \geq 0\f$. For detailed description see
  <A HREF="http://mathworld.wolfram.com/Studentst-Distribution.html">
  Mathworld</A>.

  @ingroup PdfFunc

  */

  inline double tdistribution_pdf(double x, double r, double x0 = 0) {
    // Inlined to enable clad-auto-derivation for this function.

    return (std::exp (ROOT::Math::lgamma((r + 1.0)/2.0) - ROOT::Math::lgamma(r/2.0)) / std::sqrt (M_PI * r))
    * std::pow ((1.0 + (x-x0)*(x-x0)/r), -(r + 1.0)/2.0);

  }




  /**

  Probability density function of the uniform (flat) distribution.

  \f[ p(x) = {1 \over (b-a)} \f]

  if \f$a \leq x<b\f$ and 0 otherwise. For detailed description see
  <A HREF="http://mathworld.wolfram.com/UniformDistribution.html">
  Mathworld</A>.

  @ingroup PdfFunc

  */

  inline double uniform_pdf(double x, double a, double b, double x0 = 0) {
    // Inlined to enable clad-auto-derivation for this function.

    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! when a=b

    if ((x-x0) < b && (x-x0) >= a)
      return 1.0/(b - a);
    return 0.0;

  }



} // namespace Math
} // namespace ROOT



#endif // ROOT_Math_PdfFunc
