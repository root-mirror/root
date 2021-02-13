// @(#)root/mathcore:$Id$
// Authors: L. Moneta, A. Zsenei   06/2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 , LCG ROOT MathLib Team                         *
 *                                                                    *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Math_ProbFuncMathCore
#define ROOT_Math_ProbFuncMathCore


namespace ROOT {
namespace Math {


   /** @defgroup ProbFunc Cumulative Distribution Functions (CDF)

   @ingroup StatFunc

   *  Cumulative distribution functions of various distributions.
   *  The functions with the extension <em>_cdf</em> calculate the
   *  lower tail integral of the probability density function
   *
   *  \f[ D(x) = \int_{-\infty}^{x} p(x') dx' \f]
   *
   *  while those with the <em>_cdf_c</em> extension calculate the complement of
   *  cumulative distribution function, called in statistics the survival
   *  function.
   *  It corresponds to the upper tail integral of the
   *  probability density function
   *
   *  \f[ D(x) = \int_{x}^{+\infty} p(x') dx' \f]
   *
   *
   * <strong>NOTE:</strong> In the old releases (< 5.14) the <em>_cdf</em> functions were called
   * <em>_quant</em> and the <em>_cdf_c</em> functions were called
   * <em>_prob</em>.
   * These names are currently kept for backward compatibility, but
   * their usage is deprecated.
   *
   *  These functions are defined in the header file <em>Math/ProbFunc.h</em> or in the global one
   *  including all statistical functions <em>Math/DistFunc.h</em>
   *
   */



   /**

   Complement of the cumulative distribution function of the beta distribution.
   Upper tail of the integral of the #beta_pdf

   @ingroup ProbFunc

   */

   double beta_cdf_c(double x, double a, double b);



   /**

   Cumulative distribution function of the beta distribution
   Upper tail of the integral of the #beta_pdf

   @ingroup ProbFunc

   */

   double beta_cdf(double x, double a, double b);




   /**

   Complement of the cumulative distribution function (upper tail) of the Breit_Wigner
   distribution and it is similar (just a different parameter definition) to the
   Cauchy distribution (see #cauchy_cdf_c )

   \f[ D(x) = \int_{x}^{+\infty} \frac{1}{\pi} \frac{\frac{1}{2} \Gamma}{x'^2 + (\frac{1}{2} \Gamma)^2} dx' \f]


   @ingroup ProbFunc

   */
   double breitwigner_cdf_c(double x, double gamma, double x0 = 0);


   /**

   Cumulative distribution function (lower tail) of the Breit_Wigner
   distribution and it is similar (just a different parameter definition) to the
   Cauchy distribution (see #cauchy_cdf )

   \f[ D(x) = \int_{-\infty}^{x} \frac{1}{\pi} \frac{b}{x'^2 + (\frac{1}{2} \Gamma)^2} dx' \f]


   @ingroup ProbFunc

   */
   double breitwigner_cdf(double x, double gamma, double x0 = 0);



   /**

   Complement of the cumulative distribution function (upper tail) of the
   Cauchy distribution which is also Lorentzian distribution.
   It is similar (just a different parameter definition) to the
   Breit_Wigner distribution (see #breitwigner_cdf_c )

   \f[ D(x) = \int_{x}^{+\infty} \frac{1}{\pi} \frac{ b }{ (x'-m)^2 + b^2} dx' \f]

   For detailed description see
   <A HREF="http://mathworld.wolfram.com/CauchyDistribution.html">
   Mathworld</A>.

   @ingroup ProbFunc

   */
   double cauchy_cdf_c(double x, double b, double x0 = 0);




   /**

   Cumulative distribution function (lower tail) of the
   Cauchy distribution which is also Lorentzian distribution.
   It is similar (just a different parameter definition) to the
   Breit_Wigner distribution (see #breitwigner_cdf )

   \f[ D(x) = \int_{-\infty}^{x} \frac{1}{\pi} \frac{ b }{ (x'-m)^2 + b^2} dx' \f]

   For detailed description see
   <A HREF="http://mathworld.wolfram.com/CauchyDistribution.html">
   Mathworld</A>.


   @ingroup ProbFunc

   */
   double cauchy_cdf(double x, double b, double x0 = 0);




   /**

   Complement of the cumulative distribution function of the \f$\chi^2\f$ distribution
   with \f$r\f$ degrees of freedom (upper tail).

   \f[ D_{r}(x) = \int_{x}^{+\infty} \frac{1}{\Gamma(r/2) 2^{r/2}} x'^{r/2-1} e^{-x'/2} dx' \f]

   For detailed description see
   <A HREF="http://mathworld.wolfram.com/Chi-SquaredDistribution.html">
   Mathworld</A>. It is implemented using the incomplete gamma function, ROOT::Math::inc_gamma_c,
   from <A HREF="http://www.netlib.org/cephes">Cephes</A>

   @ingroup ProbFunc

   */

   double chisquared_cdf_c(double x, double r, double x0 = 0);



   /**

   Cumulative distribution function of the \f$\chi^2\f$ distribution
   with \f$r\f$ degrees of freedom (lower tail).

   \f[ D_{r}(x) = \int_{-\infty}^{x} \frac{1}{\Gamma(r/2) 2^{r/2}} x'^{r/2-1} e^{-x'/2} dx' \f]

   For detailed description see
   <A HREF="http://mathworld.wolfram.com/Chi-SquaredDistribution.html">
   Mathworld</A>.   It is implemented using the incomplete gamma function, ROOT::Math::inc_gamma_c,
   from <A HREF="http://www.netlib.org/cephes">Cephes</A>

   @ingroup ProbFunc

   */

   double chisquared_cdf(double x, double r, double x0 = 0);


   /**

      Cumulative distribution for the Crystal Ball distribution function

      See the definition of the Crystal Ball function at
      <A HREF="http://en.wikipedia.org/wiki/Crystal_Ball_function">
      Wikipedia</A>.

      The distribution is defined only for n > 1 when the integral converges

      @ingroup ProbFunc

   */
   double crystalball_cdf(double x, double alpha, double n, double sigma, double x0 = 0);

   /**

      Complement of the Cumulative distribution for the Crystal Ball distribution

      See the definition of the Crystal Ball function at
      <A HREF="http://en.wikipedia.org/wiki/Crystal_Ball_function">
      Wikipedia</A>.

      The distribution is defined only for n > 1 when the integral converges

      @ingroup ProbFunc

   */
   double crystalball_cdf_c(double x, double alpha, double n, double sigma, double x0 = 0);

   /**
      Integral of the not-normalized Crystal Ball function

      See the definition of the Crystal Ball function at
      <A HREF="http://en.wikipedia.org/wiki/Crystal_Ball_function">
      Wikipedia</A>.

      see ROOT::Math::crystalball_function for the function evaluation.

      @ingroup ProbFunc

   */
   double crystalball_integral(double x, double alpha, double n, double sigma, double x0 = 0);

   /**

   Complement of the cumulative distribution function of the exponential distribution
   (upper tail).

   \f[ D(x) = \int_{x}^{+\infty} \lambda e^{-\lambda x'} dx' \f]

   For detailed description see
   <A HREF="http://mathworld.wolfram.com/ExponentialDistribution.html">
   Mathworld</A>.

   @ingroup ProbFunc

   */

   double exponential_cdf_c(double x, double lambda, double x0 = 0);



   /**

   Cumulative distribution function of the exponential distribution
   (lower tail).

   \f[ D(x) = \int_{-\infty}^{x} \lambda e^{-\lambda x'} dx' \f]

   For detailed description see
   <A HREF="http://mathworld.wolfram.com/ExponentialDistribution.html">
   Mathworld</A>.

   @ingroup ProbFunc

   */


   double exponential_cdf(double x, double lambda, double x0 = 0);



   /**

   Complement of the cumulative distribution function of the F-distribution
   (upper tail).

   \f[ D_{n,m}(x) = \int_{x}^{+\infty} \frac{\Gamma(\frac{n+m}{2})}{\Gamma(\frac{n}{2}) \Gamma(\frac{m}{2})} n^{n/2} m^{m/2} x'^{n/2 -1} (m+nx')^{-(n+m)/2} dx' \f]

   For detailed description see
   <A HREF="http://mathworld.wolfram.com/F-Distribution.html">
   Mathworld</A>. It is implemented using the incomplete beta function, ROOT::Math::inc_beta,
   from <A HREF="http://www.netlib.org/cephes">Cephes</A>

   @ingroup ProbFunc

   */

   double fdistribution_cdf_c(double x, double n, double m, double x0 = 0);




   /**

   Cumulative distribution function of the F-distribution
   (lower tail).

   \f[ D_{n,m}(x) = \int_{-\infty}^{x} \frac{\Gamma(\frac{n+m}{2})}{\Gamma(\frac{n}{2}) \Gamma(\frac{m}{2})} n^{n/2} m^{m/2} x'^{n/2 -1} (m+nx')^{-(n+m)/2} dx' \f]

   For detailed description see
   <A HREF="http://mathworld.wolfram.com/F-Distribution.html">
   Mathworld</A>. It is implemented using the incomplete beta function, ROOT::Math::inc_beta,
   from <A HREF="http://www.netlib.org/cephes">Cephes</A>

   @ingroup ProbFunc

   */

   double fdistribution_cdf(double x, double n, double m, double x0 = 0);



   /**

   Complement of the cumulative distribution function of the gamma distribution
   (upper tail).

   \f[ D(x) = \int_{x}^{+\infty} {1 \over \Gamma(\alpha) \theta^{\alpha}} x'^{\alpha-1} e^{-x'/\theta} dx' \f]

   For detailed description see
   <A HREF="http://mathworld.wolfram.com/GammaDistribution.html">
   Mathworld</A>. It is implemented using the incomplete gamma function, ROOT::Math::inc_gamma,
   from <A HREF="http://www.netlib.org/cephes">Cephes</A>

   @ingroup ProbFunc

   */

   double gamma_cdf_c(double x, double alpha, double theta, double x0 = 0);




   /**

   Cumulative distribution function of the gamma distribution
   (lower tail).

   \f[ D(x) = \int_{-\infty}^{x} {1 \over \Gamma(\alpha) \theta^{\alpha}} x'^{\alpha-1} e^{-x'/\theta} dx' \f]

   For detailed description see
   <A HREF="http://mathworld.wolfram.com/GammaDistribution.html">
   Mathworld</A>. It is implemented using the incomplete gamma function, ROOT::Math::inc_gamma,
   from <A HREF="http://www.netlib.org/cephes">Cephes</A>

   @ingroup ProbFunc

   */

   double gamma_cdf(double x, double alpha, double theta, double x0 = 0);



  /**

   Cumulative distribution function of the Landau
   distribution (lower tail).

   \f[ D(x) = \int_{-\infty}^{x} p(x) dx  \f]

   where \f$p(x)\f$ is the Landau probability density function :
  \f[ p(x) = \frac{1}{\xi} \phi (\lambda) \f]
   with
   \f[  \phi(\lambda) = \frac{1}{2 \pi i}\int_{c-i\infty}^{c+i\infty} e^{\lambda s + s \log{s}} ds\f]
   with \f$\lambda = (x-x_0)/\xi\f$. For a detailed description see
   K.S. K&ouml;lbig and B. Schorr, A program package for the Landau distribution,
   <A HREF="http://dx.doi.org/10.1016/0010-4655(84)90085-7">Computer Phys. Comm. 31 (1984) 97-111</A>
   <A HREF="http://dx.doi.org/10.1016/j.cpc.2008.03.002">[Erratum-ibid. 178 (2008) 972]</A>.
   The same algorithms as in
   <A HREF="https://cern-tex.web.cern.ch/cern-tex/shortwrupsdir/g110/top.html">
   CERNLIB</A> (DISLAN) is used.

   @param x The argument \f$x\f$
   @param xi The width parameter \f$\xi\f$
   @param x0 The location parameter \f$x_0\f$

   @ingroup ProbFunc

   */

   double landau_cdf(double x, double xi = 1, double x0 = 0);

  /**

     Complement of the distribution function of the Landau
     distribution (upper tail).

     \f[ D(x) = \int_{x}^{+\infty} p(x) dx  \f]

     where p(x) is the Landau probability density function.
     It is implemented simply as 1. - #landau_cdf

   @param x The argument \f$x\f$
   @param xi The width parameter \f$\xi\f$
   @param x0 The location parameter \f$x_0\f$

    @ingroup ProbFunc

  */
   inline double landau_cdf_c(double x, double xi = 1, double x0 = 0) {
      return 1. - landau_cdf(x,xi,x0);
   }

   /**

   Complement of the cumulative distribution function of the lognormal distribution
   (upper tail).

   \f[ D(x) = \int_{x}^{+\infty} {1 \over x' \sqrt{2 \pi s^2} } e^{-(\ln{x'} - m)^2/2 s^2} dx' \f]

   For detailed description see
   <A HREF="http://mathworld.wolfram.com/LogNormalDistribution.html">
   Mathworld</A>.

   @ingroup ProbFunc

   */

   double lognormal_cdf_c(double x, double m, double s, double x0 = 0);




   /**

   Cumulative distribution function of the lognormal distribution
   (lower tail).

   \f[ D(x) = \int_{-\infty}^{x} {1 \over x' \sqrt{2 \pi s^2} } e^{-(\ln{x'} - m)^2/2 s^2} dx' \f]

   For detailed description see
   <A HREF="http://mathworld.wolfram.com/LogNormalDistribution.html">
   Mathworld</A>.

   @ingroup ProbFunc

   */

   double lognormal_cdf(double x, double m, double s, double x0 = 0);




   /**

   Complement of the cumulative distribution function of the normal (Gaussian)
   distribution (upper tail).

   \f[ D(x) = \int_{x}^{+\infty} {1 \over \sqrt{2 \pi \sigma^2}} e^{-x'^2 / 2\sigma^2} dx' \f]

   For detailed description see
   <A HREF="http://mathworld.wolfram.com/NormalDistribution.html">
   Mathworld</A>.

   @ingroup ProbFunc

   */

   double normal_cdf_c(double x, double sigma = 1, double x0 = 0);
   /// Alternative name for same function
   inline double gaussian_cdf_c(double x, double sigma = 1, double x0 = 0) {
      return normal_cdf_c(x,sigma,x0);
   }



   /**

   Cumulative distribution function of the normal (Gaussian)
   distribution (lower tail).

   \f[ D(x) = \int_{-\infty}^{x} {1 \over \sqrt{2 \pi \sigma^2}} e^{-x'^2 / 2\sigma^2} dx' \f]

   For detailed description see
   <A HREF="http://mathworld.wolfram.com/NormalDistribution.html">
   Mathworld</A>.
   @ingroup ProbFunc

   */

   double normal_cdf(double x, double sigma = 1, double x0 = 0);
   /// Alternative name for same function
   inline double gaussian_cdf(double x, double sigma = 1, double x0 = 0) {
      return normal_cdf(x,sigma,x0);
   }



   /**

   Complement of the cumulative distribution function of Student's
   t-distribution (upper tail).

   \f[ D_{r}(x) = \int_{x}^{+\infty} \frac{\Gamma(\frac{r+1}{2})}{\sqrt{r \pi}\Gamma(\frac{r}{2})} \left( 1+\frac{x'^2}{r}\right)^{-(r+1)/2} dx' \f]

   For detailed description see
   <A HREF="http://mathworld.wolfram.com/Studentst-Distribution.html">
   Mathworld</A>. It is implemented using the incomplete beta function, ROOT::Math::inc_beta,
   from <A HREF="http://www.netlib.org/cephes">Cephes</A>

   @ingroup ProbFunc

   */

   double tdistribution_cdf_c(double x, double r, double x0 = 0);




   /**

   Cumulative distribution function of Student's
   t-distribution (lower tail).

   \f[ D_{r}(x) = \int_{-\infty}^{x} \frac{\Gamma(\frac{r+1}{2})}{\sqrt{r \pi}\Gamma(\frac{r}{2})} \left( 1+\frac{x'^2}{r}\right)^{-(r+1)/2} dx' \f]

   For detailed description see
   <A HREF="http://mathworld.wolfram.com/Studentst-Distribution.html">
   Mathworld</A>. It is implemented using the incomplete beta function, ROOT::Math::inc_beta,
   from <A HREF="http://www.netlib.org/cephes">Cephes</A>

   @ingroup ProbFunc

   */

   double tdistribution_cdf(double x, double r, double x0 = 0);


   /**

   Complement of the cumulative distribution function of the uniform (flat)
   distribution (upper tail).

   \f[ D(x) = \int_{x}^{+\infty} {1 \over (b-a)} dx' \f]

   For detailed description see
   <A HREF="http://mathworld.wolfram.com/UniformDistribution.html">
   Mathworld</A>.

   @ingroup ProbFunc

   */

   double uniform_cdf_c(double x, double a, double b, double x0 = 0);




   /**

   Cumulative distribution function of the uniform (flat)
   distribution (lower tail).

   \f[ D(x) = \int_{-\infty}^{x} {1 \over (b-a)} dx' \f]

   For detailed description see
   <A HREF="http://mathworld.wolfram.com/UniformDistribution.html">
   Mathworld</A>.

   @ingroup ProbFunc

   */

   double uniform_cdf(double x, double a, double b, double x0 = 0);




   /**

   Complement of the cumulative distribution function of the Poisson distribution.
   Upper tail of the integral of the #poisson_pdf

   @ingroup ProbFunc

   */

   double poisson_cdf_c(unsigned int n, double mu);

   /**

   Cumulative distribution function of the Poisson distribution
   Lower tail of the integral of the #poisson_pdf

   @ingroup ProbFunc

   */

   double poisson_cdf(unsigned int n, double mu);

   /**

   Complement of the cumulative distribution function of the Binomial distribution.
   Upper tail of the integral of the #binomial_pdf

   @ingroup ProbFunc

   */

   double binomial_cdf_c(unsigned int k, double p, unsigned int n);

   /**

   Cumulative distribution function of the Binomial distribution
   Lower tail of the integral of the #binomial_pdf

   @ingroup ProbFunc

   */

   double binomial_cdf(unsigned int k, double p, unsigned int n);


   /**

   Complement of the cumulative distribution function of the Negative Binomial distribution.
   Upper tail of the integral of the #negative_binomial_pdf

   @ingroup ProbFunc

   */

   double negative_binomial_cdf_c(unsigned int k, double p, double n);

   /**

   Cumulative distribution function of the Negative Binomial distribution
   Lower tail of the integral of the #negative_binomial_pdf

   @ingroup ProbFunc

   */

   double negative_binomial_cdf(unsigned int k, double p, double n);



#ifdef HAVE_OLD_STAT_FUNC

   /** @name Backward compatible MathCore CDF functions */


   inline double breitwigner_prob(double x, double gamma, double x0 = 0) {
      return  breitwigner_cdf_c(x,gamma,x0);
   }
   inline double breitwigner_quant(double x, double gamma, double x0 = 0) {
      return  breitwigner_cdf(x,gamma,x0);
   }

   inline double cauchy_prob(double x, double b, double x0 = 0) {
      return cauchy_cdf_c(x,b,x0);
   }
   inline double cauchy_quant(double x, double b, double x0 = 0) {
      return cauchy_cdf  (x,b,x0);
   }
   inline double chisquared_prob(double x, double r, double x0 = 0) {
      return chisquared_cdf_c(x, r, x0);
   }
   inline double chisquared_quant(double x, double r, double x0 = 0) {
      return chisquared_cdf  (x, r, x0);
   }
   inline double exponential_prob(double x, double lambda, double x0 = 0) {
      return exponential_cdf_c(x, lambda, x0 );
   }
   inline double exponential_quant(double x, double lambda, double x0 = 0) {
      return exponential_cdf  (x, lambda, x0 );
   }

   inline double gaussian_prob(double x, double sigma, double x0 = 0) {
      return  gaussian_cdf_c( x, sigma, x0 );
   }
   inline double gaussian_quant(double x, double sigma, double x0 = 0) {
      return  gaussian_cdf  ( x, sigma, x0 );
   }

   inline double lognormal_prob(double x, double m, double s, double x0 = 0) {
      return lognormal_cdf_c( x, m, s, x0 );
   }
   inline double lognormal_quant(double x, double m, double s, double x0 = 0) {
      return lognormal_cdf  ( x, m, s, x0 );
   }

   inline double normal_prob(double x, double sigma, double x0 = 0) {
      return  normal_cdf_c( x, sigma, x0 );
   }
   inline double normal_quant(double x, double sigma, double x0 = 0) {
      return  normal_cdf  ( x, sigma, x0 );
   }

   inline double uniform_prob(double x, double a, double b, double x0 = 0) {
      return uniform_cdf_c( x, a, b, x0 );
   }
   inline double uniform_quant(double x, double a, double b, double x0 = 0) {
      return uniform_cdf  ( x, a, b, x0 );
   }
   inline double fdistribution_prob(double x, double n, double m, double x0 = 0) {
      return fdistribution_cdf_c  (x, n, m, x0);
   }
   inline double fdistribution_quant(double x, double n, double m, double x0 = 0) {
      return fdistribution_cdf    (x, n, m, x0);
   }

   inline double gamma_prob(double x, double alpha, double theta, double x0 = 0) {
      return gamma_cdf_c (x, alpha, theta, x0);
   }
   inline double gamma_quant(double x, double alpha, double theta, double x0 = 0) {
      return gamma_cdf   (x, alpha, theta, x0);
   }

   inline double tdistribution_prob(double x, double r, double x0 = 0) {
      return tdistribution_cdf_c  (x, r, x0);
   }

   inline double tdistribution_quant(double x, double r, double x0 = 0) {
      return tdistribution_cdf    (x, r, x0);
   }

#endif

   /** @defgroup TruncFunc Statistical functions from truncated distributions

   @ingroup StatFunc

   Statistical functions for the truncated distributions. Examples of such functions are the
   first or the second momentum of the truncated distribution.
   In the case of the Landau, first and second momentum functions are provided for the Landau
   distribution truncated only on the right side.
   These functions are defined in the header file <em>Math/ProbFunc.h</em> or in the global one
   including all statistical functions <em>Math/StatFunc.h</em>

   */

   /**

   First moment (mean) of the truncated Landau distribution.
   \f[ \frac{1}{D (x)} \int_{-\infty}^{x} t\, p(t) d t \f]
   where \f$p(x)\f$ is the Landau distribution
   and \f$D(x)\f$ its cumulative distribution function.

   For detailed description see
   K.S. K&ouml;lbig and B. Schorr, A program package for the Landau distribution,
   <A HREF="http://dx.doi.org/10.1016/0010-4655(84)90085-7">Computer Phys. Comm. 31 (1984) 97-111</A>
   <A HREF="http://dx.doi.org/10.1016/j.cpc.2008.03.002">[Erratum-ibid. 178 (2008) 972]</A>.
   The same algorithms as in
   <A HREF="https://cern-tex.web.cern.ch/cern-tex/shortwrupsdir/g110/top.html">
   CERNLIB</A> (XM1LAN)  is used

   @param x The argument \f$x\f$
   @param xi The width parameter \f$\xi\f$
   @param x0 The location parameter \f$x_0\f$

   @ingroup TruncFunc

   */

   double landau_xm1(double x, double xi = 1, double x0 = 0);



   /**

   Second moment of the truncated Landau distribution.
   \f[ \frac{1}{D (x)} \int_{-\infty}^{x} t^2\, p(t) d t \f]
   where \f$p(x)\f$ is the Landau distribution
   and \f$D(x)\f$ its cumulative distribution function.

   For detailed description see
   K.S. K&ouml;lbig and B. Schorr, A program package for the Landau distribution,
   <A HREF="http://dx.doi.org/10.1016/0010-4655(84)90085-7">Computer Phys. Comm. 31 (1984) 97-111</A>
   <A HREF="http://dx.doi.org/10.1016/j.cpc.2008.03.002">[Erratum-ibid. 178 (2008) 972]</A>.
   The same algorithms as in
   <A HREF="https://cern-tex.web.cern.ch/cern-tex/shortwrupsdir/g110/top.html">
   CERNLIB</A> (XM1LAN)  is used

   @param x The argument \f$x\f$
   @param xi The width parameter \f$\xi\f$
   @param x0 The location parameter \f$x_0\f$

   @ingroup TruncFunc

   */

   double landau_xm2(double x, double xi = 1, double x0 = 0);



} // namespace Math
} // namespace ROOT


#endif // ROOT_Math_ProbFuncMathCore
