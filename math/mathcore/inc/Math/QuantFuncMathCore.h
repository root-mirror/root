// @(#)root/mathcore:$Id$
// Authors: L. Moneta, A. Zsenei   08/2005


// Authors: Andras Zsenei & Lorenzo Moneta   08/2005


/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 , LCG ROOT MathLib Team                         *
 *                                                                    *
 *                                                                    *
 **********************************************************************/


#ifndef ROOT_Math_QuantFuncMathCore
#define ROOT_Math_QuantFuncMathCore


namespace ROOT {
namespace Math {



  /** @defgroup QuantFunc Quantile Functions
   *  @ingroup StatFunc
   *
   *  Inverse functions of the cumulative distribution functions
   *  and the inverse of the complement of the cumulative distribution functions
   *  for various distributions.
   *  The functions with the extension <em>_quantile</em> calculate the
   *  inverse of the <em>_cdf</em> function, the
   *  lower tail integral of the probability density function
   *  \f$D^{-1}(z)\f$ where
   *
   *  \f[ D(x) = \int_{-\infty}^{x} p(x') dx' \f]
   *
   *  while those with the <em>_quantile_c</em> extension calculate the
   *  inverse of the <em>_cdf_c</em> functions, the upper tail integral of the probability
   *  density function \f$D^{-1}(z) \f$ where
   *
   *  \f[ D(x) = \int_{x}^{+\infty} p(x') dx' \f]
   *
   *  These functions are defined in the header file <em>Math/ProbFunc.h</em> or in the global one
   *  including all statistical functions <em>Math/DistFunc.h</em>
   *
   *
   * <strong>NOTE:</strong> In the old releases (< 5.14) the <em>_quantile</em> functions were called
   * <em>_quant_inv</em> and the <em>_quantile_c</em> functions were called
   * <em>_prob_inv</em>.
   * These names are currently kept for backward compatibility, but
   * their usage is deprecated.
   *
   */

   /** @name Quantile Functions from MathCore
   * The implementation is provided in MathCore and for the majority of the function comes from
   * <A HREF="http://www.netlib.org/cephes">Cephes</A>.

   */

  //@{



  /**

     Inverse (\f$D^{-1}(z)\f$) of the cumulative distribution
     function of the upper tail of the beta distribution
     (#beta_cdf_c).
     It is implemented using the function incbi from <A HREF="http://www.netlib.org/cephes">Cephes</A>.


     @ingroup QuantFunc

  */
   double beta_quantile(double x, double a, double b);

   /**

      Inverse (\f$D^{-1}(z)\f$) of the cumulative distribution
      function of the lower tail of the beta distribution
      (#beta_cdf).
      It is implemented using
      the function incbi from <A HREF="http://www.netlib.org/cephes">Cephes</A>.

      @ingroup QuantFunc

   */
   double beta_quantile_c(double x, double a, double b);



   /**

      Inverse (\f$D^{-1}(z)\f$) of the cumulative distribution
      function of the upper tail of the Cauchy distribution (#cauchy_cdf_c)
      which is also called Lorentzian distribution. For
      detailed description see
      <A HREF="http://mathworld.wolfram.com/CauchyDistribution.html">
      Mathworld</A>.

      @ingroup QuantFunc

   */

   double cauchy_quantile_c(double z, double b);




   /**

      Inverse (\f$D^{-1}(z)\f$) of the cumulative distribution
      function of the lower tail of the Cauchy distribution (#cauchy_cdf)
      which is also called Breit-Wigner or Lorentzian distribution. For
      detailed description see
      <A HREF="http://mathworld.wolfram.com/CauchyDistribution.html">
      Mathworld</A>. The implementation used is that of
      <A HREF="http://www.gnu.org/software/gsl/manual/gsl-ref_19.html#SEC294">GSL</A>.

      @ingroup QuantFunc

   */

   double cauchy_quantile(double z, double b);




   /**

      Inverse (\f$D^{-1}(z)\f$) of the cumulative distribution
      function of the upper tail of the Breit-Wigner distribution (#breitwigner_cdf_c)
      which is similar to the Cauchy distribution. For
      detailed description see
      <A HREF="http://mathworld.wolfram.com/CauchyDistribution.html">
      Mathworld</A>. It is evaluated using the same implementation of
      #cauchy_quantile_c.

      @ingroup QuantFunc

   */

   inline double breitwigner_quantile_c(double z, double gamma) {
      return cauchy_quantile_c(z, gamma/2.0);
   }




   /**

      Inverse (\f$D^{-1}(z)\f$) of the cumulative distribution
      function of the lower tail of the Breit_Wigner distribution (#breitwigner_cdf)
      which is similar to the Cauchy distribution. For
      detailed description see
      <A HREF="http://mathworld.wolfram.com/CauchyDistribution.html">
      Mathworld</A>. It is evaluated using the same implementation of
      #cauchy_quantile.


      @ingroup QuantFunc

   */

   inline double breitwigner_quantile(double z, double gamma) {
      return cauchy_quantile(z, gamma/2.0);
   }






   /**

      Inverse (\f$D^{-1}(z)\f$) of the cumulative distribution
      function of the upper tail of the \f$\chi^2\f$ distribution
      with \f$r\f$ degrees of freedom (#chisquared_cdf_c). For detailed description see
      <A HREF="http://mathworld.wolfram.com/Chi-SquaredDistribution.html">
      Mathworld</A>. It is implemented using the inverse of the incomplete complement gamma function, using
      the function igami from <A HREF="http://www.netlib.org/cephes">Cephes</A>.

      @ingroup QuantFunc

   */

   double chisquared_quantile_c(double z, double r);



   /**

      Inverse (\f$D^{-1}(z)\f$) of the cumulative distribution
      function of the lower tail of the \f$\chi^2\f$ distribution
      with \f$r\f$ degrees of freedom (#chisquared_cdf). For detailed description see
      <A HREF="http://mathworld.wolfram.com/Chi-SquaredDistribution.html">
      Mathworld</A>.
      It is implemented using  chisquared_quantile_c, therefore is not very precise for small z.
      It is recommended to use the MathMore function (ROOT::MathMore::chisquared_quantile )implemented using GSL

      @ingroup QuantFunc

   */

   double chisquared_quantile(double z, double r);



   /**

      Inverse (\f$D^{-1}(z)\f$) of the cumulative distribution
      function of the upper tail of the exponential distribution
      (#exponential_cdf_c). For detailed description see
      <A HREF="http://mathworld.wolfram.com/ExponentialDistribution.html">
      Mathworld</A>.

      @ingroup QuantFunc

   */

   double exponential_quantile_c(double z, double lambda);




   /**

      Inverse (\f$D^{-1}(z)\f$) of the cumulative distribution
      function of the lower tail of the exponential distribution
      (#exponential_cdf). For detailed description see
      <A HREF="http://mathworld.wolfram.com/ExponentialDistribution.html">
      Mathworld</A>.

      @ingroup QuantFunc

   */

   double exponential_quantile(double z, double lambda);



   /**

      Inverse (\f$D^{-1}(z)\f$) of the cumulative distribution
      function of the lower tail of the f distribution
      (#fdistribution_cdf). For detailed description see
      <A HREF="http://mathworld.wolfram.com/F-Distribution.html">
      Mathworld</A>.
      It is implemented using the inverse of the incomplete beta function,
      function incbi from <A HREF="http://www.netlib.org/cephes">Cephes</A>.

      @ingroup QuantFunc

   */
   double fdistribution_quantile(double z, double n, double m);

   /**

      Inverse (\f$D^{-1}(z)\f$) of the cumulative distribution
      function of the upper tail of the f distribution
      (#fdistribution_cdf_c). For detailed description see
      <A HREF="http://mathworld.wolfram.com/F-Distribution.html">
      Mathworld</A>.
      It is implemented using the inverse of the incomplete beta function,
      function incbi from <A HREF="http://www.netlib.org/cephes">Cephes</A>.

      @ingroup QuantFunc
   */

   double fdistribution_quantile_c(double z, double n, double m);


   /**

      Inverse (\f$D^{-1}(z)\f$) of the cumulative distribution
      function of the upper tail of the gamma distribution
      (#gamma_cdf_c). For detailed description see
      <A HREF="http://mathworld.wolfram.com/GammaDistribution.html">
      Mathworld</A>. The implementation used is that of
      <A HREF="http://www.gnu.org/software/gsl/manual/gsl-ref_19.html#SEC300">GSL</A>.
      It is implemented using the function igami taken
      from <A HREF="http://www.netlib.org/cephes">Cephes</A>.

      @ingroup QuantFunc

   */

   double gamma_quantile_c(double z, double alpha, double theta);




   /**

      Inverse (\f$D^{-1}(z)\f$) of the cumulative distribution
      function of the lower tail of the gamma distribution
      (#gamma_cdf). For detailed description see
      <A HREF="http://mathworld.wolfram.com/GammaDistribution.html">
      Mathworld</A>.
      It is implemented using  chisquared_quantile_c, therefore is not very precise for small z.
      For this special cases it is recommended to use the MathMore function ROOT::MathMore::gamma_quantile
      implemented using GSL


      @ingroup QuantFunc

   */

   double gamma_quantile(double z, double alpha, double theta);



   /**

      Inverse (\f$D^{-1}(z)\f$) of the cumulative distribution
      function of the upper tail of the normal (Gaussian) distribution
      (#gaussian_cdf_c). For detailed description see
      <A HREF="http://mathworld.wolfram.com/NormalDistribution.html">
      Mathworld</A>. It can also be evaluated using #normal_quantile_c which will
      call the same implementation.

      @ingroup QuantFunc

   */

   double gaussian_quantile_c(double z, double sigma);




   /**

      Inverse (\f$D^{-1}(z)\f$) of the cumulative distribution
      function of the lower tail of the normal (Gaussian) distribution
      (#gaussian_cdf). For detailed description see
      <A HREF="http://mathworld.wolfram.com/NormalDistribution.html">
      Mathworld</A>. It can also be evaluated using #normal_quantile which will
      call the same implementation.
      It is implemented using the function  ROOT::Math::Cephes::ndtri taken from
      <A HREF="http://www.netlib.org/cephes">Cephes</A>.

      @ingroup QuantFunc

   */

   double gaussian_quantile(double z, double sigma);




   /**

      Inverse (\f$D^{-1}(z)\f$) of the cumulative distribution
      function of the upper tail of the lognormal distribution
      (#lognormal_cdf_c). For detailed description see
      <A HREF="http://mathworld.wolfram.com/LogNormalDistribution.html">
      Mathworld</A>. The implementation used is that of
      <A HREF="http://www.gnu.org/software/gsl/manual/gsl-ref_19.html#SEC302">GSL</A>.

      @ingroup QuantFunc

   */

   double lognormal_quantile_c(double x, double m, double s);




   /**

      Inverse (\f$D^{-1}(z)\f$) of the cumulative distribution
      function of the lower tail of the lognormal distribution
      (#lognormal_cdf). For detailed description see
      <A HREF="http://mathworld.wolfram.com/LogNormalDistribution.html">
      Mathworld</A>. The implementation used is that of
      <A HREF="http://www.gnu.org/software/gsl/manual/gsl-ref_19.html#SEC302">GSL</A>.

      @ingroup QuantFunc

   */

   double lognormal_quantile(double x, double m, double s);




   /**

      Inverse (\f$D^{-1}(z)\f$) of the cumulative distribution
      function of the upper tail of the normal (Gaussian) distribution
      (#normal_cdf_c). For detailed description see
      <A HREF="http://mathworld.wolfram.com/NormalDistribution.html">
      Mathworld</A>. It can also be evaluated using #gaussian_quantile_c which will
      call the same implementation.
      It is implemented using the function  ROOT::Math::Cephes::ndtri taken from
      <A HREF="http://www.netlib.org/cephes">Cephes</A>.

      @ingroup QuantFunc

   */

   double normal_quantile_c(double z, double sigma);
   /// alternative name for same function
   inline double gaussian_quantile_c(double z, double sigma) {
      return normal_quantile_c(z,sigma);
   }




   /**

      Inverse (\f$D^{-1}(z)\f$) of the cumulative distribution
      function of the lower tail of the normal (Gaussian) distribution
      (#normal_cdf). For detailed description see
      <A HREF="http://mathworld.wolfram.com/NormalDistribution.html">
      Mathworld</A>. It can also be evaluated using #gaussian_quantile which will
      call the same implementation.
      It is implemented using the function  ROOT::Math::Cephes::ndtri taken from
      <A HREF="http://www.netlib.org/cephes">Cephes</A>.


      @ingroup QuantFunc

   */

   double normal_quantile(double z, double sigma);
   /// alternative name for same function
   inline double gaussian_quantile(double z, double sigma) {
      return normal_quantile(z,sigma);
   }



#ifdef LATER // t quantiles are still in MathMore

   /**

      Inverse (\f$D^{-1}(z)\f$) of the cumulative distribution
      function of the upper tail of Student's t-distribution
      (#tdistribution_cdf_c). For detailed description see
      <A HREF="http://mathworld.wolfram.com/Studentst-Distribution.html">
      Mathworld</A>. The implementation used is that of
      <A HREF="http://www.gnu.org/software/gsl/manual/gsl-ref_19.html#SEC305">GSL</A>.

      @ingroup QuantFunc

   */

   double tdistribution_quantile_c(double z, double r);




   /**

      Inverse (\f$D^{-1}(z)\f$) of the cumulative distribution
      function of the lower tail of Student's t-distribution
      (#tdistribution_cdf). For detailed description see
      <A HREF="http://mathworld.wolfram.com/Studentst-Distribution.html">
      Mathworld</A>. The implementation used is that of
      <A HREF="http://www.gnu.org/software/gsl/manual/gsl-ref_19.html#SEC305">GSL</A>.

      @ingroup QuantFunc

   */

   double tdistribution_quantile(double z, double r);

#endif


   /**

      Inverse (\f$D^{-1}(z)\f$) of the cumulative distribution
      function of the upper tail of the uniform (flat) distribution
      (#uniform_cdf_c). For detailed description see
      <A HREF="http://mathworld.wolfram.com/UniformDistribution.html">
      Mathworld</A>.

      @ingroup QuantFunc

   */

   double uniform_quantile_c(double z, double a, double b);




   /**

      Inverse (\f$D^{-1}(z)\f$) of the cumulative distribution
      function of the lower tail of the uniform (flat) distribution
      (#uniform_cdf). For detailed description see
      <A HREF="http://mathworld.wolfram.com/UniformDistribution.html">
      Mathworld</A>.

      @ingroup QuantFunc

   */

   double uniform_quantile(double z, double a, double b);




   /**

      Inverse (\f$D^{-1}(z)\f$) of the cumulative distribution
      function of the lower tail of the Landau distribution
      (#landau_cdf).

   For detailed description see
   K.S. K&ouml;lbig and B. Schorr, A program package for the Landau distribution,
   <A HREF="http://dx.doi.org/10.1016/0010-4655(84)90085-7">Computer Phys. Comm. 31 (1984) 97-111</A>
   <A HREF="http://dx.doi.org/10.1016/j.cpc.2008.03.002">[Erratum-ibid. 178 (2008) 972]</A>.
   The same algorithms as in
   <A HREF="https://cern-tex.web.cern.ch/cern-tex/shortwrupsdir/g110/top.html">
   CERNLIB</A> (RANLAN) is used.

   @param z The argument \f$z\f$
   @param xi The width parameter \f$\xi\f$

      @ingroup QuantFunc

   */

   double landau_quantile(double z, double xi = 1);


   /**

      Inverse (\f$D^{-1}(z)\f$) of the cumulative distribution
      function of the upper tail of the landau distribution
      (#landau_cdf_c).
      Implemented using #landau_quantile

   @param z The argument \f$z\f$
   @param xi The width parameter \f$\xi\f$

      @ingroup QuantFunc

   */

   double landau_quantile_c(double z, double xi = 1);


#ifdef HAVE_OLD_STAT_FUNC

   //@}
   /** @name Backward compatible functions */


   inline double breitwigner_prob_inv(double x, double gamma) {
      return  breitwigner_quantile_c(x,gamma);
   }
   inline double breitwigner_quant_inv(double x, double gamma) {
      return  breitwigner_quantile(x,gamma);
   }

   inline double cauchy_prob_inv(double x, double b) {
      return cauchy_quantile_c(x,b);
   }
   inline double cauchy_quant_inv(double x, double b) {
      return cauchy_quantile  (x,b);
   }

   inline double exponential_prob_inv(double x, double lambda) {
      return exponential_quantile_c(x, lambda );
   }
   inline double exponential_quant_inv(double x, double lambda) {
      return exponential_quantile  (x, lambda );
   }

   inline double gaussian_prob_inv(double x, double sigma) {
      return  gaussian_quantile_c( x, sigma );
   }
   inline double gaussian_quant_inv(double x, double sigma) {
      return  gaussian_quantile  ( x, sigma );
   }

   inline double lognormal_prob_inv(double x, double m, double s) {
      return lognormal_quantile_c( x, m, s );
   }
   inline double lognormal_quant_inv(double x, double m, double s) {
      return lognormal_quantile  ( x, m, s );
   }

   inline double normal_prob_inv(double x, double sigma) {
      return  normal_quantile_c( x, sigma );
   }
   inline double normal_quant_inv(double x, double sigma) {
      return  normal_quantile  ( x, sigma );
   }

   inline double uniform_prob_inv(double x, double a, double b) {
      return uniform_quantile_c( x, a, b );
   }
   inline double uniform_quant_inv(double x, double a, double b) {
      return uniform_quantile  ( x, a, b );
   }

   inline double chisquared_prob_inv(double x, double r) {
      return chisquared_quantile_c(x, r );
   }

   inline double gamma_prob_inv(double x, double alpha, double theta) {
      return gamma_quantile_c (x, alpha, theta );
   }


#endif


} // namespace Math
} // namespace ROOT



#endif // ROOT_Math_QuantFuncMathCore
