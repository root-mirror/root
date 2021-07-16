// @(#)root/mathmore:$Id$
// Authors: L. Moneta, A. Zsenei   08/2005

// Authors: Andras Zsenei & Lorenzo Moneta   06/2005

 /**********************************************************************
  *                                                                    *
  * Copyright (c) 2004 ROOT Foundation,  CERN/PH-SFT                   *
  *                                                                    *
  * This library is free software; you can redistribute it and/or      *
  * modify it under the terms of the GNU General Public License        *
  * as published by the Free Software Foundation; either version 2     *
  * of the License, or (at your option) any later version.             *
  *                                                                    *
  * This library is distributed in the hope that it will be useful,    *
  * but WITHOUT ANY WARRANTY; without even the implied warranty of     *
  * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU   *
  * General Public License for more details.                           *
  *                                                                    *
  * You should have received a copy of the GNU General Public License  *
  * along with this library (see file COPYING); if not, write          *
  * to the Free Software Foundation, Inc., 59 Temple Place, Suite      *
  * 330, Boston, MA 02111-1307 USA, or contact the author.             *
  *                                                                    *
  **********************************************************************/

/**

Special mathematical functions.
The naming and numbering of the functions is taken from
<A HREF="http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2004/n1687.pdf">
Matt Austern,
(Draft) Technical Report on Standard Library Extensions,
N1687=04-0127, September 10, 2004</A>

@author Created by Andras Zsenei on Mon Nov 8 2004

@defgroup SpecFunc Special functions

*/





#ifndef ROOT_Math_SpecFuncMathMore
#define ROOT_Math_SpecFuncMathMore




namespace ROOT {
namespace Math {

   /** @name Special Functions from MathMore */


  /**


  Computes the generalized Laguerre polynomials for
  \f$ n \geq 0, m > -1 \f$.
  They are defined in terms of the confluent hypergeometric function.
  For integer values of m they can be defined in terms of the Laguerre polynomials \f$L_n(x)\f$:

  \f[ L_{n}^{m}(x) = (-1)^{m} \frac{d^m}{dx^m} L_{n+m}(x) \f]


  For detailed description see
  <A HREF="http://mathworld.wolfram.com/LaguerrePolynomial.html">Mathworld</A>.
  The implementation used is that of
  <A HREF="http://www.gnu.org/software/gsl/manual/html_node/Laguerre-Functions.html">GSL</A>.

  This function is an extension of C++0x, also consistent in GSL,
  Abramowitz and Stegun 1972 and MatheMathica that uses non-integer values for m.
  C++0x calls for 'int m', more restrictive than necessary.
  The definition for was incorrect in 'n1687.pdf', but fixed in
  <A HREF="http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2005/n1836.pdf">n1836.pdf</A>,
  the most recent draft as of 2007-07-01


  @ingroup SpecFunc

  */
  // [5.2.1.1] associated Laguerre polynomials

  double assoc_laguerre(unsigned n, double m, double x);




  /**

  Computes the associated Legendre polynomials.

  \f[ P_{l}^{m}(x) = (1-x^2)^{m/2} \frac{d^m}{dx^m} P_{l}(x) \f]

  with \f$m \geq 0\f$, \f$ l \geq m \f$ and \f$ |x|<1 \f$.
  There are two sign conventions for associated Legendre polynomials.
  As is the case with the above formula, some authors (e.g., Arfken
  1985, pp. 668-669) omit the Condon-Shortley phase \f$(-1)^m\f$,
  while others include it (e.g., Abramowitz and Stegun 1972).
  One possible way to distinguish the two conventions is due to
  Abramowitz and Stegun (1972, p. 332), who use the notation

  \f[ P_{lm} (x) = (-1)^m P_{l}^{m} (x)\f]

  to distinguish the two. For detailed description see
  <A HREF="http://mathworld.wolfram.com/LegendrePolynomial.html">
  Mathworld</A>. The implementation used is that of
  <A HREF="http://www.gnu.org/software/gsl/manual/html_node/Associated-Legendre-Polynomials-and-Spherical-Harmonics.html">GSL</A>.

  The definition uses is the one of C++0x, \f$ P_{lm}\f$, while GSL and MatheMatica use the \f$P_{l}^{m}\f$ definition. Note that C++0x and GSL definitions agree instead for the normalized associated Legendre polynomial,
  sph_legendre(l,m,theta).

  @ingroup SpecFunc

  */
  // [5.2.1.2] associated Legendre functions

  double assoc_legendre(unsigned l, unsigned m, double x);

  // Shortcut for RooFit to call the gsl legendre functions without the branches in the above implementation.
  namespace internal{
    double legendre(unsigned l, unsigned m, double x);
  }


  /**

  Calculates the complete elliptic integral of the first kind.

  \f[ K(k) = F(k, \pi / 2) = \int_{0}^{\pi /2} \frac{d \theta}{\sqrt{1 - k^2 \sin^2{\theta}}} \f]

  with \f$0 \leq k^2 \leq 1\f$. For detailed description see
  <A HREF="http://mathworld.wolfram.com/CompleteEllipticIntegraloftheFirstKind.html">
  Mathworld</A>. The implementation used is that of
  <A HREF="http://www.gnu.org/software/gsl/manual/gsl-ref_7.html#SEC100">GSL</A>.

  @ingroup SpecFunc

  */
  // [5.2.1.4] (complete) elliptic integral of the first kind

  double comp_ellint_1(double k);




  /**

  Calculates the complete elliptic integral of the second kind.

  \f[ E(k) = E(k , \pi / 2) = \int_{0}^{\pi /2} \sqrt{1 - k^2 \sin^2{\theta}} d \theta  \f]

  with \f$0 \leq k^2 \leq 1\f$. For detailed description see
  <A HREF="http://mathworld.wolfram.com/CompleteEllipticIntegraloftheSecondKind.html">
  Mathworld</A>. The implementation used is that of
  <A HREF="http://www.gnu.org/software/gsl/manual/gsl-ref_7.html#SEC100">GSL</A>.

  @ingroup SpecFunc

  */
  // [5.2.1.5] (complete) elliptic integral of the second kind

  double comp_ellint_2(double k);




  /**

  Calculates the complete elliptic integral of the third kind.

  \f[ \Pi (n, k, \pi / 2) = \int_{0}^{\pi /2} \frac{d \theta}{(1 - n \sin^2{\theta})\sqrt{1 - k^2 \sin^2{\theta}}}  \f]

  with \f$0 \leq k^2 \leq 1\f$. There are two sign conventions
  for elliptic integrals of the third kind. Some authors (Abramowitz
  and Stegun,
  <A HREF="http://mathworld.wolfram.com/EllipticIntegraloftheThirdKind.html">
  Mathworld</A>,
  <A HREF="http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2004/n1687.pdf">
  C++ standard proposal</A>) use the above formula, while others
  (<A HREF="http://www.gnu.org/software/gsl/manual/gsl-ref_7.html#SEC95">
  GSL</A>, <A HREF="http://planetmath.org/encyclopedia/EllipticIntegralsAndJacobiEllipticFunctions.html">
  Planetmath</A> and
  <A HREF="https://cern-tex.web.cern.ch/cern-tex/shortwrupsdir/c346/top.html">
  CERNLIB</A>) use the + sign in front of n in the denominator. In
  order to be C++ compliant, the present library uses the former
  convention. The implementation used is that of
  <A HREF="http://www.gnu.org/software/gsl/manual/gsl-ref_7.html#SEC101">GSL</A>.

  @ingroup SpecFunc

  */
  // [5.2.1.6] (complete) elliptic integral of the third kind
  double comp_ellint_3(double n, double k);




  /**

  Calculates the confluent hypergeometric functions of the first kind.

  \f[ _{1}F_{1}(a;b;z) = \frac{\Gamma(b)}{\Gamma(a)} \sum_{n=0}^{\infty} \frac{\Gamma(a+n)}{\Gamma(b+n)} \frac{z^n}{n!}  \f]

  For detailed description see
  <A HREF="http://mathworld.wolfram.com/ConfluentHypergeometricFunctionoftheFirstKind.html">
  Mathworld</A>. The implementation used is that of
  <A HREF="http://www.gnu.org/software/gsl/manual/gsl-ref_7.html#SEC125">GSL</A>.

  @ingroup SpecFunc

  */
  // [5.2.1.7] confluent hypergeometric functions

  double conf_hyperg(double a, double b, double z);


  /**

  Calculates the confluent hypergeometric functions of the second kind, known also as Kummer function of the second kind,
  it is related to the confluent hypergeometric functions of the first kind.

  \f[ U(a,b,z)  = \frac{ \pi}{ \sin{\pi b } } \left[ \frac{ _{1}F_{1}(a,b,z) } {\Gamma(a-b+1) }
            - \frac{ z^{1-b} { _{1}F_{1}}(a-b+1,2-b,z)}{\Gamma(a)} \right]  \f]

  For detailed description see
  <A HREF="http://mathworld.wolfram.com/ConfluentHypergeometricFunctionoftheSecondKind.html">
  Mathworld</A>. The implementation used is that of
  <A HREF="http://www.gnu.org/software/gsl/manual/gsl-ref_7.html#SEC125">GSL</A>.
  This function is not part of the C++ standard proposal

  @ingroup SpecFunc

  */
  // confluent hypergeometric functions of second type

  double conf_hypergU(double a, double b, double z);



  /**

  Calculates the modified Bessel function of the first kind
  (also called regular modified (cylindrical) Bessel function).

  \f[ I_{\nu} (x) = i^{-\nu} J_{\nu}(ix) = \sum_{k=0}^{\infty} \frac{(\frac{1}{2}x)^{\nu + 2k}}{k! \Gamma(\nu + k + 1)} \f]

  for \f$x>0, \nu > 0\f$. For detailed description see
  <A HREF="http://mathworld.wolfram.com/ModifiedBesselFunctionoftheFirstKind.html">
  Mathworld</A>. The implementation used is that of
  <A HREF="http://www.gnu.org/software/gsl/manual/gsl-ref_7.html#SEC71">GSL</A>.

  @ingroup SpecFunc

  */
  // [5.2.1.8] regular modified cylindrical Bessel functions

  double cyl_bessel_i(double nu, double x);




  /**

  Calculates the (cylindrical) Bessel functions of the first kind (also
  called regular (cylindrical) Bessel functions).

  \f[ J_{\nu} (x) = \sum_{k=0}^{\infty} \frac{(-1)^k(\frac{1}{2}x)^{\nu + 2k}}{k! \Gamma(\nu + k + 1)} \f]

  For detailed description see
  <A HREF="http://mathworld.wolfram.com/BesselFunctionoftheFirstKind.html">
  Mathworld</A>. The implementation used is that of
  <A HREF="http://www.gnu.org/software/gsl/manual/gsl-ref_7.html#SEC69">GSL</A>.

  @ingroup SpecFunc

  */
  // [5.2.1.9] cylindrical Bessel functions (of the first kind)

  double cyl_bessel_j(double nu, double x);





  /**

  Calculates the modified Bessel functions of the second kind
  (also called irregular modified (cylindrical) Bessel functions).

  \f[ K_{\nu} (x) = \frac{\pi}{2} i^{\nu + 1} (J_{\nu} (ix) + iN(ix)) = \left\{ \begin{array}{cl} \frac{\pi}{2} \frac{I_{-\nu}(x) - I_{\nu}(x)}{\sin{\nu \pi}} & \mbox{for non-integral $\nu$} \\  \frac{\pi}{2} \lim{\mu \to \nu} \frac{I_{-\mu}(x) - I_{\mu}(x)}{\sin{\mu \pi}}
& \mbox{for integral $\nu$} \end{array}  \right.  \f]

  for \f$x>0, \nu > 0\f$. For detailed description see
  <A HREF="http://mathworld.wolfram.com/ModifiedBesselFunctionoftheSecondKind.html">
  Mathworld</A>. The implementation used is that of
  <A HREF="http://www.gnu.org/software/gsl/manual/gsl-ref_7.html#SEC72">GSL</A>.

  @ingroup SpecFunc

  */
  // [5.2.1.10] irregular modified cylindrical Bessel functions

  double cyl_bessel_k(double nu, double x);




  /**

  Calculates the (cylindrical) Bessel functions of the second kind
  (also called irregular (cylindrical) Bessel functions or
  (cylindrical) Neumann functions).

  \f[ N_{\nu} (x) = Y_{\nu} (x) = \left\{ \begin{array}{cl} \frac{J_{\nu} \cos{\nu \pi}-J_{-\nu}(x)}{\sin{\nu \pi}}  & \mbox{for non-integral $\nu$} \\ \lim{\mu \to \nu} \frac{J_{\mu} \cos{\mu \pi}-J_{-\mu}(x)}{\sin{\mu \pi}}  & \mbox{for integral $\nu$} \end{array} \right.  \f]

   For detailed description see
  <A HREF="http://mathworld.wolfram.com/BesselFunctionoftheSecondKind.html">
  Mathworld</A>. The implementation used is that of
  <A HREF="http://www.gnu.org/software/gsl/manual/gsl-ref_7.html#SEC70">GSL</A>.

  @ingroup SpecFunc

  */
  // [5.2.1.11] cylindrical Neumann functions;
  // cylindrical Bessel functions (of the second kind)

  double cyl_neumann(double nu, double x);




  /**

  Calculates the incomplete elliptic integral of the first kind.

  \f[ F(k, \phi) = \int_{0}^{\phi} \frac{d \theta}{\sqrt{1 - k^2 \sin^2{\theta}}} \f]

  with \f$0 \leq k^2 \leq 1\f$. For detailed description see
  <A HREF="http://mathworld.wolfram.com/EllipticIntegraloftheFirstKind.html">
  Mathworld</A>. The implementation used is that of
  <A HREF="http://www.gnu.org/software/gsl/manual/gsl-ref_7.html#SEC101">GSL</A>.

  @param k
  @param phi angle in radians

  @ingroup SpecFunc

  */
  // [5.2.1.12] (incomplete) elliptic integral of the first kind
  // phi in radians

  double ellint_1(double k, double phi);




  /**

  Calculates the complete elliptic integral of the second kind.

  \f[ E(k , \phi) = \int_{0}^{\phi} \sqrt{1 - k^2 \sin^2{\theta}} d \theta  \f]

  with \f$0 \leq k^2 \leq 1\f$. For detailed description see
  <A HREF="http://mathworld.wolfram.com/EllipticIntegraloftheSecondKind.html">
  Mathworld</A>. The implementation used is that of
  <A HREF="http://www.gnu.org/software/gsl/manual/gsl-ref_7.html#SEC101">GSL</A>.

  @param k
  @param phi angle in radians

  @ingroup SpecFunc

  */
  // [5.2.1.13] (incomplete) elliptic integral of the second kind
  // phi in radians

  double ellint_2(double k, double phi);




  /**

  Calculates the complete elliptic integral of the third kind.

  \f[ \Pi (n, k, \phi) = \int_{0}^{\phi} \frac{d \theta}{(1 - n \sin^2{\theta})\sqrt{1 - k^2 \sin^2{\theta}}}  \f]

  with \f$0 \leq k^2 \leq 1\f$. There are two sign conventions
  for elliptic integrals of the third kind. Some authors (Abramowitz
  and Stegun,
  <A HREF="http://mathworld.wolfram.com/EllipticIntegraloftheThirdKind.html">
  Mathworld</A>,
  <A HREF="http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2004/n1687.pdf">
  C++ standard proposal</A>) use the above formula, while others
  (<A HREF="http://www.gnu.org/software/gsl/manual/gsl-ref_7.html#SEC95">
  GSL</A>, <A HREF="http://planetmath.org/encyclopedia/EllipticIntegralsAndJacobiEllipticFunctions.html">
  Planetmath</A> and
  <A HREF="https://cern-tex.web.cern.ch/cern-tex/shortwrupsdir/c346/top.html">
  CERNLIB</A>) use the + sign in front of n in the denominator. In
  order to be C++ compliant, the present library uses the former
  convention. The implementation used is that of
  <A HREF="http://www.gnu.org/software/gsl/manual/gsl-ref_7.html#SEC101">GSL</A>.

  @param n
  @param k
  @param phi angle in radians

  @ingroup SpecFunc

  */
  // [5.2.1.14] (incomplete) elliptic integral of the third kind
  // phi in radians

  double ellint_3(double n, double k, double phi);




  /**

  Calculates the exponential integral.

  \f[ Ei(x) = - \int_{-x}^{\infty} \frac{e^{-t}}{t} dt \f]

  For detailed description see
  <A HREF="http://mathworld.wolfram.com/ExponentialIntegral.html">
  Mathworld</A>. The implementation used is that of
  <A HREF="http://www.gnu.org/software/gsl/manual/gsl-ref_7.html#SEC115">GSL</A>.

  @ingroup SpecFunc

  */
  // [5.2.1.15] exponential integral

  double expint(double x);
  double expint_n(int n, double x);



  // [5.2.1.16] Hermite polynomials

  //double hermite(unsigned n, double x);





  /**

  Calculates Gauss' hypergeometric function.

  \f[ _{2}F_{1}(a,b;c;x) = \frac{\Gamma(c)}{\Gamma(a) \Gamma(b)} \sum_{n=0}^{\infty} \frac{\Gamma(a+n)\Gamma(b+n)}{\Gamma(c+n)} \frac{x^n}{n!} \f]

  For detailed description see
  <A HREF="http://mathworld.wolfram.com/HypergeometricFunction.html">
  Mathworld</A>. The implementation used is that of
  <A HREF="http://www.gnu.org/software/gsl/manual/gsl-ref_7.html#SEC125">GSL</A>.

  @ingroup SpecFunc

  */
  // [5.2.1.17] hypergeometric functions

  double hyperg(double a, double b, double c, double x);



  /**

  Calculates the Laguerre polynomials

  \f[ P_{l}(x) = \frac{ e^x}{n!} \frac{d^n}{dx^n}  (x^n - e^{-x}) \f]

  for \f$x \geq 0 \f$ in the Rodrigues representation.
  They corresponds to the associated Laguerre polynomial of order m=0.
  See Abramowitz and Stegun, (22.5.16)
  For detailed description see
  <A HREF="http://mathworld.wolfram.com/LaguerrePolynomial.html">
  Mathworld</A>.
  The are implemented using the associated Laguerre polynomial of order m=0.

  @ingroup SpecFunc

  */
  // [5.2.1.18] Laguerre polynomials

  double laguerre(unsigned n, double x);


  /**

  Calculates the Lambert W function on branch 0

  The Lambert W functions are defined to be the solution of the equation

  \f[ W(x) \exp(W(x)) = x \f]

  For detailed description see
  <A HREF="https://mathworld.wolfram.com/LambertW-Function.html">Mathworld</A>
  or <A HREF="https://en.wikipedia.org/wiki/Lambert_W_function">Wikipedia</A>.

  This function implements the Lambert W function on branch 0, which is real
  valued and defined for \f$ x \geq -1/e \f$ with \f$ W_0(x) \geq -1 \f$.

  @ingroup SpecFunc

  */

  double lambert_W0(double x);


  /**

  Calculates the Lambert W function on branch -1

  The Lambert W functions are defined to be the solution of the equation

  \f[ W(x) \exp(W(x)) = x \f]

  For detailed description see
  <A HREF="https://mathworld.wolfram.com/LambertW-Function.html">Mathworld</A>
  or <A HREF="https://en.wikipedia.org/wiki/Lambert_W_function">Wikipedia</A>.

  This function implements the Lambert W function on branch -1, which is real
  valued and defined for \f$ -1/e \seq x < 0 \f$ with
  \f$ W_{-1}(x) \seq -1 \f$.

  @ingroup SpecFunc

  */

  double lambert_Wm1(double x);


  /**

  Calculates the Legendre polynomials.

  \f[ P_{l}(x) = \frac{1}{2^l l!} \frac{d^l}{dx^l}  (x^2 - 1)^l \f]

  for \f$l \geq 0, |x|\leq1\f$ in the Rodrigues representation.
  For detailed description see
  <A HREF="http://mathworld.wolfram.com/LegendrePolynomial.html">
  Mathworld</A>. The implementation used is that of
  <A HREF="http://www.gnu.org/software/gsl/manual/gsl-ref_7.html#SEC129">GSL</A>.

  @ingroup SpecFunc

  */
  // [5.2.1.19] Legendre polynomials

  double legendre(unsigned l, double x);




  /**

  Calculates the Riemann zeta function.

  \f[ \zeta (x) = \left\{ \begin{array}{cl} \sum_{k=1}^{\infty}k^{-x} & \mbox{for $x > 1$} \\ 2^x \pi^{x-1} \sin{(\frac{1}{2}\pi x)} \Gamma(1-x) \zeta (1-x) & \mbox{for $x < 1$} \end{array} \right. \f]

  For detailed description see
  <A HREF="http://mathworld.wolfram.com/RiemannZetaFunction.html">
  Mathworld</A>. The implementation used is that of
  <A HREF="http://www.gnu.org/software/gsl/manual/gsl-ref_7.html#SEC149">GSL</A>.

  CHECK WHETHER THE IMPLEMENTATION CALCULATES X<1

  @ingroup SpecFunc

  */
  // [5.2.1.20] Riemann zeta function

  double riemann_zeta(double x);


  /**

  Calculates the spherical Bessel functions of the first kind
  (also called regular spherical Bessel functions).

  \f[ j_{n}(x) = \sqrt{\frac{\pi}{2x}} J_{n+1/2}(x) \f]

  For detailed description see
  <A HREF="http://mathworld.wolfram.com/SphericalBesselFunctionoftheFirstKind.html">
  Mathworld</A>. The implementation used is that of
  <A HREF="http://www.gnu.org/software/gsl/manual/gsl-ref_7.html#SEC73">GSL</A>.

  @ingroup SpecFunc

  */
  // [5.2.1.21] spherical Bessel functions of the first kind

  double sph_bessel(unsigned n, double x);


  /**

  Computes the spherical (normalized)  associated Legendre polynomials,
  or spherical harmonic without azimuthal dependence (\f$e^(im\phi)\f$).

  \f[ Y_l^m(theta,0) = \sqrt{(2l+1)/(4\pi)} \sqrt{(l-m)!/(l+m)!} P_l^m(cos \theta) \f]

  for \f$m \geq 0, l \geq m\f$,
  where the Condon-Shortley phase  \f$(-1)^m\f$ is included in P_l^m(x)
  This function is consistent with both C++0x and GSL,
  even though there is a discrepancy in where to include the phase.
  There is no reference in Abramowitz and Stegun.


  @ingroup SpecFunc

  */

  // [5.2.1.22] spherical associated Legendre functions

  double sph_legendre(unsigned l, unsigned m, double theta);


  /**

  Calculates the spherical Bessel functions of the second kind
  (also called irregular spherical Bessel functions or
  spherical Neumann functions).

  \f[ n_n(x) = y_n(x) = \sqrt{\frac{\pi}{2x}} N_{n+1/2}(x) \f]

  For detailed description see
  <A HREF="http://mathworld.wolfram.com/SphericalBesselFunctionoftheSecondKind.html">
  Mathworld</A>. The implementation used is that of
  <A HREF="http://www.gnu.org/software/gsl/manual/gsl-ref_7.html#SEC74">GSL</A>.

  @ingroup SpecFunc

  */
  // [5.2.1.23] spherical Neumann functions

  double sph_neumann(unsigned n, double x);

  /**

  Calculates the Airy function Ai

  \f[ Ai(x) = \frac{1}{\pi} \int\limits_{0}^{\infty} \cos(xt + t^3/3) dt \f]

  For detailed description see
  <A HREF="http://mathworld.wolfram.com/AiryFunctions.html">
  Mathworld</A>
  and <A HREF="http://www.nrbook.com/abramowitz_and_stegun/page_446.htm">Abramowitz&Stegun, Sect. 10.4</A>.
  The implementation used is that of
  <A HREF="http://www.gnu.org/software/gsl/manual/html_node/Airy-Functions.html">GSL</A>.

  @ingroup SpecFunc

  */
  // Airy function Ai

  double airy_Ai(double x);

  /**

  Calculates the Airy function Bi

  \f[ Bi(x) = \frac{1}{\pi} \int\limits_{0}^{\infty} [\exp(xt-t^3/3) + \cos(xt + t^3/3)] dt \f]

  For detailed description see
  <A HREF="http://mathworld.wolfram.com/AiryFunctions.html">
  Mathworld</A>
  and <A HREF="http://www.nrbook.com/abramowitz_and_stegun/page_446.htm">Abramowitz&Stegun, Sect. 10.4</A>.
  The implementation used is that of
  <A HREF="http://www.gnu.org/software/gsl/manual/html_node/Airy-Functions.html">GSL</A>.

  @ingroup SpecFunc

  */
  // Airy function Bi

  double airy_Bi(double x);

  /**

  Calculates the derivative of the Airy function Ai

  \f[ Ai(x) = \frac{1}{\pi} \int\limits_{0}^{\infty} \cos(xt + t^3/3) dt \f]

  For detailed description see
  <A HREF="http://mathworld.wolfram.com/AiryFunctions.html">
  Mathworld</A>
  and <A HREF="http://www.nrbook.com/abramowitz_and_stegun/page_446.htm">Abramowitz&Stegun, Sect. 10.4</A>.
  The implementation used is that of
  <A HREF="http://www.gnu.org/software/gsl/manual/html_node/Derivatives-of-Airy-Functions.html">GSL</A>.

  @ingroup SpecFunc

  */
  // Derivative of the Airy function Ai

  double airy_Ai_deriv(double x);

  /**

  Calculates the derivative of the Airy function Bi

  \f[ Bi(x) = \frac{1}{\pi} \int\limits_{0}^{\infty} [\exp(xt-t^3/3) + \cos(xt + t^3/3)] dt \f]

  For detailed description see
  <A HREF="http://mathworld.wolfram.com/AiryFunctions.html">
  Mathworld</A>
  and <A HREF="http://www.nrbook.com/abramowitz_and_stegun/page_446.htm">Abramowitz&Stegun, Sect. 10.4</A>.
  The implementation used is that of
  <A HREF="http://www.gnu.org/software/gsl/manual/html_node/Derivatives-of-Airy-Functions.html">GSL</A>.

  @ingroup SpecFunc

  */
  // Derivative of the Airy function Bi

  double airy_Bi_deriv(double x);

  /**

  Calculates the zeroes of the Airy function Ai

  \f[ Ai(x) = \frac{1}{\pi} \int\limits_{0}^{\infty} \cos(xt + t^3/3) dt \f]

  For detailed description see
  <A HREF="http://mathworld.wolfram.com/AiryFunctionZeros.html">
  Mathworld</A>
  and <A HREF="http://www.nrbook.com/abramowitz_and_stegun/page_446.htm">Abramowitz&Stegun, Sect. 10.4</A>.
  The implementation used is that of
  <A HREF="http://www.gnu.org/software/gsl/manual/html_node/Zeros-of-Airy-Functions.html">GSL</A>.

  @ingroup SpecFunc

  */
  // s-th zero of the Airy function Ai

  double airy_zero_Ai(unsigned int s);

  /**

  Calculates the zeroes of the Airy function Bi

  \f[ Bi(x) = \frac{1}{\pi} \int\limits_{0}^{\infty} [\exp(xt-t^3/3) + \cos(xt + t^3/3)] dt \f]

  For detailed description see
  <A HREF="http://mathworld.wolfram.com/AiryFunctionZeros.html">
  Mathworld</A>
  and <A HREF="http://www.nrbook.com/abramowitz_and_stegun/page_446.htm">Abramowitz&Stegun, Sect. 10.4</A>.
  The implementation used is that of
  <A HREF="http://www.gnu.org/software/gsl/manual/html_node/Zeros-of-Airy-Functions.html">GSL</A>.

  @ingroup SpecFunc

  */
  // s-th zero of the Airy function Bi

  double airy_zero_Bi(unsigned int s);

  /**

  Calculates the zeroes of the derivative of the Airy function Ai

  \f[ Ai(x) = \frac{1}{\pi} \int\limits_{0}^{\infty} \cos(xt + t^3/3) dt \f]

  For detailed description see
  <A HREF="http://mathworld.wolfram.com/AiryFunctionZeros.html">
  Mathworld</A>
  and <A HREF="http://www.nrbook.com/abramowitz_and_stegun/page_446.htm">Abramowitz&Stegun, Sect. 10.4</A>.
  The implementation used is that of
  <A HREF="http://www.gnu.org/software/gsl/manual/html_node/Zeros-of-Derivatives-of-Airy-Functions.html">GSL</A>.

  @ingroup SpecFunc

  */
  // s-th zero of the derivative of the Airy function Ai

  double airy_zero_Ai_deriv(unsigned int s);

  /**

  Calculates the zeroes of the derivative of the Airy function Bi

  \f[ Bi(x) = \frac{1}{\pi} \int\limits_{0}^{\infty} [\exp(xt-t^3/3) + \cos(xt + t^3/3)] dt \f]

  For detailed description see
  <A HREF="http://mathworld.wolfram.com/AiryFunctionZeros.html">
  Mathworld</A>
  and <A HREF="http://www.nrbook.com/abramowitz_and_stegun/page_446.htm">Abramowitz&Stegun, Sect. 10.4</A>.
  The implementation used is that of
  <A HREF="http://www.gnu.org/software/gsl/manual/html_node/Zeros-of-Derivatives-of-Airy-Functions.html">GSL</A>.

  @ingroup SpecFunc

  */
  // s-th zero of the derivative of the Airy function Bi

  double airy_zero_Bi_deriv(unsigned int s);

  /**

   Calculates the Wigner 3j coupling coefficients

           (ja jb jc
           ma mb mc)

   where ja,ma,...etc are integers or half integers.
   The function takes as input arguments only integers which corresponds
   to half integer units, e.g two_ja = 2 * ja

   For detailed description see
   <A HREF="http://mathworld.wolfram.com/Wigner3j-Symbol.html.html">
   Mathworld</A>.
   The implementation used is that of
   <A HREF="http://www.gnu.org/software/gsl/manual/html_node/3_002dj-Symbols.html#g_t3_002dj-Symbols">GSL</A>.

   @ingroup SpecFunc

  */

  double wigner_3j(int two_ja, int two_jb, int two_jc, int two_ma, int two_mb, int two_mc);

  /**

  Calculates the Wigner 6j coupling coefficients

          (ja jb jc
           jd je jf)

   where ja,jb,...etc are integers or half integers.
   The function takes as input arguments only integers which corresponds
   to half integer units, e.g two_ja = 2 * ja

   For detailed description see
   <A HREF="http://mathworld.wolfram.com/Wigner6j-Symbol.html">
   Mathworld</A>.
   The implementation used is that of
   <A HREF="http://www.gnu.org/software/gsl/manual/html_node/6_002dj-Symbols.html#g_t6_002dj-Symbols">GSL</A>.

   @ingroup SpecFunc

  */

  double wigner_6j(int two_ja, int two_jb, int two_jc, int two_jd, int two_je, int two_jf);

  /**

  Calculates the Wigner 9j coupling coefficients

          (ja jb jc
           jd je jf
           jg jh ji)

   where ja,jb...etc are integers or half integers.
   The function takes as input arguments only integers which corresponds
   to half integer units, e.g two_ja = 2 * ja


   For detailed description see
   <A HREF="http://mathworld.wolfram.com/Wigner9j-Symbol.html">
   Mathworld</A>.
   The implementation used is that of
   <A HREF="http://www.gnu.org/software/gsl/manual/html_node/9_002dj-Symbols.html#g_t9_002dj-Symbols">GSL</A>.

   @ingroup SpecFunc

  */

   double wigner_9j(int two_ja, int two_jb, int two_jc, int two_jd, int two_je, int two_jf, int two_jg, int two_jh, int two_ji);



} // namespace Math
} // namespace ROOT


#endif //ROOT_Math_SpecFuncMathMore
