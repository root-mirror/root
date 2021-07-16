// @(#)root/mathmore:$Id$
// Authors: L. Moneta, A. Zsenei   08/2005

// Authors: Andras Zsenei & Lorenzo Moneta   06/2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 , LCG ROOT MathLib Team                         *
 *                                                                    *
 *                                                                    *
 **********************************************************************/

#include <cmath>

#ifndef PI
#define PI       3.14159265358979323846264338328      /* pi */
#endif


#include "gsl/gsl_sf_bessel.h"
#include "gsl/gsl_sf_legendre.h"
#include "gsl/gsl_sf_lambert.h"
#include "gsl/gsl_sf_laguerre.h"
#include "gsl/gsl_sf_hyperg.h"
#include "gsl/gsl_sf_ellint.h"
//#include "gsl/gsl_sf_gamma.h"
#include "gsl/gsl_sf_expint.h"
#include "gsl/gsl_sf_zeta.h"
#include "gsl/gsl_sf_airy.h"
#include "gsl/gsl_sf_coupling.h"


namespace ROOT {
namespace Math {




// [5.2.1.1] associated Laguerre polynomials
// (26.x.12)

double assoc_laguerre(unsigned n, double m, double x) {

   return gsl_sf_laguerre_n(n, m, x);

}



// [5.2.1.2] associated Legendre functions
// (26.x.8)

double assoc_legendre(unsigned l, unsigned m, double x) {
   // add  (-1)^m
   return (m%2 == 0) ? gsl_sf_legendre_Plm(l, m, x) : -gsl_sf_legendre_Plm(l, m, x);

}

// Shortcut for RooFit to call the gsl legendre functions without the branches in the above implementation.
namespace internal{
  double legendre(unsigned l, unsigned m, double x) {
    return gsl_sf_legendre_Plm(l, m, x);
  }
}


// [5.2.1.4] (complete) elliptic integral of the first kind
// (26.x.15.2)

double comp_ellint_1(double k) {

   return gsl_sf_ellint_Kcomp(k, GSL_PREC_DOUBLE);

}



// [5.2.1.5] (complete) elliptic integral of the second kind
// (26.x.16.2)

double comp_ellint_2(double k) {

   return gsl_sf_ellint_Ecomp(k, GSL_PREC_DOUBLE);

}



// [5.2.1.6] (complete) elliptic integral of the third kind
// (26.x.17.2)
/**
Complete elliptic integral of the third kind.

There are two different definitions used for the elliptic
integral of the third kind:

\f[
P(\phi,k,n) = \int_0^\phi \frac{dt}{(1 + n \sin^2{t})\sqrt{1 - k^2 \sin^2{t}}}
\f]

and

\f[
P(\phi,k,n) = \int_0^\phi \frac{dt}{(1 - n \sin^2{t})\sqrt{1 - k^2 \sin^2{t}}}
\f]

the former is adopted by

- GSL
     http://www.gnu.org/software/gsl/manual/gsl-ref_7.html#SEC95

- Planetmath
     http://planetmath.org/encyclopedia/EllipticIntegralsAndJacobiEllipticFunctions.html

- CERNLIB
     https://cern-tex.web.cern.ch/cern-tex/shortwrupsdir/c346/top.html

   while the latter is used by

- Abramowitz and Stegun

- Mathematica
     http://mathworld.wolfram.com/EllipticIntegraloftheThirdKind.html

- C++ standard
     http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2004/n1687.pdf

   in order to be C++ compliant, we decided to use the latter, hence the
   change of the sign in the function call to GSL.

   */

double comp_ellint_3(double n, double k) {

   return gsl_sf_ellint_P(PI/2.0, k, -n, GSL_PREC_DOUBLE);

}



// [5.2.1.7] confluent hypergeometric functions
// (26.x.14)

double conf_hyperg(double a, double b, double z) {

   return gsl_sf_hyperg_1F1(a, b, z);

}

//  confluent hypergeometric functions of second type

double conf_hypergU(double a, double b, double z) {

   return gsl_sf_hyperg_U(a, b, z);

}



// [5.2.1.8] regular modified cylindrical Bessel functions
// (26.x.4.1)

double cyl_bessel_i(double nu, double x) {

   return gsl_sf_bessel_Inu(nu, x);

}



// [5.2.1.9] cylindrical Bessel functions (of the first kind)
// (26.x.2)

double cyl_bessel_j(double nu, double x) {

   return gsl_sf_bessel_Jnu(nu, x);

}



// [5.2.1.10] irregular modified cylindrical Bessel functions
// (26.x.4.2)

double cyl_bessel_k(double nu, double x) {

   return gsl_sf_bessel_Knu(nu, x);

}



// [5.2.1.11] cylindrical Neumann functions;
// cylindrical Bessel functions (of the second kind)
// (26.x.3)

double cyl_neumann(double nu, double x) {

   return gsl_sf_bessel_Ynu(nu, x);

}



// [5.2.1.12] (incomplete) elliptic integral of the first kind
// phi in radians
// (26.x.15.1)

double ellint_1(double k, double phi) {

   return gsl_sf_ellint_F(phi, k, GSL_PREC_DOUBLE);

}



// [5.2.1.13] (incomplete) elliptic integral of the second kind
// phi in radians
// (26.x.16.1)

double ellint_2(double k, double phi) {

   return gsl_sf_ellint_E(phi, k, GSL_PREC_DOUBLE);

}



// [5.2.1.14] (incomplete) elliptic integral of the third kind
// phi in radians
// (26.x.17.1)
/**

Incomplete elliptic integral of the third kind.

There are two different definitions used for the elliptic
integral of the third kind:

\f[
P(\phi,k,n) = \int_0^\phi \frac{dt}{(1 + n \sin^2{t})\sqrt{1 - k^2 \sin^2{t}}}
\f]

and

\f[
P(\phi,k,n) = \int_0^\phi \frac{dt}{(1 - n \sin^2{t})\sqrt{1 - k^2 \sin^2{t}}}
\f]

the former is adopted by

- GSL
     http://www.gnu.org/software/gsl/manual/gsl-ref_7.html#SEC95

- Planetmath
     http://planetmath.org/encyclopedia/EllipticIntegralsAndJacobiEllipticFunctions.html

- CERNLIB
     https://cern-tex.web.cern.ch/cern-tex/shortwrupsdir/c346/top.html

   while the latter is used by

- Abramowitz and Stegun

- Mathematica
     http://mathworld.wolfram.com/EllipticIntegraloftheThirdKind.html

- C++ standard
     http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2004/n1687.pdf

   in order to be C++ compliant, we decided to use the latter, hence the
   change of the sign in the function call to GSL.

   */

double ellint_3(double n, double k, double phi) {

   return gsl_sf_ellint_P(phi, k, -n, GSL_PREC_DOUBLE);

}



// [5.2.1.15] exponential integral
// (26.x.20)

double expint(double x) {

   return gsl_sf_expint_Ei(x);

}


// Generalization of expint(x)
//
double expint_n(int n, double x) {

   return gsl_sf_expint_En(n, x);

}



// [5.2.1.16] Hermite polynomials
// (26.x.10)

//double hermite(unsigned n, double x) {
//}




// [5.2.1.17] hypergeometric functions
// (26.x.13)

double hyperg(double a, double b, double c, double x) {

   return gsl_sf_hyperg_2F1(a, b, c, x);

}



// [5.2.1.18] Laguerre polynomials
// (26.x.11)

double laguerre(unsigned n, double x) {
   return gsl_sf_laguerre_n(n, 0, x);
}


// Lambert W function on branch 0

double lambert_W0(double x) {
   return gsl_sf_lambert_W0(x);
}


// Lambert W function on branch -1

double lambert_Wm1(double x) {
   return gsl_sf_lambert_Wm1(x);
}


// [5.2.1.19] Legendre polynomials
// (26.x.7)

double legendre(unsigned l, double x) {

   return gsl_sf_legendre_Pl(l, x);

}



// [5.2.1.20] Riemann zeta function
// (26.x.22)

double riemann_zeta(double x) {

   return gsl_sf_zeta(x);

}



// [5.2.1.21] spherical Bessel functions of the first kind
// (26.x.5)

double sph_bessel(unsigned n, double x) {

   return gsl_sf_bessel_jl(n, x);

}



// [5.2.1.22] spherical associated Legendre functions
// (26.x.9) ??????????

double sph_legendre(unsigned l, unsigned m, double theta) {

   return gsl_sf_legendre_sphPlm(l, m, std::cos(theta));

}




// [5.2.1.23] spherical Neumann functions
// (26.x.6)

double sph_neumann(unsigned n, double x) {

   return gsl_sf_bessel_yl(n, x);

}

// Airy function Ai

double airy_Ai(double x) {

   return gsl_sf_airy_Ai(x, GSL_PREC_DOUBLE);

}

// Airy function Bi

double airy_Bi(double x) {

   return gsl_sf_airy_Bi(x, GSL_PREC_DOUBLE);

}

// Derivative of the Airy function Ai

double airy_Ai_deriv(double x) {

   return gsl_sf_airy_Ai_deriv(x, GSL_PREC_DOUBLE);

}

// Derivative of the Airy function Bi

double airy_Bi_deriv(double x) {

   return gsl_sf_airy_Bi_deriv(x, GSL_PREC_DOUBLE);

}

// s-th zero of the Airy function Ai

double airy_zero_Ai(unsigned int s) {

   return gsl_sf_airy_zero_Ai(s);

}

// s-th zero of the Airy function Bi

double airy_zero_Bi(unsigned int s) {

   return gsl_sf_airy_zero_Bi(s);

}

// s-th zero of the derivative of the Airy function Ai

double airy_zero_Ai_deriv(unsigned int s) {

   return gsl_sf_airy_zero_Ai_deriv(s);

}

// s-th zero of the derivative of the Airy function Bi

double airy_zero_Bi_deriv(unsigned int s) {

   return gsl_sf_airy_zero_Bi_deriv(s);

}

// wigner coefficient

double wigner_3j(int ja, int jb, int jc, int ma, int mb, int mc) {
   return gsl_sf_coupling_3j(ja,jb,jc,ma,mb,mc);
}

double wigner_6j(int ja, int jb, int jc, int jd, int je, int jf) {
   return gsl_sf_coupling_6j(ja,jb,jc,jd,je,jf);
}

double wigner_9j(int ja, int jb, int jc, int jd, int je, int jf, int jg, int jh, int ji) {
   return gsl_sf_coupling_9j(ja,jb,jc,jd,je,jf,jg,jh,ji);
}

} // namespace Math
} // namespace ROOT
