
/* MathCore/tests/test_SpecFunc.cpp
 *
 * Copyright (C) 2004 Andras Zsenei
 *
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

 */


/**

Test file for the special functions implemented in MathMore. For
the moment nothing exceptional. Evaluates the functions and checks
if the value is right against values copied from the GSL tests.

*/





#include <iostream>
#include <iomanip>
#include <string>
#include <limits>

#include "Math/SpecFuncMathMore.h"
#ifndef NO_MATHCORE
#include "Math/SpecFuncMathCore.h"
#endif

#ifdef CHECK_WITH_GSL
#include <gsl/gsl_sf.h>
#endif

#ifndef PI
#define PI       3.14159265358979323846264338328      /* pi */
#endif


int compare( const std::string & name, double v1, double v2, double scale = 2.0) {
  //  ntest = ntest + 1;

   std::cout << std::setw(50) << std::left << name << ":\t";

  // numerical double limit for epsilon
   double eps = scale* std::numeric_limits<double>::epsilon();
   int iret = 0;
   double delta = v2 - v1;
   double d = 0;
   if (delta < 0 ) delta = - delta;
   if (v1 == 0 || v2 == 0) {
      if  (delta > eps ) {
         iret = 1;
      }
   }
   // skip case v1 or v2 is infinity
   else {
      d = v1;

      if ( v1 < 0) d = -d;
      // add also case when delta is small by default
      if ( delta/d  > eps && delta > eps )
         iret =  1;
   }

   if (iret == 0)
      std::cout <<" OK" << std::endl;
   else {
      std::cout <<" FAILED " << std::endl;
      int pr = std::cout.precision (18);
      std::cout << "\nDiscrepancy in " << name << "() :\n  " << v1 << " != " << v2 << " discr = " << int(delta/d/eps)
                << "   (Allowed discrepancy is " << eps  << ")\n\n";
      std::cout.precision (pr);
      //nfail = nfail + 1;
   }
   return iret;
}



// void showDiff(std::string name, double calculatedValue, double expectedValue, double scale = 1.0) {

//    compare(calculatedValue, expectedValue, name, scale)

//   std::cout << name << calculatedValue << " expected value: " << expectedValue;
//   int prec = std::cout.precision();
//   std::cout.precision(5);
//   std::cout << " diff: " << (calculatedValue-expectedValue) << " reldiff: " <<
//     (calculatedValue-expectedValue)/expectedValue << std::endl;
//   std::cout.precision(prec);

// }



int testSpecFunc() {

   using namespace ROOT::Math;

   int iret = 0;

   std::cout.precision(20);

#ifndef NO_MATHCORE

   // explicit put namespace to be sure to use right ones

   iret |= compare("tgamma(9.0) ", ROOT::Math::tgamma(9.0), 40320.0, 4);

   iret |= compare("lgamma(0.1) ", ROOT::Math::lgamma(0.1),  2.252712651734205, 4);

   iret |= compare("inc_gamma(1,0.001) ", ROOT::Math::inc_gamma(1.0,0.001), 0.0009995001666250083319, 1);

   // increase tolerance when using Cephes (test values are correctly checked with Mathematica
   // GSL was more precise in this case
   // Adapt also to 32 bits architectures
#if defined(R__B64)
      const int inc_gamma_scale = 100;
#else
      const int inc_gamma_scale = 200;
#endif

   iret |= compare("inc_gamma(100,99) ", ROOT::Math::inc_gamma(100.,99.), 0.4733043303994607, inc_gamma_scale);

   iret |= compare("inc_gamma_c(100,99) ", ROOT::Math::inc_gamma_c(100.,99.), 0.5266956696005394, inc_gamma_scale);

   // need to increase here by a further factor of 5 for Windows
   iret |= compare("inc_gamma_c(1000,1000.1) ", ROOT::Math::inc_gamma_c(1000.,1000.1),  0.4945333598559338247, 5000);

   iret |= compare("erf(0.5) ", ROOT::Math::erf(0.5), 0.5204998778130465377);

   iret |= compare("erfc(-1.0) ", ROOT::Math::erfc(-1.0), 1.8427007929497148693);

   iret |= compare("beta(1.0, 5.0) ", ROOT::Math::beta(1.0, 5.0), 0.2);

   iret |= compare("inc_beta(1,1,1) ", ROOT::Math::inc_beta(1.0, 1.0, 1.0), 1.0);

   iret |= compare("inc_beta(0.5,0.1,1.0) ", ROOT::Math::inc_beta(0.5, 0.1, 1.0), 0.9330329915368074160 );


#endif

   iret |= compare("assoc_laguerre(4,  2, 0.5) ", assoc_laguerre(4, 2, 0.5), 6.752604166666666667,8);

   iret |= compare("assoc_legendre(10, 1, -0.5) ", assoc_legendre(10, 1, -0.5), 2.0066877394361256516);

   iret |= compare("comp_ellint_1(0.50) ", comp_ellint_1(0.50), 1.6857503548125960429);

   iret |= compare("comp_ellint_2(0.50) ", comp_ellint_2(0.50), 1.4674622093394271555);

   iret |= compare("comp_ellint_3(0.5, 0.5) ", comp_ellint_3(0.5, 0.5), 2.41367150420119, 16);

   iret |= compare("conf_hyperg(1, 1.5, 1) ", conf_hyperg(1, 1.5, 1), 2.0300784692787049755);

   iret |= compare("cyl_bessel_i(1.0, 1.0) ", cyl_bessel_i(1.0, 1.0), 0.5651591039924850272);

   iret |= compare("cyl_bessel_j(0.75, 1.0) ", cyl_bessel_j(0.75, 1.0), 0.5586524932048917478, 16);

   iret |= compare("cyl_bessel_k(1.0, 1.0) ", cyl_bessel_k(1.0, 1.0), 0.6019072301972345747);

   iret |= compare("cyl_neumann(0.75, 1.0) ", cyl_neumann(0.75, 1.0), -0.6218694174429746383 );

   iret |= compare("ellint_1(0.50, PI/3.0) ", ellint_1(0.50, PI/3.0), 1.0895506700518854093);

   iret |= compare("ellint_2(0.50, PI/3.0) ", ellint_2(0.50, PI/3.0), 1.0075555551444720293);

   iret |= compare("ellint_3(-0.50, 0.5, PI/3.0) ", ellint_3(-0.50, 0.5, PI/3.0), 0.9570574331323584890);

   iret |= compare("expint(1.0) ", expint(1.0), 1.8951178163559367555);

   iret |= compare("expint_n(3, 0.4) ", expint_n(3, 0.4), 0.2572864233199447237);

   // std::cout << "Hermite polynomials: to do!" << std::endl;

   iret |= compare("hyperg(8, -8, 1, 0.5) ", hyperg(8, -8, 1, 0.5), 0.13671875);

   iret |= compare("laguerre(4, 1.) ", laguerre(4, 1.), -0.6250, 4); // need to find more precise value

   iret |= compare("lambert_W0(-0.1) ", lambert_W0(-0.1), -0.1118325591589629648);

   iret |= compare("lambert_Wm1(-0.1) ", lambert_Wm1(-0.1), -3.5771520639572972184);

   iret |= compare("legendre(10, -0.5) ", legendre(10, -0.5), -0.1882286071777345);

   iret |= compare("riemann_zeta(-0.5) ", riemann_zeta(-0.5), -0.207886224977354566017307, 16);

   iret |= compare("sph_bessel(1, 10.0) ", sph_bessel(1, 10.0), 0.07846694179875154709000);

   iret |= compare("sph_legendre(3, 1, PI/2.) ", sph_legendre(3, 1, PI/2.), 0.323180184114150653007);

   iret |= compare("sph_neumann(0, 1.0) ", sph_neumann(0, 1.0), -0.54030230586813972);

   iret |= compare("airy_Ai(-0.5) ", airy_Ai(-0.5), 0.475728091610539583);           // wolfram alpha:  0.47572809161053958880

   iret |= compare("airy_Bi(0.5) ", airy_Bi(0.5), 0.854277043103155553);             // wolfram alpha:  0.85427704310315549330

   iret |= compare("airy_Ai_deriv(-2) ", airy_Ai_deriv(-2), 0.618259020741691145);   // wolfram alpha:  0.61825902074169104141

   iret |= compare("airy_Bi_deriv(-3) ", airy_Bi_deriv(-3), -0.675611222685258639);  // wolfram alpha: -0.67561122268525853767

   iret |= compare("airy_zero_Ai(2) ", airy_zero_Ai(2), -4.08794944413097028, 8);    // mathworld: -4.08795

   iret |= compare("airy_zero_Bi(2) ", airy_zero_Bi(2), -3.27109330283635291, 8);    // mathworld: -3.27109

   iret |= compare("airy_zero_Ai_deriv(2) ", airy_zero_Ai_deriv(2), -3.24819758217983656, 8);

   iret |= compare("airy_zero_Bi_deriv(2) ", airy_zero_Bi_deriv(2), -4.07315508907182799, 8);

   if (iret != 0) {
      std::cout << "\n\nError:  Special Functions Test FAILED !!!!!" << std::endl;
   }
   return iret;

}

void getGSLErrors() {

#ifdef CHECK_WITH_GSL
   gsl_sf_result r;
   int iret;

   iret = gsl_sf_ellint_P_e(PI/2.0, 0.5, -0.5, GSL_PREC_DOUBLE, &r);
   std::cout << "comp_ellint_3(0.50, 0.5) : " << r.val << " err:  " << r.err << "  iret:  " << iret << std::endl;

   iret = gsl_sf_ellint_P_e(PI/3.0, 0.5, 0.5, GSL_PREC_DOUBLE, &r);
   std::cout << "ellint_3(0.50, 0.5, PI/3.0) : " << r.val << " err:  " << r.err << "  iret:  " << iret << std::endl;

   iret = gsl_sf_zeta_e(-0.5, &r);
   std::cout << "riemann_zeta(-0.5) : " << r.val << " err:  " << r.err << "  iret:  " << iret << std::endl;
#endif


}


int main() {

   int iret = 0;
   iret |=  testSpecFunc();

   getGSLErrors();
   return iret;
}


