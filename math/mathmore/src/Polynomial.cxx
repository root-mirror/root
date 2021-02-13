// @(#)root/mathmore:$Id$
// Authors: L. Moneta, A. Zsenei   08/2005

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

// Implementation file for class Polynomial
//
// Created by: Lorenzo Moneta  at Wed Nov 10 17:46:19 2004
//
// Last update: Wed Nov 10 17:46:19 2004
//

#include "Math/Polynomial.h"


#include "gsl/gsl_math.h"
#include "gsl/gsl_errno.h"
#include "gsl/gsl_poly.h"
#include "complex_quartic.h"

#define DEBUG
#ifdef DEBUG
#include <iostream>
#endif

namespace ROOT {
namespace Math {


Polynomial::Polynomial(unsigned int n) :
  // number of par is order + 1
  ParFunc( n+1 ),
  fOrder(n),
  fDerived_params(std::vector<double>(n) )
{
}

  //order 1
Polynomial::Polynomial(double a, double b) :
  ParFunc( 2 ),
  fOrder(1),
  fDerived_params(std::vector<double>(1) )
{
  fParams[0] = b;
  fParams[1] = a;
}

// order 2
Polynomial::Polynomial(double a, double b, double c) :
  ParFunc( 3 ),
  fOrder(2),
  fDerived_params(std::vector<double>(2) )
{
  fParams[0] = c;
  fParams[1] = b;
  fParams[2] = a;
}

// order 3 (cubic)
Polynomial::Polynomial(double a, double b, double c, double d) :
  // number of par is order + 1
  ParFunc( 4 ),
  fOrder(3),
  fDerived_params(std::vector<double>(3) )
{
  fParams[0] = d;
  fParams[1] = c;
  fParams[2] = b;
  fParams[3] = a;
}

// order 3 (quartic)
Polynomial::Polynomial(double a, double b, double c, double d, double e) :
  // number of par is order + 1
  ParFunc( 5 ),
  fOrder(4),
  fDerived_params(std::vector<double>(4) )
{
  fParams[0] = e;
  fParams[1] = d;
  fParams[2] = c;
  fParams[3] = b;
  fParams[4] = a;
}



// Polynomial::Polynomial(const Polynomial &)
// {
// }

// Polynomial & Polynomial::operator = (const Polynomial &rhs)
// {
//    if (this == &rhs) return *this;  // time saving self-test

//    return *this;
// }


double  Polynomial::DoEvalPar (double x, const double * p) const {

    return gsl_poly_eval( p, fOrder + 1, x);

}



double  Polynomial::DoDerivative(double x) const{

   for ( unsigned int i = 0; i < fOrder; ++i )
    fDerived_params[i] =  (i + 1) * Parameters()[i+1];

   return gsl_poly_eval( &fDerived_params.front(), fOrder, x);

}

double Polynomial::DoParameterDerivative (double x, const double *, unsigned int ipar) const {

      return gsl_pow_int(x, ipar);
}



IGenFunction * Polynomial::Clone() const {
    Polynomial * f =  new Polynomial(fOrder);
    f->fDerived_params = fDerived_params;
    f->SetParameters( Parameters() );
    return f;
}


const std::vector< std::complex <double> > &  Polynomial::FindRoots(){


    // check if order is correct
    unsigned int n = fOrder;
    while ( Parameters()[n] == 0 ) {
      n--;
    }

    fRoots.clear();
    fRoots.reserve(n);


    if (n == 0) {
      return fRoots;
    }
    else if (n == 1 ) {
      if ( Parameters()[1] == 0) return fRoots;
      double r = - Parameters()[0]/ Parameters()[1];
      fRoots.push_back( std::complex<double> ( r, 0.0) );
    }
    // quadratic equations
    else if (n == 2 ) {
      gsl_complex z1, z2;
      int nr = gsl_poly_complex_solve_quadratic(Parameters()[2], Parameters()[1], Parameters()[0], &z1, &z2);
      if ( nr != 2) {
#ifdef DEBUG
         std::cout << "Polynomial quadratic ::-  FAILED to find roots" << std::endl;
#endif
         return fRoots;
      }
      fRoots.push_back(std::complex<double>(z1.dat[0],z1.dat[1]) );
      fRoots.push_back(std::complex<double>(z2.dat[0],z2.dat[1]) );
    }
    // cubic equations
    else if (n == 3 ) {
      gsl_complex  z1, z2, z3;
      // renormmalize params in this case
      double w = Parameters()[3];
      double a = Parameters()[2]/w;
      double b = Parameters()[1]/w;
      double c = Parameters()[0]/w;
      int nr = gsl_poly_complex_solve_cubic(a, b, c, &z1, &z2, &z3);
      if (nr != 3) {
#ifdef DEBUG
         std::cout << "Polynomial  cubic::-  FAILED to find roots" << std::endl;
#endif
         return fRoots;

      }
      fRoots.push_back(std::complex<double> (z1.dat[0],z1.dat[1]) );
      fRoots.push_back(std::complex<double> (z2.dat[0],z2.dat[1]) );
      fRoots.push_back(std::complex<double> (z3.dat[0],z3.dat[1]) );
    }
    // quartic equations
    else if (n == 4 ) {
      // quartic eq.
      gsl_complex  z1, z2, z3, z4;
      // renormalize params in this case
      double w = Parameters()[4];
      double a = Parameters()[3]/w;
      double b = Parameters()[2]/w;
      double c = Parameters()[1]/w;
      double d = Parameters()[0]/w;
      int nr = gsl_poly_complex_solve_quartic(a, b, c, d, &z1, &z2, &z3, & z4);
      if (nr != 4) {
#ifdef DEBUG
         std::cout << "Polynomial quartic ::-  FAILED to find roots" << std::endl;
#endif
         return fRoots;
      }
      fRoots.push_back(std::complex<double> (z1.dat[0],z1.dat[1]) );
      fRoots.push_back(std::complex<double> (z2.dat[0],z2.dat[1]) );
      fRoots.push_back(std::complex<double> (z3.dat[0],z3.dat[1]) );
      fRoots.push_back(std::complex<double> (z4.dat[0],z4.dat[1]) );
    }
    else {
    // for higher order polynomial use numerical fRoots
      FindNumRoots();
    }

    return fRoots;

  }


std::vector< double >  Polynomial::FindRealRoots(){
  FindRoots();
  std::vector<double> roots;
  roots.reserve(fOrder);
  for (unsigned int i = 0; i < fOrder; ++i) {
    if (fRoots[i].imag() == 0)
      roots.push_back( fRoots[i].real() );
  }
  return roots;
}
const std::vector< std::complex <double> > &  Polynomial::FindNumRoots(){


    // check if order is correct
    unsigned int n = fOrder;
    while ( Parameters()[n] == 0 ) {
      n--;
    }
    fRoots.clear();
    fRoots.reserve(n);


    if (n == 0) {
      return fRoots;
    }

    gsl_poly_complex_workspace * w = gsl_poly_complex_workspace_alloc( n + 1);
    std::vector<double> z(2*n);
    int status = gsl_poly_complex_solve ( Parameters(), n+1, w, &z.front() );
    gsl_poly_complex_workspace_free(w);
    if (status != GSL_SUCCESS) return fRoots;
    for (unsigned int i = 0; i < n; ++i)
      fRoots.push_back(std::complex<double> (z[2*i],z[2*i+1] ) );

    return fRoots;
}


} // namespace Math
} // namespace ROOT
