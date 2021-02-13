// @(#)root/mathcore:$Id$
// Authors: David Gonzalez Maline    01/2008

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2006 , LCG ROOT MathLib Team                         *
 *                                                                    *
 *                                                                    *
 **********************************************************************/

#include "Math/RichardsonDerivator.h"
#include "Math/IFunctionfwd.h"
#include <cmath>
#include <limits>
#include <algorithm>

#include "Math/Error.h"


namespace ROOT {
namespace Math {

RichardsonDerivator::RichardsonDerivator(double h) :
   fFunctionCopied(false),
   fStepSize(h),
   fLastError(0),
   fFunction(0)
{
   // Default Constructor.
}
RichardsonDerivator::RichardsonDerivator(const ROOT::Math::IGenFunction & f, double h, bool copyFunc) :
   fFunctionCopied(copyFunc),
   fStepSize(h),
   fLastError(0),
   fFunction(0)
{
   // Constructor from a function and step size
   if (copyFunc) fFunction = f.Clone();
   else fFunction = &f;
}


RichardsonDerivator::~RichardsonDerivator()
{
   // destructor
   if ( fFunction != 0 && fFunctionCopied )
      delete fFunction;
}

RichardsonDerivator::RichardsonDerivator(const RichardsonDerivator & rhs)
{
    // copy constructor
    // copy constructor (deep copy or not depending on fFunctionCopied)
   fStepSize = rhs.fStepSize;
   fLastError = rhs.fLastError;
   fFunctionCopied = rhs.fFunctionCopied;
   SetFunction(*rhs.fFunction);
 }

RichardsonDerivator &  RichardsonDerivator::operator= ( const RichardsonDerivator & rhs)
{
   // Assignment operator
   if (&rhs == this) return *this;
   fFunctionCopied = rhs.fFunctionCopied;
   fStepSize = rhs.fStepSize;
   fLastError = rhs.fLastError;
   SetFunction(*rhs.fFunction);
   return *this;
}

void RichardsonDerivator::SetFunction(const ROOT::Math::IGenFunction & f)
{
   // set function
   if (fFunctionCopied) {
      if (fFunction) delete fFunction;
      fFunction = f.Clone();
   }
   else fFunction = &f;
}

double RichardsonDerivator::Derivative1 (const IGenFunction & function, double x, double h)
{
   const double keps = std::numeric_limits<double>::epsilon();


   double xx;
   xx = x+h;     double f1 = (function)(xx);

   xx = x-h;     double f2 = (function)(xx);

   xx = x+h/2;   double g1 = (function)(xx);

   xx = x-h/2;   double g2 = (function)(xx);

   //compute the central differences
   double h2    = 1/(2.*h);
   double d0    = f1 - f2;
   double d2    = g1 - g2;
   double deriv = h2*(8*d2 - d0)/3.;
   // // compute the error ( to be improved ) this is just a simple truncation error
   // fLastError   = kC1*h2*0.5*(f1+f2);  //compute the error

   // compute the error ( from GSL deriv implementation)

   double e0 = (std::abs( f1) + std::abs(f2)) * keps;
   double e2 = 2* (std::abs( g1) + std::abs(g2)) * keps + e0;
   double delta = std::max( std::abs( h2*d0), std::abs( deriv) ) * std::abs( x)/h * keps;

   // estimate the truncation error from d2-d0 which is O(h^2)
   double err_trunc = std::abs( deriv - h2*d0 );
   // rounding error due to cancellation
   double err_round = std::abs( e2/h) + delta;

   fLastError = err_trunc + err_round;
   return deriv;
}

double RichardsonDerivator::DerivativeForward (const IGenFunction & function, double x, double h)
{
   const double keps = std::numeric_limits<double>::epsilon();


   double xx;
   xx = x+h/4.0;         double f1 = (function)(xx);
   xx = x+h/2.0;         double f2 = (function)(xx);

   xx = x+(3.0/4.0)*h;   double f3 = (function)(xx);
   xx = x+h;             double f4 = (function)(xx);

   //compute the forward differences

   double r2 = 2.0*(f4 - f2);
   double r4 = (22.0 / 3.0) * (f4 - f3) - (62.0 / 3.0) * (f3 - f2) +
    (52.0 / 3.0) * (f2 - f1);

   // Estimate the rounding error for r4

   double e4 = 2 * 20.67 * (fabs (f4) + fabs (f3) + fabs (f2) + fabs (f1)) * keps;

   // The next term is due to finite precision in x+h = O (eps * x)

   double dy = std::max (fabs (r2 / h), fabs (r4 / h)) * fabs (x / h) * keps;

   // The truncation error in the r4 approximation itself is O(h^3).
   //  However, for safety, we estimate the error from r4-r2, which is
   //  O(h).  By scaling h we will minimise this estimated error, not
   //  the actual truncation error in r4.

   double result = r4 / h;
   double abserr_trunc = fabs ((r4 - r2) / h); // Estimated truncation error O(h)
   double abserr_round = fabs (e4 / h) + dy;

   fLastError = abserr_trunc + abserr_round;

   return result;
}


   double RichardsonDerivator::Derivative2 (const IGenFunction & function, double x, double h)
{
   const double kC1 = 4*std::numeric_limits<double>::epsilon();

   double xx;
   xx = x+h;     double f1 = (function)(xx);
   xx = x;       double f2 = (function)(xx);
   xx = x-h;     double f3 = (function)(xx);

   xx = x+h/2;   double g1 = (function)(xx);
   xx = x-h/2;   double g3 = (function)(xx);

   //compute the central differences
   double hh    = 1/(h*h);
   double d0    = f3 - 2*f2 + f1;
   double d2    = 4*g3 - 8*f2 +4*g1;
   fLastError   = kC1*hh*f2;  //compute the error
   double deriv = hh*(4*d2 - d0)/3.;
   return deriv;
}

double RichardsonDerivator::Derivative3 (const IGenFunction & function, double x, double h)
{
   const double kC1 = 4*std::numeric_limits<double>::epsilon();

   double xx;
   xx = x+2*h;   double f1 = (function)(xx);
   xx = x+h;     double f2 = (function)(xx);
   xx = x-h;     double f3 = (function)(xx);
   xx = x-2*h;   double f4 = (function)(xx);
   xx = x;       double fx = (function)(xx);
   xx = x+h/2;   double g2 = (function)(xx);
   xx = x-h/2;   double g3 = (function)(xx);

   //compute the central differences
   double hhh  = 1/(h*h*h);
   double d0   = 0.5*f1 - f2 +f3 - 0.5*f4;
   double d2   = 4*f2 - 8*g2 +8*g3 - 4*f3;
   fLastError    = kC1*hhh*fx;   //compute the error
   double deriv = hhh*(4*d2 - d0)/3.;
   return deriv;
}

} // end namespace Math

} // end namespace ROOT
