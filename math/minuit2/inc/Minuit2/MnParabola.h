// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Minuit2_MnParabola
#define ROOT_Minuit2_MnParabola

#include <cmath>

namespace ROOT {

namespace Minuit2 {

/**

This class defines a parabola of the form a*x*x + b*x + c

@author Fred James and Matthias Winkler; comments added by Andras Zsenei
and Lorenzo Moneta

@ingroup Minuit

 */

class MnParabola {

public:
   /**

   Constructor that initializes the parabola with its three parameters.

   @param a the coefficient of the quadratic term
   @param b the coefficient of the linear term
   @param c the constant

   */

   MnParabola(double a, double b, double c) : fA(a), fB(b), fC(c) {}

   ~MnParabola() {}

   /**

   Evaluates the parabola a the point x.

   @param x the coordinate where the parabola needs to be evaluated.

   @return the y coordinate of the parabola corresponding to x.

   */

   double Y(double x) const { return (fA * x * x + fB * x + fC); }

   /**

   Calculates the bigger of the two x values corresponding to the
   given y Value.

   <p>

   ???????!!!!!!!!! And when there is none?? it looks like it will
   crash?? what is sqrt (-1.0) ?

   @param y the y Value for which the x Value is to be calculated.

   @return the bigger one of the two corresponding values.

   */

   // ok, at first glance it does not look like the formula for the quadratic
   // equation, but it is!  ;-)
   double X_pos(double y) const { return (std::sqrt(y / fA + Min() * Min() - fC / fA) + Min()); }
   // maybe it is worth to check the performance improvement with the below formula??
   //   double X_pos(double y) const {return (std::sqrt(y/fA + fB*fB/(4.*fA*fA) - fC/fA)  - fB/(2.*fA));}

   /**

   Calculates the smaller of the two x values corresponding to the
   given y Value.

   <p>

   ???????!!!!!!!!! And when there is none?? it looks like it will
   crash?? what is sqrt (-1.0) ?

   @param y the y Value for which the x Value is to be calculated.

   @return the smaller one of the two corresponding values.

   */

   double X_neg(double y) const { return (-std::sqrt(y / fA + Min() * Min() - fC / fA) + Min()); }

   /**

   Calculates the x coordinate of the Minimum of the parabola.

   @return x coordinate of the Minimum.

   */

   double Min() const { return -fB / (2. * fA); }

   /**

   Calculates the y coordinate of the Minimum of the parabola.

   @return y coordinate of the Minimum.

   */

   double YMin() const { return (-fB * fB / (4. * fA) + fC); }

   /**

   Accessor to the coefficient of the quadratic term.

   @return the coefficient of the quadratic term.

    */

   double A() const { return fA; }

   /**

   Accessor to the coefficient of the linear term.

   @return the coefficient of the linear term.

   */

   double B() const { return fB; }

   /**

   Accessor to the coefficient of the constant term.

   @return the coefficient of the constant term.

   */

   double C() const { return fC; }

private:
   double fA;
   double fB;
   double fC;
};

} // namespace Minuit2

} // namespace ROOT

#endif // ROOT_Minuit2_MnParabola
