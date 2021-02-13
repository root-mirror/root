// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Minuit2_MnParabolaPoint
#define ROOT_Minuit2_MnParabolaPoint

namespace ROOT {

namespace Minuit2 {

/**

A point of a parabola.

<p>

????!!!! in reality it is just a general point in two dimensional space,
there is nothing that would indicate, that it belongs to a parabola.
This class defines simpy an (x,y) pair!!!!

@author Fred James and Matthias Winkler; comments added by Andras Zsenei
and Lorenzo Moneta

@ingroup Minuit

\todo Should it be called MnParabolaPoint or just Point?

 */

class MnParabolaPoint {

public:
   /**

   Initializes the point with its coordinates.

   @param x the x (first) coordinate of the point.
   @param y the y (second) coordinate of the point.

   */

   MnParabolaPoint(double x, double y) : fX(x), fY(y) {}

   ~MnParabolaPoint() {}

   /**

   Accessor to the x (first) coordinate.

   @return the x (first) coordinate of the point.

   */

   double X() const { return fX; }

   /**

   Accessor to the y (second) coordinate.

   @return the y (second) coordinate of the point.

   */

   double Y() const { return fY; }

private:
   double fX;
   double fY;
};

} // namespace Minuit2

} // namespace ROOT

#endif // ROOT_Minuit2_MnParabolaPoint
