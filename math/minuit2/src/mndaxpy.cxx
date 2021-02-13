// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

/* daxpy.f -- translated by f2c (version 20010320).
   You must link the resulting object file with the libraries:
      -lf2c -lm   (in that order)
*/

namespace ROOT {

namespace Minuit2 {

int Mndaxpy(unsigned int n, double da, const double *dx, int incx, double *dy, int incy)
{
   /* System generated locals */
   int i__1;

   /* Local variables */
   int i__, m, ix, iy, mp1;

   /*     constant times a vector plus a vector. */
   /*     uses unrolled loops for increments equal to one. */
   /*     jack dongarra, linpack, 3/11/78. */
   /*     modified 12/3/93, array(1) declarations changed to array(*) */

   /* Parameter adjustments */
   --dy;
   --dx;

   /* Function Body */
   if (n <= 0) {
      return 0;
   }
   if (da == 0.) {
      return 0;
   }
   if (incx == 1 && incy == 1) {
      goto L20;
   }

   /*        code for unequal increments or equal increments */
   /*          not equal to 1 */

   ix = 1;
   iy = 1;
   if (incx < 0) {
      ix = (-static_cast<int>(n) + 1) * incx + 1;
   }
   if (incy < 0) {
      iy = (-static_cast<int>(n) + 1) * incy + 1;
   }
   i__1 = n;
   for (i__ = 1; i__ <= i__1; ++i__) {
      dy[iy] += da * dx[ix];
      ix += incx;
      iy += incy;
      /* L10: */
   }
   return 0;

   /*        code for both increments equal to 1 */

   /*        clean-up loop */

L20:
   m = n % 4;
   if (m == 0) {
      goto L40;
   }
   i__1 = m;
   for (i__ = 1; i__ <= i__1; ++i__) {
      dy[i__] += da * dx[i__];
      /* L30: */
   }
   if (n < 4) {
      return 0;
   }
L40:
   mp1 = m + 1;
   i__1 = n;
   for (i__ = mp1; i__ <= i__1; i__ += 4) {
      dy[i__] += da * dx[i__];
      dy[i__ + 1] += da * dx[i__ + 1];
      dy[i__ + 2] += da * dx[i__ + 2];
      dy[i__ + 3] += da * dx[i__ + 3];
      /* L50: */
   }
   return 0;
} /* daxpy_ */

} // namespace Minuit2

} // namespace ROOT
