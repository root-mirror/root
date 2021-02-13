// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

/* dasum.f -- translated by f2c (version 20010320).
   You must link the resulting object file with the libraries:
      -lf2c -lm   (in that order)
*/

#include <cmath>

namespace ROOT {

namespace Minuit2 {

double mndasum(unsigned int n, const double *dx, int incx)
{
   /* System generated locals */
   int i__1, i__2;
   double ret_val, d__1, d__2, d__3, d__4, d__5, d__6;

   /* Local variables */
   int i__, m;
   double dtemp;
   int nincx, mp1;

   /*     takes the sum of the absolute values. */
   /*     jack dongarra, linpack, 3/11/78. */
   /*     modified 3/93 to return if incx .le. 0. */
   /*     modified 12/3/93, array(1) declarations changed to array(*) */

   /* Parameter adjustments */
   --dx;

   /* Function Body */
   ret_val = 0.;
   dtemp = 0.;
   if (n <= 0 || incx <= 0) {
      return ret_val;
   }
   if (incx == 1) {
      goto L20;
   }

   /*        code for increment not equal to 1 */

   nincx = n * incx;
   i__1 = nincx;
   i__2 = incx;
   for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
      dtemp += (d__1 = dx[i__], std::fabs(d__1));
      /* L10: */
   }
   ret_val = dtemp;
   return ret_val;

   /*        code for increment equal to 1 */

   /*        clean-up loop */

L20:
   m = n % 6;
   if (m == 0) {
      goto L40;
   }
   i__2 = m;
   for (i__ = 1; i__ <= i__2; ++i__) {
      dtemp += (d__1 = dx[i__], std::fabs(d__1));
      /* L30: */
   }
   if (n < 6) {
      goto L60;
   }
L40:
   mp1 = m + 1;
   i__2 = n;
   for (i__ = mp1; i__ <= i__2; i__ += 6) {
      dtemp = dtemp + (d__1 = dx[i__], std::fabs(d__1)) + (d__2 = dx[i__ + 1], std::fabs(d__2)) +
              (d__3 = dx[i__ + 2], std::fabs(d__3)) + (d__4 = dx[i__ + 3], std::fabs(d__4)) +
              (d__5 = dx[i__ + 4], std::fabs(d__5)) + (d__6 = dx[i__ + 5], std::fabs(d__6));
      /* L50: */
   }
L60:
   ret_val = dtemp;
   return ret_val;
} /* dasum_ */

} // namespace Minuit2

} // namespace ROOT
