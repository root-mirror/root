// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

/* mneig.F -- translated by f2c (version 20010320).
   You must link the resulting object file with the libraries:
        -lf2c -lm   (in that order)
*/

#include <cmath>

namespace ROOT {

namespace Minuit2 {

int mneigen(double *a, unsigned int ndima, unsigned int n, unsigned int mits, double *work, double precis)
{
   // compute matrix eignevalues (transaltion from mneig.F of Minuit)

   /* System generated locals */
   unsigned int a_dim1, a_offset, i__1, i__2, i__3;
   double r__1, r__2;

   /* Local variables */
   double b, c__, f, h__;
   unsigned int i__, j, k, l, m = 0;
   double r__, s;
   unsigned int i0, i1, j1, m1, n1;
   double hh, gl, pr, pt;

   /*          PRECIS is the machine precision EPSMAC */
   /* Parameter adjustments */
   a_dim1 = ndima;
   a_offset = 1 + a_dim1 * 1;
   a -= a_offset;
   --work;

   /* Function Body */
   int ifault = 1;

   i__ = n;
   i__1 = n;
   for (i1 = 2; i1 <= i__1; ++i1) {
      l = i__ - 2;
      f = a[i__ + (i__ - 1) * a_dim1];
      gl = (double)0.;

      if (l < 1) {
         goto L25;
      }

      i__2 = l;
      for (k = 1; k <= i__2; ++k) {
         /* Computing 2nd power */
         r__1 = a[i__ + k * a_dim1];
         gl += r__1 * r__1;
      }
   L25:
      /* Computing 2nd power */
      r__1 = f;
      h__ = gl + r__1 * r__1;

      if (gl > (double)1e-35) {
         goto L30;
      }

      work[i__] = (double)0.;
      work[n + i__] = f;
      goto L65;
   L30:
      ++l;

      gl = std::sqrt(h__);

      if (f >= (double)0.) {
         gl = -gl;
      }

      work[n + i__] = gl;
      h__ -= f * gl;
      a[i__ + (i__ - 1) * a_dim1] = f - gl;
      f = (double)0.;
      i__2 = l;
      for (j = 1; j <= i__2; ++j) {
         a[j + i__ * a_dim1] = a[i__ + j * a_dim1] / h__;
         gl = (double)0.;
         i__3 = j;
         for (k = 1; k <= i__3; ++k) {
            gl += a[j + k * a_dim1] * a[i__ + k * a_dim1];
         }

         if (j >= l) {
            goto L47;
         }

         j1 = j + 1;
         i__3 = l;
         for (k = j1; k <= i__3; ++k) {
            gl += a[k + j * a_dim1] * a[i__ + k * a_dim1];
         }
      L47:
         work[n + j] = gl / h__;
         f += gl * a[j + i__ * a_dim1];
      }
      hh = f / (h__ + h__);
      i__2 = l;
      for (j = 1; j <= i__2; ++j) {
         f = a[i__ + j * a_dim1];
         gl = work[n + j] - hh * f;
         work[n + j] = gl;
         i__3 = j;
         for (k = 1; k <= i__3; ++k) {
            a[j + k * a_dim1] = a[j + k * a_dim1] - f * work[n + k] - gl * a[i__ + k * a_dim1];
         }
      }
      work[i__] = h__;
   L65:
      --i__;
   }
   work[1] = (double)0.;
   work[n + 1] = (double)0.;
   i__1 = n;
   for (i__ = 1; i__ <= i__1; ++i__) {
      l = i__ - 1;

      if (work[i__] == (double)0. || l == 0) {
         goto L100;
      }

      i__3 = l;
      for (j = 1; j <= i__3; ++j) {
         gl = (double)0.;
         i__2 = l;
         for (k = 1; k <= i__2; ++k) {
            gl += a[i__ + k * a_dim1] * a[k + j * a_dim1];
         }
         i__2 = l;
         for (k = 1; k <= i__2; ++k) {
            a[k + j * a_dim1] -= gl * a[k + i__ * a_dim1];
         }
      }
   L100:
      work[i__] = a[i__ + i__ * a_dim1];
      a[i__ + i__ * a_dim1] = (double)1.;

      if (l == 0) {
         goto L110;
      }

      i__2 = l;
      for (j = 1; j <= i__2; ++j) {
         a[i__ + j * a_dim1] = (double)0.;
         a[j + i__ * a_dim1] = (double)0.;
      }
   L110:;
   }

   n1 = n - 1;
   i__1 = n;
   for (i__ = 2; i__ <= i__1; ++i__) {
      i0 = n + i__ - 1;
      work[i0] = work[i0 + 1];
   }
   work[n + n] = (double)0.;
   b = (double)0.;
   f = (double)0.;
   i__1 = n;
   for (l = 1; l <= i__1; ++l) {
      j = 0;
      h__ = precis * ((r__1 = work[l], std::fabs(r__1)) + (r__2 = work[n + l], std::fabs(r__2)));

      if (b < h__) {
         b = h__;
      }

      i__2 = n;
      for (m1 = l; m1 <= i__2; ++m1) {
         m = m1;

         if ((r__1 = work[n + m], std::fabs(r__1)) <= b) {
            goto L150;
         }
      }

   L150:
      if (m == l) {
         goto L205;
      }

   L160:
      if (j == mits) {
         return ifault;
      }

      ++j;
      pt = (work[l + 1] - work[l]) / (work[n + l] * (double)2.);
      r__ = std::sqrt(pt * pt + (double)1.);
      pr = pt + r__;

      if (pt < (double)0.) {
         pr = pt - r__;
      }

      h__ = work[l] - work[n + l] / pr;
      i__2 = n;
      for (i__ = l; i__ <= i__2; ++i__) {
         work[i__] -= h__;
      }
      f += h__;
      pt = work[m];
      c__ = (double)1.;
      s = (double)0.;
      m1 = m - 1;
      i__ = m;
      i__2 = m1;
      for (i1 = l; i1 <= i__2; ++i1) {
         j = i__;
         --i__;
         gl = c__ * work[n + i__];
         h__ = c__ * pt;

         if (std::fabs(pt) >= (r__1 = work[n + i__], std::fabs(r__1))) {
            goto L180;
         }

         c__ = pt / work[n + i__];
         r__ = std::sqrt(c__ * c__ + (double)1.);
         work[n + j] = s * work[n + i__] * r__;
         s = (double)1. / r__;
         c__ /= r__;
         goto L190;
      L180:
         c__ = work[n + i__] / pt;
         r__ = std::sqrt(c__ * c__ + (double)1.);
         work[n + j] = s * pt * r__;
         s = c__ / r__;
         c__ = (double)1. / r__;
      L190:
         pt = c__ * work[i__] - s * gl;
         work[j] = h__ + s * (c__ * gl + s * work[i__]);
         i__3 = n;
         for (k = 1; k <= i__3; ++k) {
            h__ = a[k + j * a_dim1];
            a[k + j * a_dim1] = s * a[k + i__ * a_dim1] + c__ * h__;
            a[k + i__ * a_dim1] = c__ * a[k + i__ * a_dim1] - s * h__;
         }
      }
      work[n + l] = s * pt;
      work[l] = c__ * pt;

      if ((r__1 = work[n + l], std::fabs(r__1)) > b) {
         goto L160;
      }

   L205:
      work[l] += f;
   }
   i__1 = n1;
   for (i__ = 1; i__ <= i__1; ++i__) {
      k = i__;
      pt = work[i__];
      i1 = i__ + 1;
      i__3 = n;
      for (j = i1; j <= i__3; ++j) {

         if (work[j] >= pt) {
            goto L220;
         }

         k = j;
         pt = work[j];
      L220:;
      }

      if (k == i__) {
         goto L240;
      }

      work[k] = work[i__];
      work[i__] = pt;
      i__3 = n;
      for (j = 1; j <= i__3; ++j) {
         pt = a[j + i__ * a_dim1];
         a[j + i__ * a_dim1] = a[j + k * a_dim1];
         a[j + k * a_dim1] = pt;
      }
   L240:;
   }
   ifault = 0;

   return ifault;
} /* mneig_ */

} // namespace Minuit2

} // namespace ROOT
