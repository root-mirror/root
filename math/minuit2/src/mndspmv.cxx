// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

/* dspmv.f -- translated by f2c (version 20010320).
   You must link the resulting object file with the libraries:
   -lf2c -lm   (in that order)
*/

namespace ROOT {

namespace Minuit2 {

bool mnlsame(const char *, const char *);
int mnxerbla(const char *, int);

int Mndspmv(const char *uplo, unsigned int n, double alpha, const double *ap, const double *x, int incx, double beta,
            double *y, int incy)
{
   /* System generated locals */
   int i__1, i__2;

   /* Local variables */
   int info;
   double temp1, temp2;
   int i__, j, k;
   int kk, ix, iy, jx, jy, kx, ky;

   /*     .. Scalar Arguments .. */
   /*     .. Array Arguments .. */
   /*     .. */

   /*  Purpose */
   /*  ======= */

   /*  DSPMV  performs the matrix-vector operation */

   /*     y := alpha*A*x + beta*y, */

   /*  where alpha and beta are scalars, x and y are n element vectors and */
   /*  A is an n by n symmetric matrix, supplied in packed form. */

   /*  Parameters */
   /*  ========== */

   /*  UPLO   - CHARACTER*1. */
   /*           On entry, UPLO specifies whether the Upper or Lower */
   /*           triangular part of the matrix A is supplied in the packed */
   /*           array AP as follows: */

   /*              UPLO = 'U' or 'u'   The Upper triangular part of A is */
   /*                                  supplied in AP. */

   /*              UPLO = 'L' or 'l'   The Lower triangular part of A is */
   /*                                  supplied in AP. */

   /*           Unchanged on exit. */

   /*  N      - INTEGER. */
   /*           On entry, N specifies the order of the matrix A. */
   /*           N must be at least zero. */
   /*           Unchanged on exit. */

   /*  ALPHA  - DOUBLE PRECISION. */
   /*           On entry, ALPHA specifies the scalar alpha. */
   /*           Unchanged on exit. */

   /*  AP     - DOUBLE PRECISION array of DIMENSION at least */
   /*           ( ( n*( n + 1 ) )/2 ). */
   /*           Before entry with UPLO = 'U' or 'u', the array AP must */
   /*           contain the Upper triangular part of the symmetric matrix */
   /*           packed sequentially, column by column, so that AP( 1 ) */
   /*           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 ) */
   /*           and a( 2, 2 ) respectively, and so on. */
   /*           Before entry with UPLO = 'L' or 'l', the array AP must */
   /*           contain the Lower triangular part of the symmetric matrix */
   /*           packed sequentially, column by column, so that AP( 1 ) */
   /*           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 ) */
   /*           and a( 3, 1 ) respectively, and so on. */
   /*           Unchanged on exit. */

   /*  X      - DOUBLE PRECISION array of dimension at least */
   /*           ( 1 + ( n - 1 )*abs( INCX ) ). */
   /*           Before entry, the incremented array X must contain the n */
   /*           element vector x. */
   /*           Unchanged on exit. */

   /*  INCX   - INTEGER. */
   /*           On entry, INCX specifies the increment for the Elements of */
   /*           X. INCX must not be zero. */
   /*           Unchanged on exit. */

   /*  BETA   - DOUBLE PRECISION. */
   /*           On entry, BETA specifies the scalar beta. When BETA is */
   /*           supplied as zero then Y need not be set on input. */
   /*           Unchanged on exit. */

   /*  Y      - DOUBLE PRECISION array of dimension at least */
   /*           ( 1 + ( n - 1 )*abs( INCY ) ). */
   /*           Before entry, the incremented array Y must contain the n */
   /*           element vector y. On exit, Y is overwritten by the updated */
   /*           vector y. */

   /*  INCY   - INTEGER. */
   /*           On entry, INCY specifies the increment for the Elements of */
   /*           Y. INCY must not be zero. */
   /*           Unchanged on exit. */

   /*  Level 2 Blas routine. */

   /*  -- Written on 22-October-1986. */
   /*     Jack Dongarra, Argonne National Lab. */
   /*     Jeremy Du Croz, Nag Central Office. */
   /*     Sven Hammarling, Nag Central Office. */
   /*     Richard Hanson, Sandia National Labs. */

   /*     .. Parameters .. */
   /*     .. Local Scalars .. */
   /*     .. External Functions .. */
   /*     .. External Subroutines .. */
   /*     .. */
   /*     .. Executable Statements .. */

   /*     Test the input parameters. */

   /* Parameter adjustments */
   --y;
   --x;
   --ap;

   /* Function Body */
   info = 0;
   if (!mnlsame(uplo, "U") && !mnlsame(uplo, "L")) {
      info = 1;
   }
   //     else if (n < 0) {
   //       info = 2;
   //     }
   else if (incx == 0) {
      info = 6;
   } else if (incy == 0) {
      info = 9;
   }
   if (info != 0) {
      mnxerbla("DSPMV ", info);
      return 0;
   }

   /*     Quick return if possible. */

   if ((n == 0) || (alpha == 0. && beta == 1.)) {
      return 0;
   }

   /*     Set up the start points in  X  and  Y. */

   if (incx > 0) {
      kx = 1;
   } else {
      kx = 1 - (n - 1) * incx;
   }
   if (incy > 0) {
      ky = 1;
   } else {
      ky = 1 - (n - 1) * incy;
   }

   /*     Start the operations. In this version the Elements of the array AP */
   /*     are accessed sequentially with one pass through AP. */

   /*     First form  y := beta*y. */

   if (beta != 1.) {
      if (incy == 1) {
         if (beta == 0.) {
            i__1 = n;
            for (i__ = 1; i__ <= i__1; ++i__) {
               y[i__] = 0.;
               /* L10: */
            }
         } else {
            i__1 = n;
            for (i__ = 1; i__ <= i__1; ++i__) {
               y[i__] = beta * y[i__];
               /* L20: */
            }
         }
      } else {
         iy = ky;
         if (beta == 0.) {
            i__1 = n;
            for (i__ = 1; i__ <= i__1; ++i__) {
               y[iy] = 0.;
               iy += incy;
               /* L30: */
            }
         } else {
            i__1 = n;
            for (i__ = 1; i__ <= i__1; ++i__) {
               y[iy] = beta * y[iy];
               iy += incy;
               /* L40: */
            }
         }
      }
   }
   if (alpha == 0.) {
      return 0;
   }
   kk = 1;
   if (mnlsame(uplo, "U")) {

      /*        Form  y  when AP contains the Upper triangle. */

      if (incx == 1 && incy == 1) {
         i__1 = n;
         for (j = 1; j <= i__1; ++j) {
            temp1 = alpha * x[j];
            temp2 = 0.;
            k = kk;
            i__2 = j - 1;
            for (i__ = 1; i__ <= i__2; ++i__) {
               y[i__] += temp1 * ap[k];
               temp2 += ap[k] * x[i__];
               ++k;
               /* L50: */
            }
            y[j] = y[j] + temp1 * ap[kk + j - 1] + alpha * temp2;
            kk += j;
            /* L60: */
         }
      } else {
         jx = kx;
         jy = ky;
         i__1 = n;
         for (j = 1; j <= i__1; ++j) {
            temp1 = alpha * x[jx];
            temp2 = 0.;
            ix = kx;
            iy = ky;
            i__2 = kk + j - 2;
            for (k = 0; k <= i__2 - kk; ++k) {
               y[iy] += temp1 * ap[k + kk];
               temp2 += ap[k + kk] * x[ix];
               ix += incx;
               iy += incy;
               /* L70: */
            }
            y[jy] = y[jy] + temp1 * ap[kk + j - 1] + alpha * temp2;
            jx += incx;
            jy += incy;
            kk += j;
            /* L80: */
         }
      }
   } else {

      /*        Form  y  when AP contains the Lower triangle. */

      if (incx == 1 && incy == 1) {
         i__1 = n;
         for (j = 1; j <= i__1; ++j) {
            temp1 = alpha * x[j];
            temp2 = 0.;
            y[j] += temp1 * ap[kk];
            k = kk + 1;
            i__2 = n;
            for (i__ = j + 1; i__ <= i__2; ++i__) {
               y[i__] += temp1 * ap[k];
               temp2 += ap[k] * x[i__];
               ++k;
               /* L90: */
            }
            y[j] += alpha * temp2;
            kk += n - j + 1;
            /* L100: */
         }
      } else {
         jx = kx;
         jy = ky;
         i__1 = n;
         for (j = 1; j <= i__1; ++j) {
            temp1 = alpha * x[jx];
            temp2 = 0.;
            y[jy] += temp1 * ap[kk];
            ix = jx;
            iy = jy;
            i__2 = kk + n - j;
            for (k = kk + 1; k <= i__2; ++k) {
               ix += incx;
               iy += incy;
               y[iy] += temp1 * ap[k];
               temp2 += ap[k] * x[ix];
               /* L110: */
            }
            y[jy] += alpha * temp2;
            jx += incx;
            jy += incy;
            kk += n - j + 1;
            /* L120: */
         }
      }
   }

   return 0;

   /*     End of DSPMV . */

} /* dspmv_ */

} // namespace Minuit2

} // namespace ROOT
