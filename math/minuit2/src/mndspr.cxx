// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

/* dspr.f -- translated by f2c (version 20010320).
   You must link the resulting object file with the libraries:
   -lf2c -lm   (in that order)
*/

namespace ROOT {

namespace Minuit2 {

bool mnlsame(const char *, const char *);
int mnxerbla(const char *, int);

int mndspr(const char *uplo, unsigned int n, double alpha, const double *x, int incx, double *ap)
{
   /* System generated locals */
   int i__1, i__2;

   /* Local variables */
   int info;
   double temp;
   int i__, j, k;
   int kk, ix, jx, kx = 0;

   /*     .. Scalar Arguments .. */
   /*     .. Array Arguments .. */
   /*     .. */

   /*  Purpose */
   /*  ======= */

   /*  DSPR    performs the symmetric rank 1 operation */

   /*     A := alpha*x*x' + A, */

   /*  where alpha is a real scalar, x is an n element vector and A is an */
   /*  n by n symmetric matrix, supplied in packed form. */

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

   /*  X      - DOUBLE PRECISION array of dimension at least */
   /*           ( 1 + ( n - 1 )*abs( INCX ) ). */
   /*           Before entry, the incremented array X must contain the n */
   /*           element vector x. */
   /*           Unchanged on exit. */

   /*  INCX   - INTEGER. */
   /*           On entry, INCX specifies the increment for the Elements of */
   /*           X. INCX must not be zero. */
   /*           Unchanged on exit. */

   /*  AP     - DOUBLE PRECISION array of DIMENSION at least */
   /*           ( ( n*( n + 1 ) )/2 ). */
   /*           Before entry with  UPLO = 'U' or 'u', the array AP must */
   /*           contain the Upper triangular part of the symmetric matrix */
   /*           packed sequentially, column by column, so that AP( 1 ) */
   /*           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 ) */
   /*           and a( 2, 2 ) respectively, and so on. On exit, the array */
   /*           AP is overwritten by the Upper triangular part of the */
   /*           updated matrix. */
   /*           Before entry with UPLO = 'L' or 'l', the array AP must */
   /*           contain the Lower triangular part of the symmetric matrix */
   /*           packed sequentially, column by column, so that AP( 1 ) */
   /*           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 ) */
   /*           and a( 3, 1 ) respectively, and so on. On exit, the array */
   /*           AP is overwritten by the Lower triangular part of the */
   /*           updated matrix. */

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
   --ap;
   --x;

   /* Function Body */
   info = 0;
   if (!mnlsame(uplo, "U") && !mnlsame(uplo, "L")) {
      info = 1;
   }
   //     else if (n < 0) {
   //    info = 2;
   //     }
   else if (incx == 0) {
      info = 5;
   }
   if (info != 0) {
      mnxerbla("DSPR  ", info);
      return 0;
   }

   /*     Quick return if possible. */

   if (n == 0 || alpha == 0.) {
      return 0;
   }

   /*     Set the start point in X if the increment is not unity. */

   if (incx <= 0) {
      kx = 1 - (n - 1) * incx;
   } else if (incx != 1) {
      kx = 1;
   }

   /*     Start the operations. In this version the Elements of the array AP */
   /*     are accessed sequentially with one pass through AP. */

   kk = 1;
   if (mnlsame(uplo, "U")) {

      /*        Form  A  when Upper triangle is stored in AP. */

      if (incx == 1) {
         i__1 = n;
         for (j = 1; j <= i__1; ++j) {
            if (x[j] != 0.) {
               temp = alpha * x[j];
               k = kk;
               i__2 = j;
               for (i__ = 1; i__ <= i__2; ++i__) {
                  ap[k] += x[i__] * temp;
                  ++k;
                  /* L10: */
               }
            }
            kk += j;
            /* L20: */
         }
      } else {
         jx = kx;
         i__1 = n;
         for (j = 1; j <= i__1; ++j) {
            if (x[jx] != 0.) {
               temp = alpha * x[jx];
               ix = kx;
               i__2 = kk + j - 1;
               for (k = kk; k <= i__2; ++k) {
                  ap[k] += x[ix] * temp;
                  ix += incx;
                  /* L30: */
               }
            }
            jx += incx;
            kk += j;
            /* L40: */
         }
      }
   } else {

      /*        Form  A  when Lower triangle is stored in AP. */

      if (incx == 1) {
         i__1 = n;
         for (j = 1; j <= i__1; ++j) {
            if (x[j] != 0.) {
               temp = alpha * x[j];
               k = kk;
               i__2 = n;
               for (i__ = j; i__ <= i__2; ++i__) {
                  ap[k] += x[i__] * temp;
                  ++k;
                  /* L50: */
               }
            }
            kk = kk + n - j + 1;
            /* L60: */
         }
      } else {
         jx = kx;
         i__1 = n;
         for (j = 1; j <= i__1; ++j) {
            if (x[jx] != 0.) {
               temp = alpha * x[jx];
               ix = jx;
               i__2 = kk + n - j;
               for (k = kk; k <= i__2; ++k) {
                  ap[k] += x[ix] * temp;
                  ix += incx;
                  /* L70: */
               }
            }
            jx += incx;
            kk = kk + n - j + 1;
            /* L80: */
         }
      }
   }

   return 0;

   /*     End of DSPR  . */

} /* dspr_ */

} // namespace Minuit2

} // namespace ROOT
