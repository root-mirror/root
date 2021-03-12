// @(#)root/mathmore:$Id$
// Authors: L. Moneta, A. Zsenei   08/2005

/* poly/zsolve_quartic.c
 *
 * Copyright (C) 2003 CERN and K.S. K\"{o}lbig
 *
 * Converted to C and implemented into the GSL Library - Sept. 2003
 * by Andrew W. Steiner and Andy Buckley
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

/* zsolve_quartic.c - finds the complex roots of
 *  x^4 + a x^3 + b x^2 + c x + d = 0
 */

#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_poly.h>

#define SWAP(a,b) do { gsl_complex tmp = b ; b = a ; a = tmp ; } while(0)

int
gsl_poly_complex_solve_quartic (double a, double b, double c, double d,
                                gsl_complex * z0, gsl_complex * z1,
                                gsl_complex * z2, gsl_complex * z3)
{
  gsl_complex i, zarr[4], w1, w2, w3;
  double r4 = 1.0 / 4.0;
  double q2 = 1.0 / 2.0, q4 = 1.0 / 4.0, q8 = 1.0 / 8.0;
  double q1 = 3.0 / 8.0, q3 = 3.0 / 16.0;
  double u[3], v[3], v1, v2, disc;
  double aa, pp, qq, rr, rc, sc, tc, q, h;
  int k1 = 0, k2 = 0, mt;

  GSL_SET_COMPLEX (&i, 0.0, 1.0);
  GSL_SET_COMPLEX (&zarr[0], 0.0, 0.0);
  GSL_SET_COMPLEX (&zarr[1], 0.0, 0.0);
  GSL_SET_COMPLEX (&zarr[2], 0.0, 0.0);
  GSL_SET_COMPLEX (&zarr[3], 0.0, 0.0);
  GSL_SET_COMPLEX (&w1, 0.0, 0.0);
  GSL_SET_COMPLEX (&w2, 0.0, 0.0);
  GSL_SET_COMPLEX (&w3, 0.0, 0.0);

  /* Deal easily with the cases where the quartic is degenerate. The
   * ordering of solutions is done explicitly. */
  if (0 == b && 0 == c)
    {
      if (0 == d)
        {
          if (a > 0)
            {
              GSL_SET_COMPLEX (z0, -a, 0.0);
              GSL_SET_COMPLEX (z1, 0.0, 0.0);
              GSL_SET_COMPLEX (z2, 0.0, 0.0);
              GSL_SET_COMPLEX (z3, 0.0, 0.0);
            }
          else
            {
              GSL_SET_COMPLEX (z0, 0.0, 0.0);
              GSL_SET_COMPLEX (z1, 0.0, 0.0);
              GSL_SET_COMPLEX (z2, 0.0, 0.0);
              GSL_SET_COMPLEX (z3, -a, 0.0);
            }
          return 4;
        }
      else if (0 == a)
        {
          if (d > 0)
            {
              double sqrt_d = sqrt (d);
              gsl_complex i_sqrt_d = gsl_complex_mul_real (i, sqrt_d);
              gsl_complex minus_i = gsl_complex_conjugate (i);
              *z3 = gsl_complex_sqrt (i_sqrt_d);
              *z2 = gsl_complex_mul (minus_i, *z3);
              *z1 = gsl_complex_negative (*z2);
              *z0 = gsl_complex_negative (*z3);
            }
          else
            {
              double sqrt_abs_d = sqrt (-d);
              *z3 = gsl_complex_sqrt_real (sqrt_abs_d);
              *z2 = gsl_complex_mul (i, *z3);
              *z1 = gsl_complex_negative (*z2);
              *z0 = gsl_complex_negative (*z3);
            }
          return 4;
        }
    }

  if (0.0 == c && 0.0 == d)
    {
      disc = (a * a - 4.0 * b);
      if (disc < 0.0)
        {
          mt = 3;
        }
      else
        {
          mt = 1;
        }
      *z0 = zarr[0];
      *z1 = zarr[0];
      gsl_poly_complex_solve_quadratic (1.0, a, b, z2, z3);
    }
  else
    {
      /* For non-degenerate solutions, proceed by constructing and
       * solving the resolvent cubic */
      aa = a * a;
      pp = b - q1 * aa;
      qq = c - q2 * a * (b - q4 * aa);
      rr = d - q4 * (a * c - q4 * aa * (b - q3 * aa));
      rc = q2 * pp;
      sc = q4 * (q4 * pp * pp - rr);
      tc = -(q8 * qq * q8 * qq);

      /* This code solves the resolvent cubic in a convenient fashion
       * for this implementation of the quartic. If there are three real
       * roots, then they are placed directly into u[].  If two are
       * complex, then the real root is put into u[0] and the real
       * and imaginary part of the complex roots are placed into
       * u[1] and u[2], respectively. Additionally, this
       * calculates the discriminant of the cubic and puts it into the
       * variable disc. */
      {
        double qcub = (rc * rc - 3 * sc);
        double rcub = (2 * rc * rc * rc - 9 * rc * sc + 27 * tc);

        double Q = qcub / 9;
        double R = rcub / 54;

        double Q3 = Q * Q * Q;
        double R2 = R * R;

        disc = R2 - Q3;

//       more numerical problems with this calculation of disc
//LM - use this since it fixes case reported in issue #6900 when disc is approx 0
//      (needs to add additional test cases where disc is ~ 0)
        double CR2 = 729 * rcub * rcub;
        double CQ3 = 2916 * qcub * qcub * qcub;
        disc = (CR2 - CQ3) / 2125764.0;


        if (0 == R && 0 == Q)
          {
            u[0] = -rc / 3;
            u[1] = -rc / 3;
            u[2] = -rc / 3;
          }
        else if (disc == 0)
          {
            double sqrtQ = sqrt (Q);
            if (R > 0)
              {
                u[0] = -2 * sqrtQ - rc / 3;
                u[1] = sqrtQ - rc / 3;
                u[2] = sqrtQ - rc / 3;
              }
            else
              {
                u[0] = -sqrtQ - rc / 3;
                u[1] = -sqrtQ - rc / 3;
                u[2] = 2 * sqrtQ - rc / 3;
              }
          }
        else if ( R2 < Q3)
          {
            double sqrtQ = sqrt (Q);
            double sqrtQ3 = sqrtQ * sqrtQ * sqrtQ;
            double ctheta = R / sqrtQ3;
            double theta = 0;
            // protect against numerical error can make this larger than one
            if ( fabs(ctheta) < 1.0 )
               theta = acos( ctheta);
            else if ( ctheta <= -1.0)
               theta = M_PI;

            double norm = -2 * sqrtQ;

            u[0] = norm * cos (theta / 3) - rc / 3;
            u[1] = norm * cos ((theta + 2.0 * M_PI) / 3) - rc / 3;
            u[2] = norm * cos ((theta - 2.0 * M_PI) / 3) - rc / 3;
          }
        else
          {
            double sgnR = (R >= 0 ? 1 : -1);
            double modR = fabs (R);
            double sqrt_disc = sqrt (disc);
            double A = -sgnR * pow (modR + sqrt_disc, 1.0 / 3.0);
            double B = Q / A;

            double mod_diffAB = fabs (A - B);
            u[0] = A + B - rc / 3;
            u[1] = -0.5 * (A + B) - rc / 3;
            u[2] = -(sqrt (3.0) / 2.0) * mod_diffAB;
          }
      }
      /* End of solution to resolvent cubic */

      /* Combine the square roots of the roots of the cubic
       * resolvent appropriately. Also, calculate 'mt' which
       * designates the nature of the roots:
       * mt=1 : 4 real roots
       * mt=2 : 0 real roots
       * mt=3 : 2 real roots
       */
      // when disc == 0  2 roots are identicals
      if (0 >= disc)
        {
          mt = 2;
          v[0] = fabs (u[0]);
          v[1] = fabs (u[1]);
          v[2] = fabs (u[2]);

          v1 = GSL_MAX (GSL_MAX (v[0], v[1]), v[2]);
          if (v1 == v[0])
            {
              k1 = 0;
              v2 = GSL_MAX (v[1], v[2]);
            }
          else if (v1 == v[1])
            {
              k1 = 1;
              v2 = GSL_MAX (v[0], v[2]);
            }
          else
            {
              k1 = 2;
              v2 = GSL_MAX (v[0], v[1]);
            }

          if (v2 == v[0])
            {
              k2 = 0;
            }
          else if (v2 == v[1])
            {
              k2 = 1;
            }
          else
            {
              k2 = 2;
            }
          w1 = gsl_complex_sqrt_real (u[k1]);
          w2 = gsl_complex_sqrt_real (u[k2]);
        }
      else
        {
          mt = 3;
          GSL_SET_COMPLEX (&w1, u[1], u[2]);
          GSL_SET_COMPLEX (&w2, u[1], -u[2]);
          w1 = gsl_complex_sqrt (w1);
          w2 = gsl_complex_sqrt (w2);
        }
      /* Solve the quadratic in order to obtain the roots
       * to the quartic */
      q = qq;
      gsl_complex prod_w = gsl_complex_mul (w1, w2);
      //gsl_complex mod_prod_w = gsl_complex_abs (prod_w);
      /*
       Changed from gsl_complex to double in order to make it compile.
      */
      double mod_prod_w = gsl_complex_abs (prod_w);
      if (0.0 != mod_prod_w)
        {
          gsl_complex inv_prod_w = gsl_complex_inverse (prod_w);
          w3 = gsl_complex_mul_real (inv_prod_w, -q / 8.0);
        }

      h = r4 * a;
      gsl_complex sum_w12 = gsl_complex_add (w1, w2);
      gsl_complex neg_sum_w12 = gsl_complex_negative (sum_w12);
      gsl_complex sum_w123 = gsl_complex_add (sum_w12, w3);
      gsl_complex neg_sum_w123 = gsl_complex_add (neg_sum_w12, w3);

      gsl_complex diff_w12 = gsl_complex_sub (w2, w1);
      gsl_complex neg_diff_w12 = gsl_complex_negative (diff_w12);
      gsl_complex diff_w123 = gsl_complex_sub (diff_w12, w3);
      gsl_complex neg_diff_w123 = gsl_complex_sub (neg_diff_w12, w3);

      zarr[0] = gsl_complex_add_real (sum_w123, -h);
      zarr[1] = gsl_complex_add_real (neg_sum_w123, -h);
      zarr[2] = gsl_complex_add_real (diff_w123, -h);
      zarr[3] = gsl_complex_add_real (neg_diff_w123, -h);

      /* Arrange the roots into the variables z0, z1, z2, z3 */
      if (2 == mt)
        {
          if (u[k1] >= 0 && u[k2] >= 0)
            {
              mt = 1;
              GSL_SET_COMPLEX (z0, GSL_REAL (zarr[0]), 0.0);
              GSL_SET_COMPLEX (z1, GSL_REAL (zarr[1]), 0.0);
              GSL_SET_COMPLEX (z2, GSL_REAL (zarr[2]), 0.0);
              GSL_SET_COMPLEX (z3, GSL_REAL (zarr[3]), 0.0);
            }
          else if (u[k1] >= 0 && u[k2] < 0)
            {
              *z0 = zarr[0];
              *z1 = zarr[3];
              *z2 = zarr[2];
              *z3 = zarr[1];
            }
          else if (u[k1] < 0 && u[k2] >= 0)
            {
              *z0 = zarr[0];
              *z1 = zarr[2];
              *z2 = zarr[3];
              *z3 = zarr[1];
            }
          else if (u[k1] < 0 && u[k2] < 0)
            {
              *z0 = zarr[0];
              *z1 = zarr[1];
              *z2 = zarr[3];
              *z3 = zarr[2];
            }
        }
      else if (3 == mt)
        {
          GSL_SET_COMPLEX (z0, GSL_REAL (zarr[0]), 0.0);
          GSL_SET_COMPLEX (z1, GSL_REAL (zarr[1]), 0.0);
          *z2 = zarr[3];
          *z3 = zarr[2];
        }
    }

  /*
   * Sort the roots as usual: main sorting by ascending real part, secondary
   * sorting by ascending imaginary part
   */

  if (1 == mt)
    {
      /* Roots are all real, sort them by the real part */
      if (GSL_REAL (*z0) > GSL_REAL (*z1)) SWAP (*z0, *z1);
      if (GSL_REAL (*z0) > GSL_REAL (*z2)) SWAP (*z0, *z2);
      if (GSL_REAL (*z0) > GSL_REAL (*z3)) SWAP (*z0, *z3);

      if (GSL_REAL (*z1) > GSL_REAL (*z2)) SWAP (*z1, *z2);
      if (GSL_REAL (*z2) > GSL_REAL (*z3))
        {
          SWAP (*z2, *z3);
          if (GSL_REAL (*z1) > GSL_REAL (*z2)) SWAP (*z1, *z2);
        }
    }
  else if (2 == mt)
    {
      /* Roots are all complex. z0 and z1 are conjugates
       * and z2 and z3 are conjugates. Sort the real parts first */
      if (GSL_REAL (*z0) > GSL_REAL (*z2))
        {
          SWAP (*z0, *z2);
          SWAP (*z1, *z3);
        }
      /* Then sort by the imaginary parts */
      if (GSL_IMAG (*z0) > GSL_IMAG (*z1)) SWAP (*z0, *z1);
      if (GSL_IMAG (*z2) > GSL_IMAG (*z3)) SWAP (*z2, *z3);
    }
  else
    {
      /* 2 real roots. z2 and z3 are conjugates. */

      /* Swap complex roots */
      if (GSL_IMAG (*z2) > GSL_IMAG (*z3)) SWAP (*z2, *z3);

      /* Sort real parts */
      if (GSL_REAL (*z0) > GSL_REAL (*z1)) SWAP (*z0, *z1);
      if (GSL_REAL (*z1) > GSL_REAL (*z2))
        {
          if (GSL_REAL (*z0) > GSL_REAL (*z2))
            {
              SWAP (*z0, *z2);
              SWAP (*z1, *z3);
            }
          else
            {
              SWAP (*z1, *z2);
              SWAP (*z2, *z3);
            }
        }
    }

  return 4;
}
