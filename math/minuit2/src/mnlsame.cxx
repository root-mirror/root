// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

/* lsame.f -- translated by f2c (version 20010320).
   You must link the resulting object file with the libraries:
   -lf2c -lm   (in that order)
*/

#include <cstring>

namespace ROOT {

namespace Minuit2 {

bool mnlsame(const char *ca, const char *cb)
{
   /* System generated locals */
   bool ret_val = false;

   /* Local variables */
   //     integer inta, intb, zcode;

   /*  -- LAPACK auxiliary routine (version 2.0) -- */
   /*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
   /*     Courant Institute, Argonne National Lab, and Rice University */
   /*     January 31, 1994 */

   /*     .. Scalar Arguments .. */
   /*     .. */

   /*  Purpose */
   /*  ======= */

   /*  LSAME returns .TRUE. if CA is the same letter as CB regardless of */
   /*  case. */

   /*  Arguments */
   /*  ========= */

   /*  CA      (input) CHARACTER*1 */
   /*  CB      (input) CHARACTER*1 */
   /*          CA and CB specify the single characters to be compared. */

   /* ===================================================================== */

   /*     .. Intrinsic Functions .. */
   /*     .. */
   /*     .. Local Scalars .. */
   /*     .. */
   /*     .. Executable Statements .. */

   /*     Test if the characters are equal */

   int comp = std::strcmp(ca, cb);
   if (comp == 0)
      ret_val = true;

   return ret_val;
} /* lsame_ */

} // namespace Minuit2

} // namespace ROOT
