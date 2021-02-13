// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

/* xerbla.f -- translated by f2c (version 20010320).
   You must link the resulting object file with the libraries:
   -lf2c -lm   (in that order)
*/

#include "Minuit2/MnConfig.h"
#include "Minuit2/MnPrint.h"

namespace ROOT {

namespace Minuit2 {

/* Table of constant values */

// static integer c__1 = 1;

int mnxerbla(const char *srname, int info)
{
   /* Format strings */
   //     static char fmt_9999[] = "(\002 ** On entry to \002,a6,\002 Parameter nu\// mber \002,i2,\002 had \002,\002an
   //     illegal Value\002)";

   /*  -- LAPACK auxiliary routine (version 3.0) -- */
   /*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
   /*     Courant Institute, Argonne National Lab, and Rice University */
   /*     September 30, 1994 */

   /*     .. Scalar Arguments .. */
   /*     .. */

   /*  Purpose */
   /*  ======= */

   /*  XERBLA  is an Error handler for the LAPACK routines. */
   /*  It is called by an LAPACK routine if an input Parameter has an */
   /*  invalid Value.  A message is printed and execution stops. */

   /*  Installers may consider modifying the STOP statement in order to */
   /*  call system-specific exception-handling facilities. */

   /*  Arguments */
   /*  ========= */

   /*  SRNAME  (input) CHARACTER*6 */
   /*          The Name of the routine which called XERBLA. */

   /*  INFO    (input) INTEGER */
   /*          The position of the invalid Parameter in the Parameter list */
   /*          of the calling routine. */

   /* ===================================================================== */

   /*     .. Executable Statements .. */

   MnPrint print("mnxerbla");
   print.Warn("On entry to", srname, "Parameter number", info, "had an illegal Value");

   /*     End of XERBLA */

   return 0;
} /* xerbla_ */

} // namespace Minuit2

} // namespace ROOT
