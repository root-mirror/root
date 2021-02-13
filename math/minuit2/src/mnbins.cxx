// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#include <cmath>

namespace ROOT {

namespace Minuit2 {

void mnbins(double a1, double a2, int naa, double &bl, double &bh, int &nb, double &bwid)
{

   //*-*-*-*-*-*-*-*-*-*-*Compute reasonable histogram intervals*-*-*-*-*-*-*-*-*
   //*-*                  ======================================
   //*-*        Function TO DETERMINE REASONABLE HISTOGRAM INTERVALS
   //*-*        GIVEN ABSOLUTE UPPER AND LOWER BOUNDS  A1 AND A2
   //*-*        AND DESIRED MAXIMUM NUMBER OF BINS NAA
   //*-*        PROGRAM MAKES REASONABLE BINNING FROM BL TO BH OF WIDTH BWID
   //*-*        F. JAMES,   AUGUST, 1974 , stolen for Minuit, 1988
   //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

   /* Local variables */
   double awid, ah, al, sigfig, sigrnd, alb;
   int kwid, lwid, na = 0, log_;

   al = a1 < a2 ? a1 : a2;
   //     al = std::min(a1,a2);
   //     ah = std::max(a1,a2);
   ah = a1 > a2 ? a1 : a2;
   if (al == ah)
      ah = al + 1;

   //*-*-       IF NAA .EQ. -1 , PROGRAM USES BWID INPUT FROM CALLING ROUTINE
   if (naa == -1)
      goto L150;
L10:
   na = naa - 1;
   if (na < 1)
      na = 1;

   //*-*-        GET NOMINAL BIN WIDTH IN EXPON FORM
L20:
   awid = (ah - al) / double(na);
   log_ = int(log10(awid));
   if (awid <= 1)
      --log_;
   sigfig = awid * pow(10.0, -log_);
   //*-*-       ROUND MANTISSA UP TO 2, 2.5, 5, OR 10
   if (sigfig > 2)
      goto L40;
   sigrnd = 2;
   goto L100;
L40:
   if (sigfig > 2.5)
      goto L50;
   sigrnd = 2.5;
   goto L100;
L50:
   if (sigfig > 5)
      goto L60;
   sigrnd = 5;
   goto L100;
L60:
   sigrnd = 1;
   ++log_;
L100:
   bwid = sigrnd * pow(10.0, log_);
   goto L200;
   //*-*-       GET NEW BOUNDS FROM NEW WIDTH BWID
L150:
   if (bwid <= 0)
      goto L10;
L200:
   alb = al / bwid;
   lwid = int(alb);
   if (alb < 0)
      --lwid;
   bl = bwid * double(lwid);
   alb = ah / bwid + 1;
   kwid = int(alb);
   if (alb < 0)
      --kwid;
   bh = bwid * double(kwid);
   nb = kwid - lwid;
   if (naa > 5)
      goto L240;
   if (naa == -1)
      return;
   //*-*-        REQUEST FOR ONE BIN IS DIFFICULT CASE
   if (naa > 1 || nb == 1)
      return;
   bwid *= 2;
   nb = 1;
   return;
L240:
   if (nb << 1 != naa)
      return;
   ++na;
   goto L20;
} /* mnbins_ */

} // namespace Minuit2

} // namespace ROOT
