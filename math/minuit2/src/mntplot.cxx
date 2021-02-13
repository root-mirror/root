// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

/* mnplot.F -- translated by f2c (version 20010320).
   You must link the resulting object file with the libraries:
   -lf2c -lm   (in that order)
*/

#include <cmath>
#include <cstdio>
#include <cstring>

namespace ROOT {

namespace Minuit2 {

void mnbins(double, double, int, double &, double &, int &, double &);

void mnplot(double *xpt, double *ypt, char *chpt, int nxypt, int npagwd, int npagln)
{
   //*-*-*-*Plots points in array xypt onto one page with labelled axes*-*-*-*-*
   //*-*    ===========================================================
   //*-*        NXYPT is the number of points to be plotted
   //*-*        XPT(I) = x-coord. of ith point
   //*-*        YPT(I) = y-coord. of ith point
   //*-*        CHPT(I) = character to be plotted at this position
   //*-*        the input point arrays XPT, YPT, CHPT are destroyed.
   //*-*
   //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

   //      char cdot[]   = ".";
   //      char cslash[] = "/";

   /* Local variables */
   double xmin, ymin, xmax, ymax, savx, savy, yprt;
   double bwidx, bwidy, xbest, ybest, ax, ay, bx, by;
   double xvalus[12], any, dxx, dyy;
   int iten, i, j, k, maxnx, maxny, iquit, ni, linodd;
   int nxbest, nybest, km1, ibk, isp1, nx, ny, ks, ix;
   char ctemp[120];
   bool overpr;
   char cline[120];
   char chsav, chbest;

   /* Function Body */
   //*-*  Computing MIN
   maxnx = npagwd - 20 < 100 ? npagwd - 20 : 100;
   if (maxnx < 10)
      maxnx = 10;
   maxny = npagln;
   if (maxny < 10)
      maxny = 10;
   if (nxypt <= 1)
      return;
   xbest = xpt[0];
   ybest = ypt[0];
   chbest = chpt[0];
   //*-*-        order the points by decreasing y
   km1 = nxypt - 1;
   for (i = 1; i <= km1; ++i) {
      iquit = 0;
      ni = nxypt - i;
      for (j = 1; j <= ni; ++j) {
         if (ypt[j - 1] > ypt[j])
            continue;
         savx = xpt[j - 1];
         xpt[j - 1] = xpt[j];
         xpt[j] = savx;
         savy = ypt[j - 1];
         ypt[j - 1] = ypt[j];
         ypt[j] = savy;
         chsav = chpt[j - 1];
         chpt[j - 1] = chpt[j];
         chpt[j] = chsav;
         iquit = 1;
      }
      if (iquit == 0)
         break;
   }
   //*-*-        find extreme values
   xmax = xpt[0];
   xmin = xmax;
   for (i = 1; i <= nxypt; ++i) {
      if (xpt[i - 1] > xmax)
         xmax = xpt[i - 1];
      if (xpt[i - 1] < xmin)
         xmin = xpt[i - 1];
   }
   dxx = (xmax - xmin) * .001;
   xmax += dxx;
   xmin -= dxx;
   mnbins(xmin, xmax, maxnx, xmin, xmax, nx, bwidx);
   ymax = ypt[0];
   ymin = ypt[nxypt - 1];
   if (ymax == ymin)
      ymax = ymin + 1;
   dyy = (ymax - ymin) * .001;
   ymax += dyy;
   ymin -= dyy;
   mnbins(ymin, ymax, maxny, ymin, ymax, ny, bwidy);
   any = (double)ny;
   //*-*-        if first point is blank, it is an 'origin'
   if (chbest == ' ')
      goto L50;
   xbest = (xmax + xmin) * .5;
   ybest = (ymax + ymin) * .5;
L50:
   //*-*-        find Scale constants
   ax = 1 / bwidx;
   ay = 1 / bwidy;
   bx = -ax * xmin + 2;
   by = -ay * ymin - 2;
   //*-*-        convert points to grid positions
   for (i = 1; i <= nxypt; ++i) {
      xpt[i - 1] = ax * xpt[i - 1] + bx;
      ypt[i - 1] = any - ay * ypt[i - 1] - by;
   }
   nxbest = int((ax * xbest + bx));
   nybest = int((any - ay * ybest - by));
   //*-*-        print the points
   ny += 2;
   nx += 2;
   isp1 = 1;
   linodd = 1;
   overpr = false;
   for (i = 1; i <= ny; ++i) {
      for (ibk = 1; ibk <= nx; ++ibk) {
         cline[ibk - 1] = ' ';
      }
      cline[nx] = '\0';
      cline[nx + 1] = '\0';
      cline[0] = '.';
      // not needed - but to avoid a wrongly reported compiler warning (see ROOT-6496)
      if (nx > 0)
         cline[nx - 1] = '.';
      cline[nxbest - 1] = '.';
      if (i != 1 && i != nybest && i != ny)
         goto L320;
      for (j = 1; j <= nx; ++j) {
         cline[j - 1] = '.';
      }
   L320:
      yprt = ymax - double(i - 1) * bwidy;
      if (isp1 > nxypt)
         goto L350;
      //*-*-        find the points to be plotted on this line
      for (k = isp1; k <= nxypt; ++k) {
         ks = int(ypt[k - 1]);
         if (ks > i)
            goto L345;
         ix = int(xpt[k - 1]);
         if (cline[ix - 1] == '.')
            goto L340;
         if (cline[ix - 1] == ' ')
            goto L340;
         if (cline[ix - 1] == chpt[k - 1])
            continue;
         overpr = true;
         //*-*-        OVERPR is true if one or more positions contains more than
         //*-*-           one point
         cline[ix - 1] = '&';
         continue;
      L340:
         cline[ix - 1] = chpt[k - 1];
      }
      isp1 = nxypt + 1;
      goto L350;
   L345:
      isp1 = k;
   L350:
      if (linodd == 1 || i == ny)
         goto L380;
      linodd = 1;
      memcpy(ctemp, cline, 120);
      printf("                  %s", (const char *)ctemp);
      goto L400;
   L380:
      //   ctemp = cline;
      memcpy(ctemp, cline, 120);
      printf(" %14.7g ..%s", yprt, (const char *)ctemp);
      linodd = 0;
   L400:
      printf("\n");
   }
   //*-*-        print labels on x-axis every ten columns
   for (ibk = 1; ibk <= nx; ++ibk) {
      cline[ibk - 1] = ' ';
      if (ibk % 10 == 1)
         cline[ibk - 1] = '/';
   }
   printf("                  %s", cline);
   printf("\n");

   for (ibk = 1; ibk <= 12; ++ibk) {
      xvalus[ibk - 1] = xmin + double(ibk - 1) * 10 * bwidx;
   }
   printf("           ");
   iten = (nx + 9) / 10;
   for (ibk = 1; ibk <= iten; ++ibk) {
      printf(" %9.4g", xvalus[ibk - 1]);
   }
   printf("\n");

   if (overpr) {
      char chmess[] = "   Overprint character is &";
      printf("                         ONE COLUMN=%13.7g%s", bwidx, (const char *)chmess);
   } else {
      char chmess[] = " ";
      printf("                         ONE COLUMN=%13.7g%s", bwidx, (const char *)chmess);
   }
   printf("\n");

} /* mnplot_ */

} // namespace Minuit2

} // namespace ROOT
