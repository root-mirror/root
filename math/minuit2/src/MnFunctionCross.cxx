// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#include "Minuit2/MnFunctionCross.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/FCNBase.h"
#include "Minuit2/MnParabola.h"
#include "Minuit2/MnParabolaPoint.h"
#include "Minuit2/MnParabolaFactory.h"
#include "Minuit2/MnCross.h"
#include "Minuit2/MnMachinePrecision.h"
#include "Minuit2/MnPrint.h"

namespace ROOT {

namespace Minuit2 {

MnCross MnFunctionCross::operator()(const std::vector<unsigned int> &par, const std::vector<double> &pmid,
                                    const std::vector<double> &pdir, double tlr, unsigned int maxcalls) const
{
   // evaluate crossing point where function is equal to MIN + UP,
   // with direction pdir from values pmid
   // tlr indicate tolerance and maxcalls maximum number of calls

   //   double edmmax = 0.5*0.001*toler*fFCN.Up();

   unsigned int npar = par.size();
   unsigned int nfcn = 0;
   const MnMachinePrecision &prec = fState.Precision();
   // tolerance used when calling Migrad
   double mgr_tlr = 0.5 * tlr; // to be consistent with F77 version (for default values of tlr which is 0.1)
   // other olerance values are fixed at 0.01
   tlr = 0.01;
   // convergence when F is within tlf of aim and next prediction
   // of aopt is within tla of previous value of aopt
   double up = fFCN.Up();
   // for finding the point :
   double tlf = tlr * up;
   double tla = tlr;
   unsigned int maxitr = 15;
   unsigned int ipt = 0;
   double aminsv = fFval;
   double aim = aminsv + up;
   // std::cout<<"aim= "<<aim<<std::endl;
   double aopt = 0.;
   bool limset = false;
   std::vector<double> alsb(3, 0.), flsb(3, 0.);

   MnPrint print("MnFunctionCross");

   print.Debug([&](std::ostream &os) {
      for (unsigned int i = 0; i < par.size(); ++i)
         os << "Parameter " << par[i] << " value " << pmid[i] << " dir " << pdir[i] << " function min = " << aminsv
            << " contur value aim = (fmin + up) = " << aim;
   });

   // find the largest allowed aulim

   double aulim = 100.;
   for (unsigned int i = 0; i < par.size(); i++) {
      unsigned int kex = par[i];
      if (fState.Parameter(kex).HasLimits()) {
         double zmid = pmid[i];
         double zdir = pdir[i];
         //       double zlim = 0.;
         if (zdir > 0. && fState.Parameter(kex).HasUpperLimit()) {
            double zlim = fState.Parameter(kex).UpperLimit();
            if (std::fabs(zdir) < fState.Precision().Eps()) {
               // we have a limit
               if (std::fabs(zlim - zmid) < fState.Precision().Eps())
                  limset = true;
               continue;
            }
            aulim = std::min(aulim, (zlim - zmid) / zdir);
         } else if (zdir < 0. && fState.Parameter(kex).HasLowerLimit()) {
            double zlim = fState.Parameter(kex).LowerLimit();
            if (std::fabs(zdir) < fState.Precision().Eps()) {
               // we have a limit
               if (std::fabs(zlim - zmid) < fState.Precision().Eps())
                  limset = true;
               continue;
            }
            aulim = std::min(aulim, (zlim - zmid) / zdir);
         }
      }
   }

   print.Debug("Largest allowed aulim", aulim);

   // case of a single parameter and we are at limit
   if (limset && npar == 1) {
      print.Warn("Parameter is at limit", pmid[0], "delta", pdir[0]);
      return MnCross(fState, nfcn, MnCross::CrossParLimit());
   }

   if (aulim < aopt + tla)
      limset = true;

   MnMigrad migrad(fFCN, fState, MnStrategy(std::max(0, int(fStrategy.Strategy() - 1))));

   print.Info([&](std::ostream &os) {
      os << "Run Migrad with fixed parameters:";
      for (unsigned i = 0; i < npar; ++i)
         os << "\n  Pos " << par[i] << ": " << fState.Name(par[i]) << " = " << pmid[i];
   });

   for (unsigned int i = 0; i < npar; i++)
      migrad.SetValue(par[i], pmid[i]);

   // find minimum with respect all the other parameters (n- npar) (npar are the fixed ones)

   FunctionMinimum min0 = migrad(maxcalls, mgr_tlr);
   nfcn += min0.NFcn();

   print.Info("Result after Migrad", MnPrint::Oneline(min0), min0.UserState().Parameters());

   // case a new minimum is found
   if (min0.Fval() < fFval - tlf) {
      // case of new minimum is found
      print.Warn("New minimum found while scanning parameter", par.front(), "new value =", min0.Fval(),
                 "old value =", fFval);
      return MnCross(min0.UserState(), nfcn, MnCross::CrossNewMin());
   }
   if (min0.HasReachedCallLimit())
      return MnCross(min0.UserState(), nfcn, MnCross::CrossFcnLimit());
   if (!min0.IsValid())
      return MnCross(fState, nfcn);
   if (limset == true && min0.Fval() < aim)
      return MnCross(min0.UserState(), nfcn, MnCross::CrossParLimit());

   ipt++;
   alsb[0] = 0.;
   flsb[0] = min0.Fval();
   flsb[0] = std::max(flsb[0], aminsv + 0.1 * up);
   aopt = std::sqrt(up / (flsb[0] - aminsv)) - 1.;
   if (std::fabs(flsb[0] - aim) < tlf)
      return MnCross(aopt, min0.UserState(), nfcn);

   if (aopt > 1.)
      aopt = 1.;
   if (aopt < -0.5)
      aopt = -0.5;
   limset = false;
   if (aopt > aulim) {
      aopt = aulim;
      limset = true;
   }

   print.Debug("flsb[0]", flsb[0], "aopt", aopt);

   print.Info([&](std::ostream &os) {
      os << "Run Migrad again (2nd) with fixed parameters:";
      for (unsigned i = 0; i < npar; ++i)
         os << "\n  Pos " << par[i] << ": " << fState.Name(par[i]) << " = " << pmid[i] + (aopt)*pdir[i];
   });

   for (unsigned int i = 0; i < npar; i++)
      migrad.SetValue(par[i], pmid[i] + (aopt)*pdir[i]);

   FunctionMinimum min1 = migrad(maxcalls, mgr_tlr);
   nfcn += min1.NFcn();

   print.Info("Result after 2nd Migrad", MnPrint::Oneline(min1), min1.UserState().Parameters());

   if (min1.Fval() < fFval - tlf) // case of new minimum found
      return MnCross(min1.UserState(), nfcn, MnCross::CrossNewMin());
   if (min1.HasReachedCallLimit())
      return MnCross(min1.UserState(), nfcn, MnCross::CrossFcnLimit());
   if (!min1.IsValid())
      return MnCross(fState, nfcn);
   if (limset == true && min1.Fval() < aim)
      return MnCross(min1.UserState(), nfcn, MnCross::CrossParLimit());

   ipt++;
   alsb[1] = aopt;
   flsb[1] = min1.Fval();
   double dfda = (flsb[1] - flsb[0]) / (alsb[1] - alsb[0]);

   print.Debug("aopt", aopt, "min1Val", flsb[1], "dfda", dfda);

L300:
   if (dfda < 0.) {
      // looking for slope of the right sign
      print.Debug("dfda < 0 - iterate from", ipt, "to max of", maxitr);
      // iterate (max times is maxitr) incrementing aopt

      unsigned int maxlk = maxitr - ipt;
      for (unsigned int it = 0; it < maxlk; it++) {
         alsb[0] = alsb[1];
         flsb[0] = flsb[1];
         // LM: Add + 1, looking at Fortran code it starts from 1 ( see bug #8396)
         aopt = alsb[0] + 0.2 * (it + 1);
         limset = false;
         if (aopt > aulim) {
            aopt = aulim;
            limset = true;
         }

         print.Info([&](std::ostream &os) {
            os << "Run Migrad again (iteration " << it << " ) :";
            for (unsigned i = 0; i < npar; ++i)
               os << "\n  parameter " << par[i] << " (" << fState.Name(par[i]) << ") fixed to "
                  << pmid[i] + (aopt)*pdir[i];
         });

         for (unsigned int i = 0; i < npar; i++)
            migrad.SetValue(par[i], pmid[i] + (aopt)*pdir[i]);

         min1 = migrad(maxcalls, mgr_tlr);
         nfcn += min1.NFcn();

         print.Info("Result after Migrad", MnPrint::Oneline(min1), '\n', min1.UserState().Parameters());

         if (min1.Fval() < fFval - tlf) // case of new minimum found
            return MnCross(min1.UserState(), nfcn, MnCross::CrossNewMin());
         if (min1.HasReachedCallLimit())
            return MnCross(min1.UserState(), nfcn, MnCross::CrossFcnLimit());
         if (!min1.IsValid())
            return MnCross(fState, nfcn);
         if (limset == true && min1.Fval() < aim)
            return MnCross(min1.UserState(), nfcn, MnCross::CrossParLimit());
         ipt++;
         alsb[1] = aopt;
         flsb[1] = min1.Fval();
         dfda = (flsb[1] - flsb[0]) / (alsb[1] - alsb[0]);
         //       if(dfda > 0.) goto L460;

         print.Debug("aopt", aopt, "min1Val", flsb[1], "dfda", dfda);

         if (dfda > 0.)
            break;
      }
      if (ipt > maxitr)
         return MnCross(fState, nfcn);
   } // if(dfda < 0.)

L460:

   // dfda > 0: we have two points with the right slope

   aopt = alsb[1] + (aim - flsb[1]) / dfda;

   print.Debug("dfda > 0 : aopt", aopt);

   double fdist = std::min(std::fabs(aim - flsb[0]), std::fabs(aim - flsb[1]));
   double adist = std::min(std::fabs(aopt - alsb[0]), std::fabs(aopt - alsb[1]));
   tla = tlr;
   if (std::fabs(aopt) > 1.)
      tla = tlr * std::fabs(aopt);
   if (adist < tla && fdist < tlf)
      return MnCross(aopt, min1.UserState(), nfcn);
   if (ipt > maxitr)
      return MnCross(fState, nfcn);
   double bmin = std::min(alsb[0], alsb[1]) - 1.;
   if (aopt < bmin)
      aopt = bmin;
   double bmax = std::max(alsb[0], alsb[1]) + 1.;
   if (aopt > bmax)
      aopt = bmax;

   limset = false;
   if (aopt > aulim) {
      aopt = aulim;
      limset = true;
   }

   print.Info([&](std::ostream &os) {
      os << "Run Migrad again (3rd) with fixed parameters:";
      for (unsigned i = 0; i < npar; ++i)
         os << "\n  Pos " << par[i] << ": " << fState.Name(par[i]) << " = " << pmid[i] + (aopt)*pdir[i];
   });

   for (unsigned int i = 0; i < npar; i++)
      migrad.SetValue(par[i], pmid[i] + (aopt)*pdir[i]);

   FunctionMinimum min2 = migrad(maxcalls, mgr_tlr);
   nfcn += min2.NFcn();

   print.Info("Result after Migrad (3rd):", MnPrint::Oneline(min2), min2.UserState().Parameters());

   if (min2.Fval() < fFval - tlf) // case of new minimum found
      return MnCross(min2.UserState(), nfcn, MnCross::CrossNewMin());
   if (min2.HasReachedCallLimit())
      return MnCross(min2.UserState(), nfcn, MnCross::CrossFcnLimit());
   if (!min2.IsValid())
      return MnCross(fState, nfcn);
   if (limset == true && min2.Fval() < aim)
      return MnCross(min2.UserState(), nfcn, MnCross::CrossParLimit());

   ipt++;
   alsb[2] = aopt;
   flsb[2] = min2.Fval();

   // now we have three points, ask how many < AIM

   double ecarmn = std::fabs(flsb[2] - aim);
   double ecarmx = 0.;
   unsigned int ibest = 2;
   unsigned int iworst = 0;
   unsigned int noless = 0;

   for (unsigned int i = 0; i < 3; i++) {
      double ecart = std::fabs(flsb[i] - aim);
      if (ecart > ecarmx) {
         ecarmx = ecart;
         iworst = i;
      }
      if (ecart < ecarmn) {
         ecarmn = ecart;
         ibest = i;
      }
      if (flsb[i] < aim)
         noless++;
   }

   print.Debug("have three points : noless < aim; noless", noless, "ibest", ibest, "iworst", iworst);

   // std::cout<<"480"<<std::endl;

   // at least one on each side of AIM (contour), fit a parabola
   if (noless == 1 || noless == 2)
      goto L500;
   // if all three are above AIM, third point must be the closest to AIM, return it
   if (noless == 0 && ibest != 2)
      return MnCross(fState, nfcn);
   // if all three below and third is not best then the slope has again gone negative,
   // re-iterate and look for positive slope
   if (noless == 3 && ibest != 2) {
      alsb[1] = alsb[2];
      flsb[1] = flsb[2];

      print.Debug("All three points below - look again fir positive slope");
      goto L300;
   }

   // in other case new straight line thru first two points

   flsb[iworst] = flsb[2];
   alsb[iworst] = alsb[2];
   dfda = (flsb[1] - flsb[0]) / (alsb[1] - alsb[0]);

   print.Debug("New straight line using point 1-2; dfda", dfda);

   goto L460;

L500:

   do {
      // do parabola fit
      MnParabola parbol = MnParabolaFactory()(MnParabolaPoint(alsb[0], flsb[0]), MnParabolaPoint(alsb[1], flsb[1]),
                                              MnParabolaPoint(alsb[2], flsb[2]));
      //   aopt = parbol.X_pos(aim);
      // std::cout<<"alsb1,2,3= "<<alsb[0]<<", "<<alsb[1]<<", "<<alsb[2]<<std::endl;
      // std::cout<<"flsb1,2,3= "<<flsb[0]<<", "<<flsb[1]<<", "<<flsb[2]<<std::endl;

      print.Debug("Parabola fit: iteration", ipt);

      double coeff1 = parbol.C();
      double coeff2 = parbol.B();
      double coeff3 = parbol.A();
      double determ = coeff2 * coeff2 - 4. * coeff3 * (coeff1 - aim);

      print.Debug("Parabola fit: a =", coeff3, "b =", coeff2, "c =", coeff1, "determ =", determ);

      // curvature is negative
      if (determ < prec.Eps())
         return MnCross(fState, nfcn);
      double rt = std::sqrt(determ);
      double x1 = (-coeff2 + rt) / (2. * coeff3);
      double x2 = (-coeff2 - rt) / (2. * coeff3);
      double s1 = coeff2 + 2. * x1 * coeff3;
      double s2 = coeff2 + 2. * x2 * coeff3;

      print.Debug("Parabola fit: x1", x1, "x2", x2, "s1", s1, "s2", s2);

      if (s1 * s2 > 0.)
         print.Warn("Problem 1");

      // find with root is the right one
      aopt = x1;
      double slope = s1;
      if (s2 > 0.) {
         aopt = x2;
         slope = s2;
      }

      print.Debug("Parabola fit: aopt", aopt, "slope", slope);

      // ask if converged
      tla = tlr;
      if (std::fabs(aopt) > 1.)
         tla = tlr * std::fabs(aopt);

      print.Debug("Delta(aopt)", std::fabs(aopt - alsb[ibest]), "tla", tla, "Delta(F)", std::fabs(flsb[ibest] - aim),
                  "tlf", tlf);

      if (std::fabs(aopt - alsb[ibest]) < tla && std::fabs(flsb[ibest] - aim) < tlf)
         return MnCross(aopt, min2.UserState(), nfcn);

      //     if(ipt > maxitr) return MnCross();

      // see if proposed point is in acceptable zone between L and R
      // first find ileft, iright, iout and ibest

      unsigned int ileft = 3;
      unsigned int iright = 3;
      unsigned int iout = 3;
      ibest = 0;
      ecarmx = 0.;
      ecarmn = std::fabs(aim - flsb[0]);
      for (unsigned int i = 0; i < 3; i++) {
         double ecart = std::fabs(flsb[i] - aim);
         if (ecart < ecarmn) {
            ecarmn = ecart;
            ibest = i;
         }
         if (ecart > ecarmx)
            ecarmx = ecart;
         if (flsb[i] > aim) {
            if (iright == 3)
               iright = i;
            else if (flsb[i] > flsb[iright])
               iout = i;
            else {
               iout = iright;
               iright = i;
            }
         } else if (ileft == 3)
            ileft = i;
         else if (flsb[i] < flsb[ileft])
            iout = i;
         else {
            iout = ileft;
            ileft = i;
         }
      }

      print.Debug("ileft", ileft, "iright", iright, "iout", iout, "ibest", ibest);

      // avoid keeping a bad point nest time around

      if (ecarmx > 10. * std::fabs(flsb[iout] - aim))
         aopt = 0.5 * (aopt + 0.5 * (alsb[iright] + alsb[ileft]));

      // knowing ileft and iright, get acceptable window
      double smalla = 0.1 * tla;
      if (slope * smalla > tlf)
         smalla = tlf / slope;
      double aleft = alsb[ileft] + smalla;
      double aright = alsb[iright] - smalla;

      // move proposed point AOPT into window if necessary
      if (aopt < aleft)
         aopt = aleft;
      if (aopt > aright)
         aopt = aright;
      if (aleft > aright)
         aopt = 0.5 * (aleft + aright);

      // see if proposed point outside limits (should be impossible)
      limset = false;
      if (aopt > aulim) {
         aopt = aulim;
         limset = true;
      }

      // evaluate at new point aopt
      print.Info([&](std::ostream &os) {
         os << "Run Migrad again at new point (#iter = " << ipt << " ):";
         for (unsigned i = 0; i < npar; ++i)
            os << "\n\t - parameter " << i << " fixed to " << pmid[i] + (aopt)*pdir[i];
      });

      for (unsigned int i = 0; i < npar; i++)
         migrad.SetValue(par[i], pmid[i] + (aopt)*pdir[i]);

      min2 = migrad(maxcalls, mgr_tlr);
      nfcn += min2.NFcn();

      print.Info("Result after new Migrad:", MnPrint::Oneline(min2), min2.UserState().Parameters());

      if (min2.Fval() < fFval - tlf) // case of new minimum found
         return MnCross(min2.UserState(), nfcn, MnCross::CrossNewMin());
      if (min2.HasReachedCallLimit())
         return MnCross(min2.UserState(), nfcn, MnCross::CrossFcnLimit());
      if (!min2.IsValid())
         return MnCross(fState, nfcn);
      if (limset == true && min2.Fval() < aim)
         return MnCross(min2.UserState(), nfcn, MnCross::CrossParLimit());

      ipt++;
      // replace odd point with new one (which is the best of three)
      alsb[iout] = aopt;
      flsb[iout] = min2.Fval();
      ibest = iout;
   } while (ipt < maxitr);

   // goto L500;

   return MnCross(fState, nfcn);
}

} // namespace Minuit2

} // namespace ROOT
