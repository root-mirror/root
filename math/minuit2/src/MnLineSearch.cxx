// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#include "Minuit2/MnLineSearch.h"
#include "Minuit2/MnFcn.h"
#include "Minuit2/MinimumParameters.h"
#include "Minuit2/MnMachinePrecision.h"
#include "Minuit2/MnParabola.h"
#include "Minuit2/MnParabolaPoint.h"
#include "Minuit2/MnParabolaFactory.h"
#include "Minuit2/LaSum.h"
#include "Minuit2/MnPrint.h"

#include <algorithm>
#include <cmath>

#ifdef USE_OTHER_LS

#include "Math/SMatrix.h"
#include "Math/SVector.h"

#include "Math/IFunction.h"
#include "Math/Minimizer1D.h"

#endif

namespace ROOT {

namespace Minuit2 {

/**  Perform a line search from position defined by the vector st
       along the direction step, where the length of vector step
       gives the expected position of Minimum.
       fcn is Value of function at the starting position ,
       gdel (if non-zero) is df/dx along step at st.
       Return a parabola point containing Minimum x position and y (function Value)
    - add a falg to control the debug
*/

MnParabolaPoint MnLineSearch::operator()(const MnFcn &fcn, const MinimumParameters &st, const MnAlgebraicVector &step,
                                         double gdel, const MnMachinePrecision &prec) const
{

   //*-*-*-*-*-*-*-*-*-*Perform a line search from position st along step   *-*-*-*-*-*-*-*
   //*-*                =========================================
   //*-* SLAMBG and ALPHA control the maximum individual steps allowed.
   //*-* The first step is always =1. The max length of second step is SLAMBG.
   //*-* The max size of subsequent steps is the maximum previous successful
   //*-*   step multiplied by ALPHA + the size of most recent successful step,
   //*-*   but cannot be smaller than SLAMBG.
   //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

   MnPrint print("MnLineSearch");

   print.Debug("gdel", gdel, "step", step);

   double overal = 1000.;
   double undral = -100.;
   double toler = 0.05;
   double slamin = 0.;
   double slambg = 5.;
   double alpha = 2.;
   int maxiter = 12;
   // start as in Fortran from 1 and count all the time we evaluate the function
   int niter = 1;

   for (unsigned int i = 0; i < step.size(); i++) {
      if (step(i) == 0)
         continue;
      double ratio = std::fabs(st.Vec()(i) / step(i));
      if (slamin == 0)
         slamin = ratio;
      if (ratio < slamin)
         slamin = ratio;
   }
   if (std::fabs(slamin) < prec.Eps())
      slamin = prec.Eps();
   slamin *= prec.Eps2();

   double f0 = st.Fval();
   double f1 = fcn(st.Vec() + step);
   niter++;
   double fvmin = st.Fval();
   double xvmin = 0.;

   if (f1 < f0) {
      fvmin = f1;
      xvmin = 1.;
   }
   double toler8 = toler;
   double slamax = slambg;
   double flast = f1;
   double slam = 1.;

   bool iterate = false;
   MnParabolaPoint p0(0., f0);
   MnParabolaPoint p1(slam, flast);
   double f2 = 0.;
   // quadratic interpolation using the two points p0,p1 and the slope at p0
   do {
      // cut toler8 as function goes up
      iterate = false;
      // MnParabola pb = MnParabolaFactory()(p0, gdel, p1);

      print.Debug("flast", flast, "f0", f0, "flast-f0", flast - f0, "slam", slam);
      //     double df = flast-f0;
      //     if(std::fabs(df) < prec.Eps2()) {
      //       if(flast-f0 < 0.) df = -prec.Eps2();
      //       else df = prec.Eps2();
      //     }
      //     std::cout<<"df= "<<df<<std::endl;
      //     double denom = 2.*(df-gdel*slam)/(slam*slam);
      double denom = 2. * (flast - f0 - gdel * slam) / (slam * slam);

      print.Debug("denom", denom);
      if (denom != 0) {
         slam = -gdel / denom;
      } else {
         denom = -0.1 * gdel;
         slam = 1.;
      }
      print.Debug("new slam", slam);

#ifdef TRY_OPPOSIT_DIR
      // if is less than zero indicates maximum position. Use then slamax or x = x0 - 2slam + 1
      bool slamIsNeg = false;
      double slamNeg = 0;
#endif

      if (slam < 0.) {
         print.Debug("slam is negative - set to", slamax);
#ifdef TRY_OPPOSITE_DIR
         slamNeg = 2 * slam - 1;
         slamIsNeg = true;
         print.Debug("slam is negative - compare values between", slamNeg, "and", slamax);
#endif
         slam = slamax;
      }
      //     std::cout<<"slam= "<<slam<<std::endl;
      if (slam > slamax) {
         slam = slamax;
         print.Debug("slam larger than mac value - set to", slamax);
      }

      if (slam < toler8) {
         print.Debug("slam too small - set to", toler8);
         slam = toler8;
      }
      //     std::cout<<"slam= "<<slam<<std::endl;
      if (slam < slamin) {
         print.Debug("slam smaller than", slamin, "return");
         //       std::cout<<"f1, f2= "<<p0.Y()<<", "<<p1.Y()<<std::endl;
         //       std::cout<<"x1, x2= "<<p0.X()<<", "<<p1.X()<<std::endl;
         //       std::cout<<"x, f= "<<xvmin<<", "<<fvmin<<std::endl;
         return MnParabolaPoint(xvmin, fvmin);
      }
      if (std::fabs(slam - 1.) < toler8 && p1.Y() < p0.Y()) {
         //       std::cout<<"f1, f2= "<<p0.Y()<<", "<<p1.Y()<<std::endl;
         //       std::cout<<"x1, x2= "<<p0.X()<<", "<<p1.X()<<std::endl;
         //       std::cout<<"x, f= "<<xvmin<<", "<<fvmin<<std::endl;
         return MnParabolaPoint(xvmin, fvmin);
      }
      if (std::fabs(slam - 1.) < toler8)
         slam = 1. + toler8;

      //     if(std::fabs(gdel) < prec.Eps2() && std::fabs(denom) < prec.Eps2())
      //       slam = 1000.;
      //     MnAlgebraicVector tmp = step;
      //     tmp *= slam;
      //     f2 = fcn(st.Vec()+tmp);
      f2 = fcn(st.Vec() + slam * step);

      niter++; // do as in Minuit (count all func evalu)

#ifdef TRY_OPPOSITE_DIR
      if (slamIsNeg) {
         // try alternative in the opposite direction
         double f2alt = fcn(st.Vec() + slamNeg * step);
         if (f2alt < f2) {
            slam = slamNeg;
            f2 = f2alt;
            undral += slam;
         }
      }
#endif
      if (f2 < fvmin) {
         fvmin = f2;
         xvmin = slam;
      }
      // LM : correct a bug using precision
      if (std::fabs(p0.Y() - fvmin) < std::fabs(fvmin) * prec.Eps()) {
         //   if(p0.Y()-prec.Eps() < fvmin && fvmin < p0.Y()+prec.Eps()) {
         iterate = true;
         flast = f2;
         toler8 = toler * slam;
         overal = slam - toler8;
         slamax = overal;
         p1 = MnParabolaPoint(slam, flast);
         // niter++;
      }
   } while (iterate && niter < maxiter);
   if (niter >= maxiter) {
      // exhausted max number of iterations
      return MnParabolaPoint(xvmin, fvmin);
   }

   print.Debug("after initial 2-point iter:", '\n', " x0, x1, x2:", p0.X(), p1.X(), slam, '\n', " f0, f1, f2:", p0.Y(),
               p1.Y(), f2);

   MnParabolaPoint p2(slam, f2);

   // do now the quadratic interpolation with 3 points
   do {
      slamax = std::max(slamax, alpha * std::fabs(xvmin));
      MnParabola pb = MnParabolaFactory()(p0, p1, p2);
      print.Debug("Iteration", niter, '\n', " x0, x1, x2:", p0.X(), p1.X(), p2.X(), '\n', " f0, f1, f2:", p0.Y(),
                  p1.Y(), p2.Y(), '\n', " slamax    :", slamax, '\n', " p2-p0,p1  :", p2.Y() - p0.Y(), p2.Y() - p1.Y(),
                  '\n', " a, b, c   :", pb.A(), pb.B(), pb.C());
      if (pb.A() < prec.Eps2()) {
         double slopem = 2. * pb.A() * xvmin + pb.B();
         if (slopem < 0.)
            slam = xvmin + slamax;
         else
            slam = xvmin - slamax;
         print.Debug("xvmin", xvmin, "slopem", slopem, "slam", slam);
      } else {
         slam = pb.Min();
         //      std::cout<<"pb.Min() slam= "<<slam<<std::endl;
         if (slam > xvmin + slamax)
            slam = xvmin + slamax;
         if (slam < xvmin - slamax)
            slam = xvmin - slamax;
      }
      if (slam > 0.) {
         if (slam > overal)
            slam = overal;
      } else {
         if (slam < undral)
            slam = undral;
      }

      print.Debug("slam", slam, "undral", undral, "overal", overal);

      double f3 = 0.;
      do {

         print.Debug("iterate on f3- slam", niter, "slam", slam, "xvmin", xvmin);

         iterate = false;
         double toler9 = std::max(toler8, std::fabs(toler8 * slam));
         // min. of parabola at one point
         if (std::fabs(p0.X() - slam) < toler9 || std::fabs(p1.X() - slam) < toler9 ||
             std::fabs(p2.X() - slam) < toler9) {
            //   std::cout<<"f1, f2, f3= "<<p0.Y()<<", "<<p1.Y()<<", "<<p2.Y()<<std::endl;
            //   std::cout<<"x1, x2, x3= "<<p0.X()<<", "<<p1.X()<<", "<<p2.X()<<std::endl;
            //   std::cout<<"x, f= "<<xvmin<<", "<<fvmin<<std::endl;
            return MnParabolaPoint(xvmin, fvmin);
         }

         // take the step
         //       MnAlgebraicVector tmp = step;
         //       tmp *= slam;
         f3 = fcn(st.Vec() + slam * step);
         print.Debug("f3", f3, "f3-p(2-0).Y()", f3 - p2.Y(), f3 - p1.Y(), f3 - p0.Y());
         // if latest point worse than all three previous, cut step
         if (f3 > p0.Y() && f3 > p1.Y() && f3 > p2.Y()) {
            print.Debug("f3 worse than all three previous");
            if (slam > xvmin)
               overal = std::min(overal, slam - toler8);
            if (slam < xvmin)
               undral = std::max(undral, slam + toler8);
            slam = 0.5 * (slam + xvmin);
            print.Debug("new slam", slam);
            iterate = true;
            niter++;
         }
      } while (iterate && niter < maxiter);
      if (niter >= maxiter) {
         // exhausted max number of iterations
         return MnParabolaPoint(xvmin, fvmin);
      }

      // find worst previous point out of three and replace
      MnParabolaPoint p3(slam, f3);
      if (p0.Y() > p1.Y() && p0.Y() > p2.Y())
         p0 = p3;
      else if (p1.Y() > p0.Y() && p1.Y() > p2.Y())
         p1 = p3;
      else
         p2 = p3;
      print.Debug("f3", f3, "fvmin", fvmin, "xvmin", xvmin);
      if (f3 < fvmin) {
         fvmin = f3;
         xvmin = slam;
      } else {
         if (slam > xvmin)
            overal = std::min(overal, slam - toler8);
         if (slam < xvmin)
            undral = std::max(undral, slam + toler8);
      }

      niter++;
   } while (niter < maxiter);

   print.Debug("f1, f2 =", p0.Y(), p1.Y(), '\n', "x1, x2 =", p0.X(), p1.X(), '\n', "x, f =", xvmin, fvmin);
   return MnParabolaPoint(xvmin, fvmin);
}

#ifdef USE_OTHER_LS

/**  Perform a line search using a cubic interpolation using x0, x1 , df/dx(x0) and d2/dx(x0) (second derivative)
     This is used at the beginning when the second derivative is known to be negative
*/

MnParabolaPoint MnLineSearch::CubicSearch(const MnFcn &fcn, const MinimumParameters &st, const MnAlgebraicVector &step,
                                          double gdel, double g2del, const MnMachinePrecision &prec) const
{
   MnPrint print("MnLineSearch::CubicSearch");

   print.Debug("gdel", gdel, "g2del", g2del, "step", step);

   // change ot large values
   double overal = 100.;
   double undral = -100.;
   double toler = 0.05;
   double slamin = 0.;
   double slambg = 5.;
   double alpha = 2.;

   for (unsigned int i = 0; i < step.size(); i++) {
      if (step(i) == 0)
         continue;
      double ratio = std::fabs(st.Vec()(i) / step(i));
      if (slamin == 0)
         slamin = ratio;
      if (ratio < slamin)
         slamin = ratio;
   }
   if (std::fabs(slamin) < prec.Eps())
      slamin = prec.Eps();
   slamin *= prec.Eps2();

   double f0 = st.Fval();
   double f1 = fcn(st.Vec() + step);
   double fvmin = st.Fval();
   double xvmin = 0.;
   print.Debug("f0", f0, "f1", f1);

   if (f1 < f0) {
      fvmin = f1;
      xvmin = 1.;
   }
   double toler8 = toler;
   double slamax = slambg;
   double flast = f1;
   double slam = 1.;

   //    MnParabolaPoint p0(0., f0);
   //    MnParabolaPoint p1(slam, flast);

   ROOT::Math::SMatrix<double, 3> cubicMatrix;
   ROOT::Math::SVector<double, 3> cubicCoeff; // cubic coefficients to be found
   ROOT::Math::SVector<double, 3> bVec;       // cubic coefficients to be found
   double x0 = 0;

   // cubic interpolation using the two points p0,p1 and the slope at p0 and the second derivative at p0

   // cut toler8 as function goes up
   double x1 = slam;
   cubicMatrix(0, 0) = (x0 * x0 * x0 - x1 * x1 * x1) / 3.;
   cubicMatrix(0, 1) = (x0 * x0 - x1 * x1) / 2.;
   cubicMatrix(0, 2) = (x0 - x1);
   cubicMatrix(1, 0) = x0 * x0;
   cubicMatrix(1, 1) = x0;
   cubicMatrix(1, 2) = 1;
   cubicMatrix(2, 0) = 2. * x0;
   cubicMatrix(2, 1) = 1;
   cubicMatrix(2, 2) = 0;

   bVec(0) = f0 - f1;
   bVec(1) = gdel;
   bVec(2) = g2del;

   // if (debug) std::cout << "Matrix:\n " << cubicMatrix << std::endl;
   print.Debug("Vec:\n  ", bVec);

   // find the minimum need to invert the matrix
   if (!cubicMatrix.Invert()) {
      print.Warn("Inversion failed - return");
      return MnParabolaPoint(xvmin, fvmin);
   }

   cubicCoeff = cubicMatrix * bVec;
   print.Debug("Cubic:\n   ", cubicCoeff);

   double ddd = cubicCoeff(1) * cubicCoeff(1) - 4 * cubicCoeff(0) * cubicCoeff(2); // b**2 - 4ac
   double slam1, slam2 = 0;
   // slam1 should be minimum and slam2 the maximum
   if (cubicCoeff(0) > 0) {
      slam1 = (-cubicCoeff(1) - std::sqrt(ddd)) / (2. * cubicCoeff(0));
      slam2 = (-cubicCoeff(1) + std::sqrt(ddd)) / (2. * cubicCoeff(0));
   } else if (cubicCoeff(0) < 0) {
      slam1 = (-cubicCoeff(1) + std::sqrt(ddd)) / (2. * cubicCoeff(0));
      slam2 = (-cubicCoeff(1) - std::sqrt(ddd)) / (2. * cubicCoeff(0));
   } else { // case == 0 (-b/c)
      slam1 = -gdel / g2del;
      slam2 = slam1;
   }

   print.Debug("slam1", slam1, "slam2", slam2);

   // this should be the minimum otherwise inversion failed and I should do something else

   if (slam2 < undral)
      slam2 = undral;
   if (slam2 > overal)
      slam2 = overal;

   // I am stack somewhere - take a large step
   if (std::fabs(slam2) < toler)
      slam2 = (slam2 >= 0) ? slamax : -slamax;

   double f2 = fcn(st.Vec() + slam2 * step);

   print.Debug("try with slam 2", slam2, "f2", f2);

   double fp;
   // use this as new minimum
   // bool noImpr = false;
   if (f2 < fvmin) {
      slam = slam2;
      xvmin = slam;
      fvmin = f2;
      fp = fvmin;
   } else {
      // try with slam2 if it is better

      if (slam1 < undral)
         slam1 = undral;
      if (slam1 > overal)
         slam1 = overal;

      if (std::fabs(slam1) < toler)
         slam1 = (slam1 >= 0) ? -slamax : slamax;

      double f3 = fcn(st.Vec() + slam1 * step);

      print.Debug("try with slam 1", slam1, "f3", f3);

      if (f3 < fvmin) {
         slam = slam1;
         fp = fvmin;
         xvmin = slam;
         fvmin = f3;
      } else {
         // case both f2 and f3 did not produce a better result
         if (f2 < f3) {
            slam = slam1;
            fp = f2;
         } else {
            slam = slam2;
            fp = f3;
         }
      }
   }

   bool iterate2 = false;
   int niter = 0;

   int maxiter = 10;

   do {
      iterate2 = false;

      print.Debug("iter", niter, "test approx deriv ad second deriv at", slam, "fp", fp);

      // estimate grad and second derivative at new point taking a step of 10-3
      double h = 0.001 * slam;
      double fh = fcn(st.Vec() + (slam + h) * step);
      double fl = fcn(st.Vec() + (slam - h) * step);
      double df = (fh - fl) / (2. * h);
      double df2 = (fh + fl - 2. * fp) / (h * h);

      print.Debug("deriv", df, df2);

      // if I am in a point of still negative derivative
      if (std::fabs(df) < prec.Eps() && std::fabs(df2) < prec.Eps()) {
         // try in opposite direction with an opposite value
         slam = (slam >= 0) ? -slamax : slamax;
         slamax *= 10;
         fp = fcn(st.Vec() + slam * step);
      } else if (std::fabs(df2) <= 0) { // gradient is significative different than zero then redo a cubic interpolation
                                        // from new point

         return MnParabolaPoint(slam, fp); // should redo a cubic interpol.  ??
                                           //          niter ++;
                                           //          if (niter > maxiter) break;

         //          MinimumParameters pa = MinimumParameters(st.Vec() + stepNew, fp);
         //          gdel = stepNew(i)
         //          MnParabolaPoint pp = CubicSearch(fcn, st, stepNew, df, df2

      }

      else
         return MnParabolaPoint(slam, fp);

      niter++;
   } while (niter < maxiter);

   return MnParabolaPoint(xvmin, fvmin);
}

// class describing Fcn function in one dimension
class ProjectedFcn : public ROOT::Math::IGenFunction {

public:
   ProjectedFcn(const MnFcn &fcn, const MinimumParameters &pa, const MnAlgebraicVector &step)
      : fFcn(fcn), fPar(pa), fStep(step)
   {
   }

   ROOT::Math::IGenFunction *Clone() const { return new ProjectedFcn(*this); }

private:
   double DoEval(double x) const { return fFcn(fPar.Vec() + x * fStep); }

   const MnFcn &fFcn;
   const MinimumParameters &fPar;
   const MnAlgebraicVector &fStep;
};

MnParabolaPoint MnLineSearch::BrentSearch(const MnFcn &fcn, const MinimumParameters &st, const MnAlgebraicVector &step,
                                          double gdel, double g2del, const MnMachinePrecision &prec) const
{
   MnPrint print("MnLineSearch::BrentSearch");

   print.Debug("gdel", gdel, "g2del", g2del);

   print.Debug([&](std::ostream &os) {
      for (unsigned int i = 0; i < step.size(); ++i) {
         if (step(i) != 0) {
            os << "step(i) " << step(i) << '\n';
            std::cout << "par(i) " << st.Vec()(i) << '\n';
            break;
         }
      }
   });

   ProjectedFcn func(fcn, st, step);

   // do first a cubic interpolation

   double f0 = st.Fval();
   double f1 = fcn(st.Vec() + step);
   f0 = func(0.0);
   f1 = func(1.);
   double fvmin = st.Fval();
   double xvmin = 0.;
   print.Debug("f0", f0, "f1", f1);

   if (f1 < f0) {
      fvmin = f1;
      xvmin = 1.;
   }
   //    double toler8 = toler;
   //    double slamax = slambg;
   //    double flast = f1;
   double slam = 1.;

   double undral = -1000;
   double overal = 1000;

   double x0 = 0;

//    MnParabolaPoint p0(0., f0);
//    MnParabolaPoint p1(slam, flast);
#ifdef USE_CUBIC

   ROOT::Math::SMatrix<double, 3> cubicMatrix;
   ROOT::Math::SVector<double, 3> cubicCoeff; // cubic coefficients to be found
   ROOT::Math::SVector<double, 3> bVec;       // cubic coefficients to be found

   // cubic interpolation using the two points p0,p1 and the slope at p0 and the second derivative at p0

   // cut toler8 as function goes up
   double x1 = slam;
   cubicMatrix(0, 0) = (x0 * x0 * x0 - x1 * x1 * x1) / 3.;
   cubicMatrix(0, 1) = (x0 * x0 - x1 * x1) / 2.;
   cubicMatrix(0, 2) = (x0 - x1);
   cubicMatrix(1, 0) = x0 * x0;
   cubicMatrix(1, 1) = x0;
   cubicMatrix(1, 2) = 1;
   cubicMatrix(2, 0) = 2. * x0;
   cubicMatrix(2, 1) = 1;
   cubicMatrix(2, 2) = 0;

   bVec(0) = f0 - f1;
   bVec(1) = gdel;
   bVec(2) = g2del;

   // if (debug) std::cout << "Matrix:\n " << cubicMatrix << std::endl;
   print.Debug("Vec:\n  ", bVec);

   // find the minimum need to invert the matrix
   if (!cubicMatrix.Invert()) {
      print.Warn("Inversion failed - return");
      return MnParabolaPoint(xvmin, fvmin);
   }

   cubicCoeff = cubicMatrix * bVec;
   print.Debug("Cubic:\n  ", cubicCoeff);

   double ddd = cubicCoeff(1) * cubicCoeff(1) - 4 * cubicCoeff(0) * cubicCoeff(2); // b**2 - 4ac
   double slam1, slam2 = 0;
   // slam1 should be minimum and slam2 the maximum
   if (cubicCoeff(0) > 0) {
      slam1 = (-cubicCoeff(1) - std::sqrt(ddd)) / (2. * cubicCoeff(0));
      slam2 = (-cubicCoeff(1) + std::sqrt(ddd)) / (2. * cubicCoeff(0));
   } else if (cubicCoeff(0) < 0) {
      slam1 = (-cubicCoeff(1) + std::sqrt(ddd)) / (2. * cubicCoeff(0));
      slam2 = (-cubicCoeff(1) - std::sqrt(ddd)) / (2. * cubicCoeff(0));
   } else { // case == 0 (-b/c)
      slam1 = -gdel / g2del;
      slam2 = slam1;
   }

   if (slam1 < undral)
      slam1 = undral;
   if (slam1 > overal)
      slam1 = overal;

   if (slam2 < undral)
      slam2 = undral;
   if (slam2 > overal)
      slam2 = overal;

   double fs1 = func(slam1);
   double fs2 = func(slam2);
   print.Debug("slam1", slam1, "slam2", slam2, "f(slam1)", fs1, "f(slam2)", fs2);

   if (fs1 < fs2) {
      x0 = slam1;
      f0 = fs1;
   } else {
      x0 = slam2;
      f0 = fs2;
   }

#else
   x0 = xvmin;
   f0 = fvmin;
#endif

   double astart = 100;

   // do a Brent search in a large interval
   double a = x0 - astart;
   double b = x0 + astart;
   // double x0 = 1;
   int maxiter = 20;
   double absTol = 0.5;
   double relTol = 0.05;

   ROOT::Math::Minim1D::Type minType = ROOT::Math::Minim1D::BRENT;

   bool iterate = false;
   int iter = 0;

   print.Debug("f(0)", func(0.), "f1", func(1.0), "f(3)", func(3.0), "f(5)", func(5.0));

   double fa = func(a);
   double fb = func(b);
   // double f0 = func(x0);
   double toler = 0.01;
   double delta0 = 5;
   double delta = delta0;
   bool enlarge = true;
   bool scanmin = false;
   double x2, f2 = 0;
   double dir = 1;
   int nreset = 0;

   do {

      print.Debug("iter", iter, "a", a, "b", b, "x0", x0, "fa", fa, "fb", fb, "f0", f0);
      if (fa <= f0 || fb <= f0) {
         scanmin = false;
         if (fa < fb) {
            if (fa < fvmin) {
               fvmin = fa;
               xvmin = a;
            }
            //             double f2 = fa;
            //             double x2 = a;
            if (enlarge) {
               x2 = a - 1000; // move lower
               f2 = func(x2);
            }
            if (std::fabs((fa - f2) / (a - x2)) > toler) { //  significant change in f continue to enlarge interval
               x0 = a;
               f0 = fa;
               a = x2;
               fa = f2;
               enlarge = true;
            } else {
               // move direction of [a
               // values change a little, start from central point try with x0 = x0 - delta
               if (nreset == 0)
                  dir = -1;
               enlarge = false;
               x0 = x0 + dir * delta;
               f0 = func(x0);
               // if reached limit try opposite direction , direction b]
               if (std::fabs((f0 - fa) / (x0 - a)) < toler) {
                  delta = 10 * delta0 / (10. * (nreset + 1)); // decrease the delta if still not change observed
                  a = x0;
                  fa = f0;
                  x0 = delta;
                  f0 = func(x0);
                  dir *= -1;
                  nreset++;
                  print.Info("A: Done a reset - scan in opposite direction!");
               }
               delta *= 2; // double delta at next iteration
               if (x0 < b && x0 > a)
                  scanmin = true;
               else {
                  x0 = 0;
                  f0 = st.Fval();
               }
            }
         } else { // fb < fa
            if (fb < fvmin) {
               fvmin = fb;
               xvmin = b;
            }
            if (enlarge) {
               x2 = b + 1000; // move upper
               f2 = func(x2);
            }
            if (std::fabs((fb - f2) / (x2 - b)) > toler) { //  significant change in f continue to enlarge interval
               f0 = fb;
               x0 = b;
               b = x2; //
               fb = f2;
               enlarge = true;
            } else {
               // here move direction b
               // values change a little, try with x0 = fa + delta
               if (nreset == 0)
                  dir = 1;
               enlarge = false;
               x0 = x0 + dir * delta;
               f0 = func(x0);
               // if reached limit try other side
               if (std::fabs((f0 - fb) / (x0 - b)) < toler) {
                  delta = 10 * delta0 / (10. * (nreset + 1)); // decrease the delta if still not change observed
                  b = x0;
                  fb = f0;
                  x0 = -delta;
                  f0 = func(x0);
                  dir *= -1;
                  nreset++;
                  print.Info("B: Done a reset - scan in opposite direction!");
               }
               delta *= 2; // double delta at next iteration
               if (x0 > a && x0 < b)
                  scanmin = true;
               else {
                  x0 = 0;
                  f0 = st.Fval();
               }
            }
         }

         if (f0 < fvmin) {
            x0 = xvmin;
            fvmin = f0;
         }

         print.Debug("new x0", x0, "f0", f0);

         // use golden section
         iterate = scanmin || enlarge;

      } else { // x0 is a min of [a,b]
         iterate = false;
      }

      iter++;
   } while (iterate && iter < maxiter);

   if (f0 > fa || f0 > fb)
      // skip minim 1d try Minuit LS
      // return (*this) (fcn, st, step, gdel, prec, debug);
      return MnParabolaPoint(xvmin, fvmin);

   print.Info("1D minimization using", minType);

   ROOT::Math::Minimizer1D min(minType);

   min.SetFunction(func, x0, a, b);
   int ret = min.Minimize(maxiter, absTol, relTol);

   MnPrint::info("result of GSL 1D minimization:", ret, "niter", min.Iterations(), "xmin", min.XMinimum(), "fmin",
                 min.FValMinimum());

   return MnParabolaPoint(min.XMinimum(), min.FValMinimum());
}

#endif

} // namespace Minuit2

} // namespace ROOT
