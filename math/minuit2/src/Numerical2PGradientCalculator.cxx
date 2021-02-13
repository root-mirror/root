// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#include "Minuit2/Numerical2PGradientCalculator.h"
#include "Minuit2/InitialGradientCalculator.h"
#include "Minuit2/MnFcn.h"
#include "Minuit2/MnUserTransformation.h"
#include "Minuit2/MnMachinePrecision.h"
#include "Minuit2/MinimumParameters.h"
#include "Minuit2/FunctionGradient.h"
#include "Minuit2/MnStrategy.h"
#include "Minuit2/MnPrint.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#include <cmath>
#include <cassert>
#include <iomanip>

#include "Minuit2/MPIProcess.h"

namespace ROOT {

namespace Minuit2 {

FunctionGradient Numerical2PGradientCalculator::operator()(const MinimumParameters &par) const
{
   // calculate gradient using Initial gradient calculator and from MinimumParameters object

   InitialGradientCalculator gc(fFcn, fTransformation, fStrategy);
   FunctionGradient gra = gc(par);

   return (*this)(par, gra);
}

// comment it, because it was added
FunctionGradient Numerical2PGradientCalculator::operator()(const std::vector<double> &params) const
{
   // calculate gradient from an std;:vector of paramteters

   int npar = params.size();

   MnAlgebraicVector par(npar);
   for (int i = 0; i < npar; ++i) {
      par(i) = params[i];
   }

   double fval = Fcn()(par);

   MinimumParameters minpars = MinimumParameters(par, fval);

   return (*this)(minpars);
}

FunctionGradient Numerical2PGradientCalculator::
operator()(const MinimumParameters &par, const FunctionGradient &Gradient) const
{
   // calculate numerical gradient from MinimumParameters object
   // the algorithm takes correctly care when the gradient is approximatly zero

   //    std::cout<<"########### Numerical2PDerivative"<<std::endl;
   //    std::cout<<"initial grd: "<<Gradient.Grad()<<std::endl;
   //    std::cout<<"position: "<<par.Vec()<<std::endl;
   MnPrint print("Numerical2PGradientCalculator");

   assert(par.IsValid());

   double fcnmin = par.Fval();
   //   std::cout<<"fval: "<<fcnmin<<std::endl;

   double eps2 = Precision().Eps2();
   double eps = Precision().Eps();

   print.Debug("Assumed precision eps", eps, "eps2", eps2);

   double dfmin = 8. * eps2 * (std::fabs(fcnmin) + Fcn().Up());
   double vrysml = 8. * eps * eps;
   //   double vrysml = std::max(1.e-4, eps2);
   //    std::cout<<"dfmin= "<<dfmin<<std::endl;
   //    std::cout<<"vrysml= "<<vrysml<<std::endl;
   //    std::cout << " ncycle " << Ncycle() << std::endl;

   unsigned int n = (par.Vec()).size();
   unsigned int ncycle = Ncycle();
   //   MnAlgebraicVector vgrd(n), vgrd2(n), vgstp(n);
   MnAlgebraicVector grd = Gradient.Grad();
   MnAlgebraicVector g2 = Gradient.G2();
   MnAlgebraicVector gstep = Gradient.Gstep();

   print.Debug("Calculating gradient around value", fcnmin, "at point", par.Vec());

#ifndef _OPENMP

   MPIProcess mpiproc(n, 0);

   // for serial execution this can be outside the loop
   MnAlgebraicVector x = par.Vec();

   unsigned int startElementIndex = mpiproc.StartElementIndex();
   unsigned int endElementIndex = mpiproc.EndElementIndex();

   for (unsigned int i = startElementIndex; i < endElementIndex; i++) {

#else

   // parallelize this loop using OpenMP
//#define N_PARALLEL_PAR 5
#pragma omp parallel
#pragma omp for
   //#pragma omp for schedule (static, N_PARALLEL_PAR)

   for (int i = 0; i < int(n); i++) {

#endif

#ifdef _OPENMP
      // create in loop since each thread will use its own copy
      MnAlgebraicVector x = par.Vec();
#endif

      double xtf = x(i);
      double epspri = eps2 + std::fabs(grd(i) * eps2);
      double stepb4 = 0.;
      for (unsigned int j = 0; j < ncycle; j++) {
         double optstp = std::sqrt(dfmin / (std::fabs(g2(i)) + epspri));
         double step = std::max(optstp, std::fabs(0.1 * gstep(i)));
         //       std::cout<<"step: "<<step;
         if (Trafo().Parameter(Trafo().ExtOfInt(i)).HasLimits()) {
            if (step > 0.5)
               step = 0.5;
         }
         double stpmax = 10. * std::fabs(gstep(i));
         if (step > stpmax)
            step = stpmax;
         //       std::cout<<" "<<step;
         double stpmin = std::max(vrysml, 8. * std::fabs(eps2 * x(i)));
         if (step < stpmin)
            step = stpmin;
         //       std::cout<<" "<<step<<std::endl;
         //       std::cout<<"step: "<<step<<std::endl;
         if (std::fabs((step - stepb4) / step) < StepTolerance()) {
            //    std::cout<<"(step-stepb4)/step"<<std::endl;
            //    std::cout<<"j= "<<j<<std::endl;
            //    std::cout<<"step= "<<step<<std::endl;
            break;
         }
         gstep(i) = step;
         stepb4 = step;
         //       MnAlgebraicVector pstep(n);
         //       pstep(i) = step;
         //       double fs1 = Fcn()(pstate + pstep);
         //       double fs2 = Fcn()(pstate - pstep);

         x(i) = xtf + step;
         double fs1 = Fcn()(x);
         x(i) = xtf - step;
         double fs2 = Fcn()(x);
         x(i) = xtf;

         double grdb4 = grd(i);
         grd(i) = 0.5 * (fs1 - fs2) / step;
         g2(i) = (fs1 + fs2 - 2. * fcnmin) / step / step;

#ifdef _OPENMP
#pragma omp critical
#endif
         {
#ifdef _OPENMP
            // must create thread-local MnPrint instances when printing inside threads
            MnPrint print("Numerical2PGradientCalculator[OpenMP]");
#endif
            if (i == 0 && j == 0) {
               print.Debug([&](std::ostream &os) {
                  os << std::setw(10) << "parameter" << std::setw(6) << "cycle" << std::setw(15) << "x" << std::setw(15)
                     << "step" << std::setw(15) << "f1" << std::setw(15) << "f2" << std::setw(15) << "grd"
                     << std::setw(15) << "g2" << std::endl;
               });
            }
            print.Debug([&](std::ostream &os) {
               const int pr = os.precision(13);
               const int iext = Trafo().ExtOfInt(i);
               os << std::setw(10) << Trafo().Name(iext) << std::setw(5) << j << "  " << x(i) << " " << step << " "
                  << fs1 << " " << fs2 << " " << grd(i) << " " << g2(i) << std::endl;
               os.precision(pr);
            });
         }

         if (std::fabs(grdb4 - grd(i)) / (std::fabs(grd(i)) + dfmin / step) < GradTolerance()) {
            //    std::cout<<"j= "<<j<<std::endl;
            //    std::cout<<"step= "<<step<<std::endl;
            //    std::cout<<"fs1, fs2: "<<fs1<<" "<<fs2<<std::endl;
            //    std::cout<<"fs1-fs2: "<<fs1-fs2<<std::endl;
            break;
         }
      }

      //     vgrd(i) = grd;
      //     vgrd2(i) = g2;
      //     vgstp(i) = gstep;
   }

#ifndef _OPENMP
   mpiproc.SyncVector(grd);
   mpiproc.SyncVector(g2);
   mpiproc.SyncVector(gstep);
#endif

   // print after parallel processing to avoid synchronization issues
   print.Debug([&](std::ostream &os) {
      const int pr = os.precision(13);
      os << std::endl;
      os << std::setw(14) << "Parameter" << std::setw(14) << "Gradient" << std::setw(14) << "g2 " << std::setw(14)
         << "step" << std::endl;
      for (int i = 0; i < int(n); i++) {
         const int iext = Trafo().ExtOfInt(i);
         os << std::setw(14) << Trafo().Name(iext) << " " << grd(i) << " " << g2(i) << " " << gstep(i) << std::endl;
      }
      os.precision(pr);
   });

   return FunctionGradient(grd, g2, gstep);
}

const MnMachinePrecision &Numerical2PGradientCalculator::Precision() const
{
   // return global precision (set in transformation)
   return fTransformation.Precision();
}

unsigned int Numerical2PGradientCalculator::Ncycle() const
{
   // return number of cycles for gradient calculation (set in strategy object)
   return Strategy().GradientNCycles();
}

double Numerical2PGradientCalculator::StepTolerance() const
{
   // return gradient step tolerance (set in strategy object)
   return Strategy().GradientStepTolerance();
}

double Numerical2PGradientCalculator::GradTolerance() const
{
   // return gradient tolerance (set in strategy object)
   return Strategy().GradientTolerance();
}

} // namespace Minuit2

} // namespace ROOT
