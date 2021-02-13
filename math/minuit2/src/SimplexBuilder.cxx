// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#include "Minuit2/SimplexBuilder.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnFcn.h"
#include "Minuit2/MinimumSeed.h"
#include "Minuit2/SimplexParameters.h"
#include "Minuit2/MinimumState.h"
#include "Minuit2/MnPrint.h"

namespace ROOT {

namespace Minuit2 {

class GradientCalculator;
class MnStrategy;
FunctionMinimum SimplexBuilder::Minimum(const MnFcn &mfcn, const GradientCalculator &, const MinimumSeed &seed,
                                        const MnStrategy &, unsigned int maxfcn, double minedm) const
{
   // find the minimum using the Simplex method of Nelder and Mead (does not use function gradient)
   // method to find initial simplex is slightly different than in the orginal Fortran
   // Minuit since has not been proofed that one to be better

   // synchronize print levels
   MnPrint print("SimplexBuilder", PrintLevel());

   print.Debug("Running with maxfcn", maxfcn, "minedm", minedm);

   const MnMachinePrecision &prec = seed.Precision();
   MnAlgebraicVector x = seed.Parameters().Vec();
   MnAlgebraicVector step = 10. * seed.Gradient().Gstep();

   unsigned int n = x.size();
   double wg = 1. / double(n);
   double alpha = 1., beta = 0.5, gamma = 2., rhomin = 4., rhomax = 8.;
   double rho1 = 1. + alpha;
   // double rho2 = rho1 + alpha*gamma;
   // change proposed by david sachs (fnal)
   double rho2 = 1. + alpha * gamma;

   std::vector<std::pair<double, MnAlgebraicVector>> simpl;
   simpl.reserve(n + 1);
   simpl.push_back(std::pair<double, MnAlgebraicVector>(seed.Fval(), x));

   unsigned int jl = 0, jh = 0;
   double amin = seed.Fval(), aming = seed.Fval();

   for (unsigned int i = 0; i < n; i++) {
      double dmin = 8. * prec.Eps2() * (std::fabs(x(i)) + prec.Eps2());
      if (step(i) < dmin)
         step(i) = dmin;
      x(i) += step(i);
      double tmp = mfcn(x);
      if (tmp < amin) {
         amin = tmp;
         jl = i + 1;
      }
      if (tmp > aming) {
         aming = tmp;
         jh = i + 1;
      }
      simpl.push_back(std::pair<double, MnAlgebraicVector>(tmp, x));
      x(i) -= step(i);
   }
   SimplexParameters simplex(simpl, jh, jl);

   print.Debug([&](std::ostream &os) {
      os << "Initial parameters - min  " << jl << "  " << amin << " max " << jh << "  " << aming << '\n';
      for (unsigned int i = 0; i < simplex.Simplex().size(); ++i)
         os << " i = " << i << " x = " << simplex(i).second << " fval(x) = " << simplex(i).first << '\n';
   });

   double edmPrev = simplex.Edm();
   int niterations = 0;
   do {
      jl = simplex.Jl();
      jh = simplex.Jh();
      amin = simplex(jl).first;
      edmPrev = simplex.Edm();

      print.Debug("iteration: edm =", simplex.Edm(), '\n', "--> Min Param is", jl, "pmin", simplex(jl).second,
                  "f(pmin)", amin, '\n', "--> Max param is", jh, simplex(jh).first);

      //     std::cout << "ALL SIMPLEX PARAMETERS: "<< std::endl;
      //     for (unsigned int i = 0; i < simplex.Simplex().size(); ++i)
      //       std::cout << " i = " << i << " x = " << simplex(i).second << " fval(x) = " << simplex(i).first <<
      //       std::endl;

      // trace the iterations (need to create a MinimunState although errors and gradient are not existing)
      if (TraceIter())
         TraceIteration(niterations, MinimumState(MinimumParameters(simplex(jl).second, simplex(jl).first),
                                                  simplex.Edm(), mfcn.NumOfCalls()));
      print.Info(MnPrint::Oneline(simplex(jl).first, simplex.Edm(), mfcn.NumOfCalls(), niterations));
      niterations++;

      MnAlgebraicVector pbar(n);
      for (unsigned int i = 0; i < n + 1; i++) {
         if (i == jh)
            continue;
         pbar += (wg * simplex(i).second);
      }

      MnAlgebraicVector pstar = (1. + alpha) * pbar - alpha * simplex(jh).second;
      double ystar = mfcn(pstar);

      print.Debug("pbar", pbar, "pstar", pstar, "f(pstar)", ystar);

      if (ystar > amin) {
         if (ystar < simplex(jh).first) {
            simplex.Update(ystar, pstar);
            if (jh != simplex.Jh())
               continue;
         }
         MnAlgebraicVector pstst = beta * simplex(jh).second + (1. - beta) * pbar;
         double ystst = mfcn(pstst);

         print.Debug("Reduced simplex pstst", pstst, "f(pstst)", ystst);

         if (ystst > simplex(jh).first)
            break;
         simplex.Update(ystst, pstst);
         continue;
      }

      MnAlgebraicVector pstst = gamma * pstar + (1. - gamma) * pbar;
      double ystst = mfcn(pstst);

      print.Debug("pstst", pstst, "f(pstst)", ystst);

      double y1 = (ystar - simplex(jh).first) * rho2;
      double y2 = (ystst - simplex(jh).first) * rho1;
      double rho = 0.5 * (rho2 * y1 - rho1 * y2) / (y1 - y2);
      if (rho < rhomin) {
         if (ystst < simplex(jl).first)
            simplex.Update(ystst, pstst);
         else
            simplex.Update(ystar, pstar);
         continue;
      }
      if (rho > rhomax)
         rho = rhomax;
      MnAlgebraicVector prho = rho * pbar + (1. - rho) * simplex(jh).second;
      double yrho = mfcn(prho);

      print.Debug("prho", prho, "f(prho)", yrho);

      if (yrho < simplex(jl).first && yrho < ystst) {
         simplex.Update(yrho, prho);
         continue;
      }
      if (ystst < simplex(jl).first) {
         simplex.Update(ystst, pstst);
         continue;
      }
      if (yrho > simplex(jl).first) {
         if (ystst < simplex(jl).first)
            simplex.Update(ystst, pstst);
         else
            simplex.Update(ystar, pstar);
         continue;
      }
      if (ystar > simplex(jh).first) {
         pstst = beta * simplex(jh).second + (1. - beta) * pbar;
         ystst = mfcn(pstst);
         if (ystst > simplex(jh).first)
            break;
         simplex.Update(ystst, pstst);
      }

      print.Debug("End loop : Edm", simplex.Edm(), "pstst", pstst, "f(pstst)", ystst);

   } while ((simplex.Edm() > minedm || edmPrev > minedm) && mfcn.NumOfCalls() < maxfcn);

   jl = simplex.Jl();
   jh = simplex.Jh();
   amin = simplex(jl).first;

   MnAlgebraicVector pbar(n);
   for (unsigned int i = 0; i < n + 1; i++) {
      if (i == jh)
         continue;
      pbar += (wg * simplex(i).second);
   }
   double ybar = mfcn(pbar);
   if (ybar < amin)
      simplex.Update(ybar, pbar);
   else {
      pbar = simplex(jl).second;
      ybar = simplex(jl).first;
   }

   MnAlgebraicVector dirin = simplex.Dirin();
   //   Scale to sigmas on parameters werr^2 = dirin^2 * (up/edm)
   dirin *= std::sqrt(mfcn.Up() / simplex.Edm());

   print.Debug("End simplex edm =", simplex.Edm(), "pbar =", pbar, "f(p) =", ybar);

   MinimumState st(MinimumParameters(pbar, dirin, ybar), simplex.Edm(), mfcn.NumOfCalls());

   print.Info("Final iteration", MnPrint::Oneline(st));

   if (TraceIter())
      TraceIteration(niterations, st);

   if (mfcn.NumOfCalls() > maxfcn) {
      print.Warn("Simplex did not converge, #fcn calls exhausted");

      return FunctionMinimum(seed, std::vector<MinimumState>(1, st), mfcn.Up(), FunctionMinimum::MnReachedCallLimit());
   }
   if (simplex.Edm() > minedm) {
      print.Warn("Simplex did not converge, edm > minedm");

      return FunctionMinimum(seed, std::vector<MinimumState>(1, st), mfcn.Up(), FunctionMinimum::MnAboveMaxEdm());
   }

   return FunctionMinimum(seed, std::vector<MinimumState>(1, st), mfcn.Up());
}

} // namespace Minuit2

} // namespace ROOT
