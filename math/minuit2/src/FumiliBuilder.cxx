// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#include "Minuit2/FumiliBuilder.h"
#include "Minuit2/FumiliStandardMaximumLikelihoodFCN.h"
#include "Minuit2/GradientCalculator.h"
//#include "Minuit2/Numerical2PGradientCalculator.h"
#include "Minuit2/MinimumState.h"
#include "Minuit2/MinimumError.h"
#include "Minuit2/FunctionGradient.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnLineSearch.h"
#include "Minuit2/MinimumSeed.h"
#include "Minuit2/MnFcn.h"
#include "Minuit2/MnMachinePrecision.h"
#include "Minuit2/MnPosDef.h"
#include "Minuit2/MnParabolaPoint.h"
#include "Minuit2/LaSum.h"
#include "Minuit2/LaProd.h"
#include "Minuit2/MnStrategy.h"
#include "Minuit2/MnHesse.h"
#include "Minuit2/MnPrint.h"

namespace ROOT {

namespace Minuit2 {

double inner_product(const LAVector &, const LAVector &);

FunctionMinimum FumiliBuilder::Minimum(const MnFcn &fcn, const GradientCalculator &gc, const MinimumSeed &seed,
                                       const MnStrategy &strategy, unsigned int maxfcn, double edmval) const
{
   // top level function to find minimum from a given initial seed
   // iterate on a minimum search in case of first attempt is not successful

   MnPrint print("FumiliBuilder", PrintLevel());

   edmval *= 0.0001;
   // edmval *= 0.1; // use small factor for Fumili

   print.Debug("Convergence when edm <", edmval);

   if (seed.Parameters().Vec().size() == 0) {
      return FunctionMinimum(seed, fcn.Up());
   }

   //   double edm = Estimator().Estimate(seed.Gradient(), seed.Error());
   double edm = seed.State().Edm();

   FunctionMinimum min(seed, fcn.Up());

   if (edm < 0.) {
      print.Warn("Initial matrix not pos.def.");
      return min;
   }

   std::vector<MinimumState> result;
   //   result.reserve(1);
   result.reserve(8);

   result.push_back(seed.State());

   print.Info("Start iterating until Edm is <", edmval, '\n', "Initial state", MnPrint::Oneline(seed.State()));

   if (TraceIter())
      TraceIteration(result.size() - 1, result.back());

   // do actual iterations

   // try first with a maxfxn = 50% of maxfcn
   // Fumili in principle needs much less iterations
   int maxfcn_eff = int(0.5 * maxfcn);
   int ipass = 0;
   double edmprev = 1;

   do {

      min = Minimum(fcn, gc, seed, result, maxfcn_eff, edmval);

      // second time check for validity of function Minimum
      if (ipass > 0) {
         if (!min.IsValid()) {

            print.Warn("FunctionMinimum is invalid");

            return min;
         }
      }

      // resulting edm of minimization
      edm = result.back().Edm();

      print.Debug("Approximate Edm", edm, "npass", ipass);

      // call always Hesse (error matrix from Fumili is never accurate since is approximate)

      print.Debug("FumiliBuilder will verify convergence and Error matrix; "
                  "dcov",
                  min.Error().Dcovar());

      //       // recalculate the gradient using the numerical gradient calculator
      //       Numerical2PGradientCalculator ngc(fcn, min.Seed().Trafo(), strategy);
      //       FunctionGradient ng = ngc( min.State().Parameters() );
      //       MinimumState tmp( min.State().Parameters(), min.State().Error(), ng, min.State().Edm(),
      //       min.State().NFcn() );

      MinimumState st = MnHesse(strategy)(fcn, min.State(), min.Seed().Trafo(), maxfcn);
      result.push_back(st);

      print.Info("After Hessian");
      if (TraceIter())
         TraceIteration(result.size() - 1, result.back());

      // check edm
      edm = st.Edm();

      print.Debug("Edm", edm, "State", st);

      // break the loop if edm is NOT getting smaller
      if (ipass > 0 && edm >= edmprev) {
         print.Warn("Stop iterations, no improvements after Hesse; current Edm", edm, "previous value", edmprev);
         break;
      }
      if (edm > edmval) {
         print.Debug("Tolerance not sufficient, continue minimization; Edm", edm, "Requested", edmval);
      } else {
         // Case when edm < edmval after Heasse but min is flagged eith a  bad edm:
         // make then a new Function minimum since now edm is ok
         if (min.IsAboveMaxEdm()) {
            min = FunctionMinimum(min.Seed(), min.States(), min.Up());
            break;
         }
      }

      // end loop on iterations
      // ? need a maximum here (or max of function calls is enough ? )
      // continnue iteration (re-calculate function Minimum if edm IS NOT sufficient)
      // no need to check that hesse calculation is done (if isnot done edm is OK anyway)
      // count the pass to exit second time when function Minimum is invalid
      // increase by 20% maxfcn for doing some more tests
      if (ipass == 0)
         maxfcn_eff = maxfcn;
      ipass++;
      edmprev = edm;

   } while (edm > edmval);

   // Add latest state (Hessian calculation)
   min.Add(result.back());

   return min;
}

FunctionMinimum FumiliBuilder::Minimum(const MnFcn &fcn, const GradientCalculator &gc, const MinimumSeed &seed,
                                       std::vector<MinimumState> &result, unsigned int maxfcn, double edmval) const
{

   // function performing the minimum searches using the FUMILI algorithm
   // after the modification when I iterate on this functions, so it can be called many times,
   //  the seed is used here only to get precision and construct the returned FunctionMinimum object

   /*
    Three options were possible:

    1) create two parallel and completely separate hierarchies, in which case
    the FumiliMinimizer would NOT inherit from ModularFunctionMinimizer,
    FumiliBuilder would not inherit from MinimumBuilder etc

    2) Use the inheritance (base classes of ModularFunctionMinimizer,
                            MinimumBuilder etc), but recreate the member functions Minimize() and
    Minimum() respectively (naming them for example minimize2() and
                            minimum2()) so that they can take FumiliFCNBase as Parameter instead FCNBase
    (otherwise one wouldn't be able to call the Fumili-specific methods).

    3) Cast in the daughter classes derived from ModularFunctionMinimizer,
    MinimumBuilder.

    The first two would mean to duplicate all the functionality already existent,
    which is a very bad practice and Error-prone. The third one is the most
    elegant and effective solution, where the only constraint is that the user
    must know that they have to pass a subclass of FumiliFCNBase to the FumiliMinimizer
    and not just a subclass of FCNBase.
    BTW, the first two solutions would have meant to recreate also a parallel
    structure for MnFcn...
    **/
   //  const FumiliFCNBase* tmpfcn =  dynamic_cast<const FumiliFCNBase*>(&(fcn.Fcn()));

   const MnMachinePrecision &prec = seed.Precision();

   const MinimumState &initialState = result.back();

   double edm = initialState.Edm();

   MnPrint print("FumiliBuilder");

   print.Debug("Initial State:", "\n  Parameter", initialState.Vec(), "\n  Gradient", initialState.Gradient().Vec(),
               "\n  Inv Hessian", initialState.Error().InvHessian(), "\n  edm", initialState.Edm(), "\n  maxfcn",
               maxfcn, "\n  tolerance", edmval);

   // iterate until edm is small enough or max # of iterations reached
   edm *= (1. + 3. * initialState.Error().Dcovar());
   MnAlgebraicVector step(initialState.Gradient().Vec().size());

   // initial lambda Value
   double lambda = 0.001;

   do {

      //     const MinimumState& s0 = result.back();
      MinimumState s0 = result.back();

      step = -1. * s0.Error().InvHessian() * s0.Gradient().Vec();

      print.Debug("Iteration -", result.size(), "\n  Fval", s0.Fval(), "numOfCall", fcn.NumOfCalls(),
                  "\n  Internal Parameter values", s0.Vec(), "\n  Newton step", step);

      double gdel = inner_product(step, s0.Gradient().Grad());
      if (gdel > 0.) {
         print.Warn("Matrix not pos.def, gdel =", gdel, " > 0");

         MnPosDef psdf;
         s0 = psdf(s0, prec);
         step = -1. * s0.Error().InvHessian() * s0.Gradient().Vec();
         gdel = inner_product(step, s0.Gradient().Grad());

         print.Warn("After correction, gdel =", gdel);

         if (gdel > 0.) {
            result.push_back(s0);
            return FunctionMinimum(seed, result, fcn.Up());
         }
      }

      //     MnParabolaPoint pp = lsearch(fcn, s0.Parameters(), step, gdel, prec);

      //     if(std::fabs(pp.Y() - s0.Fval()) < prec.Eps()) {
      //       std::cout<<"FumiliBuilder: no improvement"<<std::endl;
      //       break; //no improvement
      //     }

      //     MinimumParameters p(s0.Vec() + pp.X()*step, pp.Y());

      // if taking a full step

      // take a full step

      MinimumParameters p(s0.Vec() + step, fcn(s0.Vec() + step));

      // check that taking the full step does not deteriorate minimum
      // in that case do a line search
      if (p.Fval() >= s0.Fval()) {
         MnLineSearch lsearch;
         MnParabolaPoint pp = lsearch(fcn, s0.Parameters(), step, gdel, prec);

         if (std::fabs(pp.Y() - s0.Fval()) < prec.Eps()) {
            // std::cout<<"FumiliBuilder: no improvement"<<std::endl;
            break; // no improvement
         }
         p = MinimumParameters(s0.Vec() + pp.X() * step, pp.Y());
      }

      print.Debug("Before Gradient", fcn.NumOfCalls());

      FunctionGradient g = gc(p, s0.Gradient());

      print.Debug("After Gradient", fcn.NumOfCalls());

      // FunctionGradient g = gc(s0.Parameters(), s0.Gradient());

      // move Error updator after Gradient since the Value is cached inside

      MinimumError e = ErrorUpdator().Update(s0, p, gc, lambda);

      edm = Estimator().Estimate(g, s0.Error());

      print.Debug("Updated new point:", "\n  FVAL     ", p.Fval(), "\n  Parameter", p.Vec(), "\n  Gradient", g.Vec(),
                  "\n  InvHessian", e.Matrix(), "\n  Hessian", e.Hessian(), "\n  Edm", edm);

      if (edm < 0.) {
         print.Warn("Matrix not pos.def., Edm < 0");

         MnPosDef psdf;
         s0 = psdf(s0, prec);
         edm = Estimator().Estimate(g, s0.Error());
         if (edm < 0.) {
            result.push_back(s0);
            if (TraceIter())
               TraceIteration(result.size() - 1, result.back());
            return FunctionMinimum(seed, result, fcn.Up());
         }
      }

      // check lambda according to step
      if (p.Fval() < s0.Fval())
         // fcn is decreasing along the step
         lambda *= 0.1;
      else {
         lambda *= 10;
         // if close to minimum stop to avoid oscillations around minimum value
         if (edm < 0.1)
            break;
      }

      print.Debug("finish iteration -", result.size(), "lambda =", lambda, "f1 =", p.Fval(), "f0 =", s0.Fval(),
                  "num of calls =", fcn.NumOfCalls(), "edm =", edm);

      result.push_back(MinimumState(p, e, g, edm, fcn.NumOfCalls()));
      if (TraceIter())
         TraceIteration(result.size() - 1, result.back());

      print.Info(MnPrint::Oneline(result.back(), result.size()));

      edm *= (1. + 3. * e.Dcovar());

   } while (edm > edmval && fcn.NumOfCalls() < maxfcn);

   if (fcn.NumOfCalls() >= maxfcn) {
      print.Warn("Call limit exceeded");

      return FunctionMinimum(seed, result, fcn.Up(), FunctionMinimum::MnReachedCallLimit());
   }

   if (edm > edmval) {
      if (edm < std::fabs(prec.Eps2() * result.back().Fval())) {
         print.Warn("Machine accuracy limits further improvement");

         return FunctionMinimum(seed, result, fcn.Up());
      } else if (edm < 10 * edmval) {
         return FunctionMinimum(seed, result, fcn.Up());
      } else {

         print.Warn("No convergence; Edm", edm, "is above tolerance", 10 * edmval);

         return FunctionMinimum(seed, result, fcn.Up(), FunctionMinimum::MnAboveMaxEdm());
      }
   }
   //   std::cout<<"result.back().Error().Dcovar()= "<<result.back().Error().Dcovar()<<std::endl;

   print.Debug("Exiting successfully", "Ncalls", fcn.NumOfCalls(), "FCN", result.back().Fval(), "Edm", edm, "Requested",
               edmval);

   return FunctionMinimum(seed, result, fcn.Up());
}

} // namespace Minuit2

} // namespace ROOT
