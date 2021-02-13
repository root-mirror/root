// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#include "Minuit2/VariableMetricBuilder.h"
#include "Minuit2/GradientCalculator.h"
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

#include <cmath>
#include <cassert>

namespace ROOT {

namespace Minuit2 {

double inner_product(const LAVector &, const LAVector &);

void VariableMetricBuilder::AddResult(std::vector<MinimumState> &result, const MinimumState &state) const
{
   // // if (!store) store = StorageLevel();
   // // store |= (result.size() == 0);
   // if (store)
   result.push_back(state);
   //  else {
   //     result.back() = state;
   //  }
   if (TraceIter())
      TraceIteration(result.size() - 1, result.back());
   else {
      MnPrint print("VariableMetricBuilder", PrintLevel());
      print.Info(MnPrint::Oneline(result.back(), result.size() - 1));
   }
}

FunctionMinimum VariableMetricBuilder::Minimum(const MnFcn &fcn, const GradientCalculator &gc, const MinimumSeed &seed,
                                               const MnStrategy &strategy, unsigned int maxfcn, double edmval) const
{
   MnPrint print("VariableMetricBuilder", PrintLevel());

   // top level function to find minimum from a given initial seed
   // iterate on a minimum search in case of first attempt is not successful

   // to be consistent with F77 Minuit
   // in Minuit2 edm is correct and is ~ a factor of 2 smaller than F77Minuit
   // There are also a check for convergence if (edm < 0.1 edmval for exiting the loop)
   // LM: change factor to 2E-3 to be consistent with F77Minuit
   edmval *= 0.002;

   // set global printlevel to the local one so all calls to MN_INFO_MSG can be controlled in the same way
   // at exit of this function the BuilderPrintLevelConf object is destructed and automatically the
   // previous level will be restored

   //   double edm = Estimator().Estimate(seed.Gradient(), seed.Error());
   double edm = seed.State().Edm();

   FunctionMinimum min(seed, fcn.Up());

   if (seed.Parameters().Vec().size() == 0) {
      print.Warn("No free parameters.");
      return min;
   }

   if (!seed.IsValid()) {
     print.Error("Minimum seed invalid.");
     return min;
   }

   if (edm < 0.) {
      print.Error("Initial matrix not pos.def.");

      // assert(!seed.Error().IsPosDef());
      return min;
   }

   std::vector<MinimumState> result;
   if (StorageLevel() > 0)
      result.reserve(10);
   else
      result.reserve(2);

   // do actual iterations
   print.Info("Start iterating until Edm is <", edmval, "with call limit =", maxfcn);

   AddResult(result, seed.State());

   // try first with a maxfxn = 80% of maxfcn
   int maxfcn_eff = maxfcn;
   int ipass = 0;
   bool iterate = false;
   bool hessianComputed = false;

   do {

      iterate = false;
      hessianComputed = false;

      print.Debug(ipass > 0 ? "Continue" : "Start", "iterating...");

      min = Minimum(fcn, gc, seed, result, maxfcn_eff, edmval);

      // if max function call reached exits
      if (min.HasReachedCallLimit()) {
         print.Warn("FunctionMinimum is invalid, reached function call limit");
         return min;
      }

      // second time check for validity of function Minimum
      if (ipass > 0) {
         if (!min.IsValid()) {
            print.Warn("FunctionMinimum is invalid after second try");
            return min;
         }
      }

      // resulting edm of minimization
      edm = result.back().Edm();
      // need to correct again for Dcovar: edm *= (1. + 3. * e.Dcovar()) ???

      if ((strategy.Strategy() == 2) || (strategy.Strategy() == 1 && min.Error().Dcovar() > 0.05)) {

         print.Debug("MnMigrad will verify convergence and Error matrix; dcov =", min.Error().Dcovar());

         MinimumState st = MnHesse(strategy)(fcn, min.State(), min.Seed().Trafo(), maxfcn);

         hessianComputed = true;

         print.Info("After Hessian");

         AddResult(result, st);

         if (!st.IsValid()) {
           print.Warn("Invalid Hessian - exit the minimization");
           break;
         }

         // check new edm
         edm = st.Edm();

         print.Debug("New Edm", edm, "Requested", edmval);

         if (edm > edmval) {
            // be careful with machine precision and avoid too small edm
            double machineLimit = std::fabs(seed.Precision().Eps2() * result.back().Fval());
            if (edm >= machineLimit) {
               iterate = true;

               print.Info("Tolerance not sufficient, continue minimization; "
                          "Edm",
                          edm, "Required", edmval);
            } else {
               print.Warn("Reached machine accuracy limit; Edm", edm, "is smaller than machine limit", machineLimit,
                          "while", edmval, "was requested");
            }
         }
      }

      // end loop on iterations
      // ? need a maximum here (or max of function calls is enough ? )
      // continnue iteration (re-calculate function Minimum if edm IS NOT sufficient)
      // no need to check that hesse calculation is done (if isnot done edm is OK anyway)
      // count the pass to exit second time when function Minimum is invalid
      // increase by 20% maxfcn for doing some more tests
      if (ipass == 0)
         maxfcn_eff = int(maxfcn * 1.3);
      ipass++;
   } while (iterate);

   // Add latest state (Hessian calculation)
   // and check edm (add a factor of 10 in tolerance )
   if (edm > 10 * edmval) {
      min.Add(result.back(), FunctionMinimum::MnAboveMaxEdm());
      print.Warn("No convergence; Edm", edm, "is above tolerance", 10 * edmval);
   } else {
      // check if minimum has edm above max before
      if (min.IsAboveMaxEdm()) {
         print.Info("Edm has been re-computed after Hesse; Edm", edm, "is now within tolerance");
      }
      if (hessianComputed) min.Add(result.back());
   }

   print.Debug("Minimum found", min);

   return min;
}

FunctionMinimum VariableMetricBuilder::Minimum(const MnFcn &fcn, const GradientCalculator &gc, const MinimumSeed &seed,
                                               std::vector<MinimumState> &result, unsigned int maxfcn,
                                               double edmval) const
{
   // function performing the minimum searches using the Variable Metric  algorithm (MIGRAD)
   // perform first a line search in the - Vg direction and then update using the Davidon formula (Davidon Error
   // updator) stop when edm reached is less than required (edmval)

   // after the modification when I iterate on this functions, so it can be called many times,
   //  the seed is used here only to get precision and construct the returned FunctionMinimum object

   MnPrint print("VariableMetricBuilder", PrintLevel());

   const MnMachinePrecision &prec = seed.Precision();

   //   result.push_back(MinimumState(seed.Parameters(), seed.Error(), seed.Gradient(), edm, fcn.NumOfCalls()));
   const MinimumState &initialState = result.back();

   double edm = initialState.Edm();

   print.Debug("Initial State:", "\n  Parameter:", initialState.Vec(), "\n  Gradient:", initialState.Gradient().Vec(),
               "\n  InvHessian:", initialState.Error().InvHessian(), "\n  Edm:", initialState.Edm());

   // iterate until edm is small enough or max # of iterations reached
   edm *= (1. + 3. * initialState.Error().Dcovar());
   MnLineSearch lsearch;
   MnAlgebraicVector step(initialState.Gradient().Vec().size());
   // keep also prevStep
   MnAlgebraicVector prevStep(initialState.Gradient().Vec().size());

   MinimumState s0 = result.back();

   do {

      // MinimumState s0 = result.back();

      step = -1. * s0.Error().InvHessian() * s0.Gradient().Vec();

      print.Debug("Iteration", result.size(), "Fval", s0.Fval(), "numOfCall", fcn.NumOfCalls(),
                  "\n  Internal parameters", s0.Vec(), "\n  Newton step", step);

      // check if derivatives are not zero
      if (inner_product(s0.Gradient().Vec(), s0.Gradient().Vec()) <= 0) {
         print.Debug("all derivatives are zero - return current status");
         break;
      }

      double gdel = inner_product(step, s0.Gradient().Grad());

      if (gdel > 0.) {
         print.Warn("Matrix not pos.def, gdel =", gdel, "> 0");

         MnPosDef psdf;
         s0 = psdf(s0, prec);
         step = -1. * s0.Error().InvHessian() * s0.Gradient().Vec();
         // #ifdef DEBUG
         //       std::cout << "After MnPosdef - Error  " << s0.Error().InvHessian() << " Gradient " <<
         //       s0.Gradient().Vec() << " step " << step << std::endl;
         // #endif
         gdel = inner_product(step, s0.Gradient().Grad());

         print.Warn("gdel =", gdel);

         if (gdel > 0.) {
            AddResult(result, s0);

            return FunctionMinimum(seed, result, fcn.Up());
         }
      }

      MnParabolaPoint pp = lsearch(fcn, s0.Parameters(), step, gdel, prec);

      // <= needed for case 0 <= 0
      if (std::fabs(pp.Y() - s0.Fval()) <= std::fabs(s0.Fval()) * prec.Eps()) {

         print.Warn("No improvement in line search");

         // no improvement exit   (is it really needed LM ? in vers. 1.22 tried alternative )
         // add new state when only fcn changes
         if (result.size() <= 1)
            AddResult(result, MinimumState(s0.Parameters(), s0.Error(), s0.Gradient(), s0.Edm(), fcn.NumOfCalls()));
         else
            // no need to re-store the state
            AddResult(result, MinimumState(pp.Y(), s0.Edm(), fcn.NumOfCalls()));

         break;
      }

      print.Debug("Result after line search :", "\n  x =", pp.X(), "\n  Old Fval =", s0.Fval(),
                  "\n  New Fval =", pp.Y(), "\n  NFcalls =", fcn.NumOfCalls());

      MinimumParameters p(s0.Vec() + pp.X() * step, pp.Y());

      FunctionGradient g = gc(p, s0.Gradient());

      edm = Estimator().Estimate(g, s0.Error());

      if (std::isnan(edm)) {
         print.Warn("Edm is NaN; stop iterations");
         AddResult(result, s0);
         return FunctionMinimum(seed, result, fcn.Up());
      }

      if (edm < 0.) {
         print.Warn("Matrix not pos.def., try to make pos.def.");

         MnPosDef psdf;
         s0 = psdf(s0, prec);
         edm = Estimator().Estimate(g, s0.Error());
         if (edm < 0.) {
            print.Warn("Matrix still not pos.def.; stop iterations");

            AddResult(result, s0);

            return FunctionMinimum(seed, result, fcn.Up());
         }
      }
      MinimumError e = ErrorUpdator().Update(s0, p, g);

      // avoid print Hessian that will invert the matrix
      print.Debug("Updated new point:", "\n  Parameter:", p.Vec(), "\n  Gradient:", g.Vec(),
                  "\n  InvHessian:", e.Matrix(), "\n  Edm:", edm);

      // update the state
      s0 = MinimumState(p, e, g, edm, fcn.NumOfCalls());
      if (StorageLevel() || result.size() <= 1)
         AddResult(result, s0);
      else
         // use a reduced state for not-final iterations
         AddResult(result, MinimumState(p.Fval(), edm, fcn.NumOfCalls()));

      // correct edm
      edm *= (1. + 3. * e.Dcovar());

      print.Debug("Dcovar =", e.Dcovar(), "\tCorrected edm =", edm);

   } while (edm > edmval && fcn.NumOfCalls() < maxfcn); // end of iteration loop

   // save last result in case of no complete final states
   // when the result is filled above (reduced storage) the resulting state will not be valid
   // since they will not have parameter values and error
   // the line above will fill as last element a valid state
   if (!result.back().IsValid())
      result.back() = s0;

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
         print.Warn("Iterations finish without convergence; Edm", edm, "Requested", edmval);

         return FunctionMinimum(seed, result, fcn.Up(), FunctionMinimum::MnAboveMaxEdm());
      }
   }
   //   std::cout<<"result.back().Error().Dcovar()= "<<result.back().Error().Dcovar()<<std::endl;

   print.Debug("Exiting successfully;", "Ncalls", fcn.NumOfCalls(), "FCN", result.back().Fval(), "Edm", edm,
               "Requested", edmval);

   return FunctionMinimum(seed, result, fcn.Up());
}

} // namespace Minuit2

} // namespace ROOT
