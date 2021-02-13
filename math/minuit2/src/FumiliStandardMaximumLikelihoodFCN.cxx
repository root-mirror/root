// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#include "Minuit2/FumiliStandardMaximumLikelihoodFCN.h"

#include <vector>
#include <cmath>
#include <limits>

namespace ROOT {

namespace Minuit2 {

std::vector<double> FumiliStandardMaximumLikelihoodFCN::Elements(const std::vector<double> &par) const
{

   //   calculate likelihood element f(i) = pdf(x(i))
   std::vector<double> result;
   double tmp1 = 0.0;
   unsigned int fPositionsSize = fPositions.size();

   for (unsigned int i = 0; i < fPositionsSize; i++) {

      const std::vector<double> &currentPosition = fPositions[i];

      // The commented line is the object-oriented way to do it
      // but it is faster to do a single function call...
      //(*(this->getModelFunction())).SetParameters(par);
      tmp1 = (*(this->ModelFunction()))(par, currentPosition);

      // std::cout << " i = " << i << "  " << currentPosition[0] << "  " << tmp1 << std::endl;

      result.push_back(tmp1);
   }

   return result;
}

const std::vector<double> &FumiliStandardMaximumLikelihoodFCN::GetMeasurement(int index) const
{
   // Return x(i).
   return fPositions[index];
}

int FumiliStandardMaximumLikelihoodFCN::GetNumberOfMeasurements() const
{
   // return size of positions (coordinates).
   return fPositions.size();
}

void FumiliStandardMaximumLikelihoodFCN::EvaluateAll(const std::vector<double> &par)
{
   // Evaluate in one loop likelihood value, gradient and hessian

   const double minDouble = 8.0 * std::numeric_limits<double>::min();
   const double minDouble2 = std::sqrt(minDouble);
   const double maxDouble2 = 1.0 / minDouble2;
   // loop on the measurements

   int nmeas = GetNumberOfMeasurements();
   std::vector<double> &grad = Gradient();
   std::vector<double> &h = Hessian();
   int npar = par.size();
   double logLikelihood = 0;
   grad.resize(npar);
   h.resize(static_cast<unsigned int>(0.5 * npar * (npar + 1)));
   grad.assign(npar, 0.0);
   h.assign(static_cast<unsigned int>(0.5 * npar * (npar + 1)), 0.0);

   const ParametricFunction &modelFunc = *ModelFunction();

   for (int i = 0; i < nmeas; ++i) {

      // work for one-dimensional points
      const std::vector<double> &currentPosition = fPositions[i];
      modelFunc.SetParameters(currentPosition);
      double fval = modelFunc(par);
      if (fval < minDouble)
         fval = minDouble; // to avoid getting infinity and nan's
      logLikelihood -= std::log(fval);
      double invFval = 1.0 / fval;
      // this method should return a reference
      std::vector<double> mfg = modelFunc.GetGradient(par);

      // calc derivatives

      for (int j = 0; j < npar; ++j) {
         if (std::fabs(mfg[j]) < minDouble) {
            //  std::cout << "SMALL values: grad =  " << mfg[j] << "  "  << minDouble << " f(x) = " << fval
            //    << " params " << j << " p0 = " << par[0] << " p1 = " << par[1] <<  std::endl;
            if (mfg[j] < 0)
               mfg[j] = -minDouble;
            else
               mfg[j] = minDouble;
         }

         double dfj = invFval * mfg[j];
         // to avoid summing infinite and nan later when calculating the Hessian
         if (std::fabs(dfj) > maxDouble2) {
            if (dfj > 0)
               dfj = maxDouble2;
            else
               dfj = -maxDouble2;
         }

         grad[j] -= dfj;
         //       if ( ! ( dfj > 0) && ! ( dfj <= 0 ) )
         // std::cout << " nan : dfj = " << dfj << " fval =  " << fval << " invF = " << invFval << " grad = " << mfg[j]
         // << " par[j] = " << par[j] << std::endl;

         // std::cout << " x = "  << currentPosition[0] <<  " par[j] = " << par[j] << " : dfj = " << dfj << " fval =  "
         // << fval << " invF = " << invFval << " grad = " << mfg[j] << " deriv = " << grad[j] << std::endl;

         // in second derivative use Fumili approximation neglecting the term containing the
         // second derivatives of the model function
         for (int k = j; k < npar; ++k) {
            int idx = j + k * (k + 1) / 2;
            if (std::fabs(mfg[k]) < minDouble) {
               if (mfg[k] < 0)
                  mfg[k] = -minDouble;
               else
                  mfg[k] = minDouble;
            }

            double dfk = invFval * mfg[k];
            // avoid that dfk*dfj are one small and one huge so I get a nan
            // to avoid summing infinite and nan later when calculating the Hessian
            if (std::fabs(dfk) > maxDouble2) {
               if (dfk > 0)
                  dfk = maxDouble2;
               else
                  dfk = -maxDouble2;
            }

            h[idx] += dfj * dfk;
            // if ( ( ! ( h[idx] > 0) && ! ( h[idx] <= 0 ) ) )
            //   std::cout << " nan : dfj = " << dfj << " fval =  " << fval << " invF = " << invFval << " gradj = " <<
            //   mfg[j]
            //     << " dfk = " << dfk << " gradk =  "<< mfg[k]  << " hess_jk = " << h[idx] << " par[k] = " << par[k] <<
            //     std::endl;
         }

      } // end param loop

   } // end points loop

   //   std::cout <<"\nEVALUATED GRADIENT and HESSIAN " << std::endl;
   //   for (int j = 0; j < npar; ++j) {
   //     std::cout << " j = " << j << " grad = " << grad[j] << std::endl;
   //     for (int k = j; k < npar; ++k) {
   //       std::cout << " k = " << k << " hess = " << Hessian(j,k) << "  " << h[ j + k*(k+1)/2]  << std::endl;
   //     }
   //   }

   // set Value in base class
   SetFCNValue(logLikelihood);
}

} // namespace Minuit2

} // namespace ROOT
