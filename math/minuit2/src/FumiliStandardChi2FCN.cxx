// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#include "Minuit2/FumiliStandardChi2FCN.h"

#include <vector>
#include <cmath>

namespace ROOT {

namespace Minuit2 {

std::vector<double> FumiliStandardChi2FCN::Elements(const std::vector<double> &par) const
{
   // Calculate the f(i) contribution to the Chi2. Chi2 = Sum[f(i)**2]

   std::vector<double> result;
   double tmp1 = 0.0;
   unsigned int fPositionsSize = fPositions.size();

   for (unsigned int i = 0; i < fPositionsSize; i++) {

      const std::vector<double> &currentPosition = fPositions[i];

      // The commented line is the object-oriented way to do it
      // but it is faster to do a single function call...
      //(*(this->getModelFunction())).SetParameters(par);
      tmp1 = (*(this->ModelFunction()))(par, currentPosition) - fMeasurements[i];

      result.push_back(tmp1 * fInvErrors[i]);

      // std::cout << "element " << i << "  " << (*(this->getModelFunction()))(par, currentPosition) << "  " <<
      // fMeasurements[i] << "  " << result[i] << std::endl;
   }

   return result;
}

const std::vector<double> &FumiliStandardChi2FCN::GetMeasurement(int index) const
{
   // Return the coordinate (position) values.
   return fPositions[index];
}

int FumiliStandardChi2FCN::GetNumberOfMeasurements() const
{
   // Return size
   return fPositions.size();
}

void FumiliStandardChi2FCN::EvaluateAll(const std::vector<double> &par)
{
   // Evaluate chi2 value, gradient and hessian all in a single
   // loop on the measurements

   int nmeas = GetNumberOfMeasurements();
   std::vector<double> &grad = Gradient();
   std::vector<double> &h = Hessian();
   int npar = par.size();
   double chi2 = 0;
   grad.resize(npar);
   h.resize(static_cast<unsigned int>(0.5 * npar * (npar + 1)));
   // reset Elements
   grad.assign(npar, 0.0);
   h.assign(static_cast<unsigned int>(0.5 * npar * (npar + 1)), 0.0);

   const ParametricFunction &modelFunc = *ModelFunction();

   for (int i = 0; i < nmeas; ++i) {

      // work for multi-dimensional points
      const std::vector<double> &currentPosition = fPositions[i];
      modelFunc.SetParameters(currentPosition);
      double invError = fInvErrors[i];
      double fval = modelFunc(par);

      double element = (fval - fMeasurements[i]) * invError;
      chi2 += element * element;

      // calc derivatives

      // this method should return a reference
      std::vector<double> mfg = modelFunc.GetGradient(par);

      // grad is derivative of chi2 w.r.t to parameters
      for (int j = 0; j < npar; ++j) {
         double dfj = invError * mfg[j];
         grad[j] += 2.0 * element * dfj;

         // in second derivative use Fumili approximation neglecting the term containing the
         // second derivatives of the model function
         for (int k = j; k < npar; ++k) {
            int idx = j + k * (k + 1) / 2;
            h[idx] += 2.0 * dfj * invError * mfg[k];
         }

      } // end param loop

   } // end points loop

   // set Value in base class
   SetFCNValue(chi2);
}

} // namespace Minuit2

} // namespace ROOT
