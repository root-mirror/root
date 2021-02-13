// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#include "GaussFcn2.h"
#include "GaussFunction.h"

#include <cassert>

namespace ROOT {

namespace Minuit2 {

double GaussFcn2::operator()(const std::vector<double> &par) const
{

   assert(par.size() == 6);

   GaussFunction gauss1 = GaussFunction(par[0], par[1], par[2]);
   GaussFunction gauss2 = GaussFunction(par[3], par[4], par[5]);

   double chi2 = 0.;
   int nmeas = fMeasurements.size();
   for (int n = 0; n < nmeas; n++) {
      chi2 += ((gauss1(fPositions[n]) + gauss2(fPositions[n]) - fMeasurements[n]) *
               (gauss1(fPositions[n]) + gauss2(fPositions[n]) - fMeasurements[n]) / fMVariances[n]);
   }

   return chi2;
}

void GaussFcn2::Init()
{

   // calculate initial Value of chi2

   int nmeas = fMeasurements.size();
   double x = 0.;
   double x2 = 0.;
   double norm = 0.;
   double dx = fPositions[1] - fPositions[0];
   double c = 0.;
   for (int i = 0; i < nmeas; i++) {
      norm += fMeasurements[i];
      x += (fMeasurements[i] * fPositions[i]);
      x2 += (fMeasurements[i] * fPositions[i] * fPositions[i]);
      c += dx * fMeasurements[i];
   }
   double mean = x / norm;
   double rms2 = x2 / norm - mean * mean;

   //   std::cout<<"FCN initial mean: "<<mean<<std::endl;
   //   std::cout<<"FCN initial sigma: "<<std::sqrt(rms2)<<std::endl;

   std::vector<double> par;
   par.push_back(mean);
   par.push_back(std::sqrt(rms2));
   par.push_back(c);
   par.push_back(mean);
   par.push_back(std::sqrt(rms2));
   par.push_back(c);

   fMin = (*this)(par);
   //   std::cout<<"GaussFcnHistoData2 initial chi2: "<<fMin<<std::endl;
}

} // namespace Minuit2

} // namespace ROOT
