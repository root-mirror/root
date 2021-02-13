// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#include "GaussDataGen.h"
#include "GaussianModelFunction.h"
#include "Minuit2/MnFumiliMinimize.h"
#include "Minuit2/FumiliStandardChi2FCN.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnMigrad.h"

#include <iostream>

using namespace ROOT::Minuit2;

int main()
{

   // generate the data (100 data points)
   GaussDataGen gdg(100);

   std::vector<double> pos = gdg.Positions();
   std::vector<double> meas = gdg.Measurements();
   std::vector<double> var = gdg.Variances();

   // Estimate initial starting values for parameters
   double x = 0.;
   double x2 = 0.;
   double norm = 0.;
   double dx = pos[1] - pos[0];
   double area = 0.;
   for (unsigned int i = 0; i < meas.size(); i++) {
      norm += meas[i];
      x += (meas[i] * pos[i]);
      x2 += (meas[i] * pos[i] * pos[i]);
      area += dx * meas[i];
   }
   double mean = x / norm;
   double rms2 = x2 / norm - mean * mean;
   double rms = rms2 > 0. ? std::sqrt(rms2) : 1.;

   // create parameters
   MnUserParameters upar;
   upar.Add("mean", mean, 0.1);
   upar.Add("sigma", rms, 0.1);
   upar.Add("area", area, 0.1);

   // create FCN function for Fumili using model function
   GaussianModelFunction modelFunction;
   FumiliStandardChi2FCN fFCN(modelFunction, meas, pos, var);

   {

      std::cout << "Minimize using FUMILI : \n" << std::endl;
      MnFumiliMinimize fumili(fFCN, upar);

      // Minimize
      FunctionMinimum min = fumili();

      // output
      std::cout << "minimum: " << min << std::endl;
   }

   {

      std::cout << "Minimize using MIGRAD : \n" << std::endl;
      MnMigrad migrad(fFCN, upar);

      // Minimize
      FunctionMinimum min = migrad();

      // output
      std::cout << "minimum: " << min << std::endl;
   }

   return 0;
}
