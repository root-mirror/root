// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#include "GaussFcn.h"
#include "GaussDataGen.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnContours.h"
#include "Minuit2/MnPlot.h"
#include "Minuit2/MinosError.h"
#include "Minuit2/ContoursError.h"

#include <iostream>

using namespace ROOT::Minuit2;

int main()
{

   // generate the data (100 data points)
   GaussDataGen gdg(100);

   std::vector<double> pos = gdg.Positions();
   std::vector<double> meas = gdg.Measurements();
   std::vector<double> var = gdg.Variances();

   // create FCN function
   GaussFcn fFCN(meas, pos, var);

   // create initial starting values for parameters
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

   {
      // demonstrate minimal required interface for minimization
      // create Minuit parameters without names

      // starting values for parameters
      std::vector<double> init_par;
      init_par.push_back(mean);
      init_par.push_back(rms);
      init_par.push_back(area);

      // starting values for initial uncertainties
      std::vector<double> init_err;
      init_err.push_back(0.1);
      init_err.push_back(0.1);
      init_err.push_back(0.1);

      // create minimizer (default constructor)
      VariableMetricMinimizer fMinimizer;

      // Minimize
      FunctionMinimum min = fMinimizer.Minimize(fFCN, init_par, init_err);

      // output
      std::cout << "minimum: " << min << std::endl;
   }

   {
      // demonstrate standard minimization using MIGRAD
      // create Minuit parameters with names
      MnUserParameters upar;
      upar.Add("mean", mean, 0.1);
      upar.Add("sigma", rms, 0.1);
      upar.Add("area", area, 0.1);

      // create MIGRAD minimizer
      MnMigrad migrad(fFCN, upar);

      // Minimize
      FunctionMinimum min = migrad();

      // output
      std::cout << "minimum: " << min << std::endl;
   }

   {
      // demonstrate full interaction with parameters over subsequent
      // minimizations

      // create Minuit parameters with names
      MnUserParameters upar;
      upar.Add("mean", mean, 0.1);
      upar.Add("sigma", rms, 0.1);
      upar.Add("area", area, 0.1);

      // access Parameter by Name to set limits...
      upar.SetLimits("mean", mean - 0.01, mean + 0.01);

      // ... or access Parameter by Index
      upar.SetLimits(1, rms - 0.1, rms + 0.1);

      // create Migrad minimizer
      MnMigrad migrad(fFCN, upar);

      // Fix a Parameter...
      migrad.Fix("mean");

      // ... and Minimize
      FunctionMinimum min = migrad();

      // output
      std::cout << "minimum: " << min << std::endl;

      // Release a Parameter...
      migrad.Release("mean");

      // ... and Fix another one
      migrad.Fix(1);

      // and Minimize again
      FunctionMinimum min1 = migrad();

      // output
      std::cout << "minimum1: " << min1 << std::endl;

      // Release the Parameter...
      migrad.Release(1);

      // ... and Minimize with all three parameters (still with limits!)
      FunctionMinimum min2 = migrad();

      // output
      std::cout << "minimum2: " << min2 << std::endl;

      // remove all limits on parameters...
      migrad.RemoveLimits("mean");
      migrad.RemoveLimits("sigma");

      // ... and Minimize again with all three parameters (now without limits!)
      FunctionMinimum min3 = migrad();

      // output
      std::cout << "minimum3: " << min3 << std::endl;
   }

   {
      // test single sided limits
      MnUserParameters upar;
      upar.Add("mean", mean, 0.1);
      upar.Add("sigma", rms - 1., 0.1);
      upar.Add("area", area, 0.1);

      // test Lower limits
      upar.SetLowerLimit("mean", mean - 0.01);

      // test Upper limits
      upar.SetUpperLimit("sigma", rms - 0.5);

      // create MIGRAD minimizer
      MnMigrad migrad(fFCN, upar);

      // ... and Minimize
      FunctionMinimum min = migrad();
      std::cout << "test Lower limit minimim= " << min << std::endl;
   }

   {
      // demonstrate MINOS Error analysis

      // create Minuit parameters with names
      MnUserParameters upar;
      upar.Add("mean", mean, 0.1);
      upar.Add("sigma", rms, 0.1);
      upar.Add("area", area, 0.1);

      // create Migrad minimizer
      MnMigrad migrad(fFCN, upar);

      // Minimize
      FunctionMinimum min = migrad();

      // create MINOS Error factory
      MnMinos Minos(fFCN, min);

      {
         // 1-sigma MINOS errors (minimal interface)
         std::pair<double, double> e0 = Minos(0);
         std::pair<double, double> e1 = Minos(1);
         std::pair<double, double> e2 = Minos(2);

         // output
         std::cout << "1-sigma Minos errors: " << std::endl;
         std::cout << "par0: " << min.UserState().Value("mean") << " " << e0.first << " " << e0.second << std::endl;
         std::cout << "par1: " << min.UserState().Value(1) << " " << e1.first << " " << e1.second << std::endl;
         std::cout << "par2: " << min.UserState().Value("area") << " " << e2.first << " " << e2.second << std::endl;
      }

      {
         // 2-sigma MINOS errors (rich interface)
         fFCN.SetErrorDef(4.);
         MinosError e0 = Minos.Minos(0);
         MinosError e1 = Minos.Minos(1);
         MinosError e2 = Minos.Minos(2);

         // output
         std::cout << "2-sigma Minos errors: " << std::endl;
         std::cout << e0 << std::endl;
         std::cout << e1 << std::endl;
         std::cout << e2 << std::endl;
      }
   }

   {
      // demostrate MINOS Error analysis with limits

      // create Minuit parameters with names
      MnUserParameters upar;
      upar.Add("mean", mean, 0.1);
      upar.Add("sigma", rms, 0.1);
      upar.Add("area", area, 0.1);

      double meanLow = -50.03;
      double rmsUp = 1.55;
      std::cout << "sigma Limit: " << rmsUp << "\tmean limit: " << meanLow << std::endl;
      // test Lower limits
      upar.SetLowerLimit("mean", meanLow);
      // test Upper limits
      upar.SetUpperLimit("sigma", rmsUp);

      // create Migrad minimizer
      MnMigrad migrad(fFCN, upar);

      // Minimize
      FunctionMinimum min = migrad();

      // create MINOS Error factory
      MnMinos Minos(fFCN, min);

      {
         // 3-sigma MINOS errors (minimal interface)
         fFCN.SetErrorDef(9.);
         std::pair<double, double> e0 = Minos(0);
         std::pair<double, double> e1 = Minos(1);
         std::pair<double, double> e2 = Minos(2);

         // output
         std::cout << "3-sigma Minos errors with limits: " << std::endl;
         std::cout.precision(16);
         std::cout << "par0: " << min.UserState().Value("mean") << " " << e0.first << " " << e0.second << std::endl;
         std::cout << "par1: " << min.UserState().Value(1) << " " << e1.first << " " << e1.second << std::endl;
         std::cout << "par2: " << min.UserState().Value("area") << " " << e2.first << " " << e2.second << std::endl;
      }
   }

   {
      // demonstrate how to use the CONTOURs

      // create Minuit parameters with names
      MnUserParameters upar;
      upar.Add("mean", mean, 0.1);
      upar.Add("sigma", rms, 0.1);
      upar.Add("area", area, 0.1);

      // create Migrad minimizer
      MnMigrad migrad(fFCN, upar);

      // Minimize
      FunctionMinimum min = migrad();

      // create contours factory with FCN and Minimum
      MnContours contours(fFCN, min);

      // 70% confidence level for 2 parameters Contour around the Minimum
      // (minimal interface)
      fFCN.SetErrorDef(2.41);
      std::vector<std::pair<double, double>> cont = contours(0, 1, 20);

      // 95% confidence level for 2 parameters Contour
      // (rich interface)
      fFCN.SetErrorDef(5.99);
      ContoursError cont4 = contours.Contour(0, 1, 20);

      // plot the contours
      MnPlot plot;
      cont.insert(cont.end(), cont4().begin(), cont4().end());
      plot(min.UserState().Value("mean"), min.UserState().Value("sigma"), cont);

      // print out one Contour
      std::cout << cont4 << std::endl;
   }

   return 0;
}
