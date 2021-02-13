// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#include "GaussFcn.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnPrint.h"

#include <iostream>
#include <fstream>

#ifdef USE_SEALBASE
#include "SealBase/Filename.h"
#include "SealBase/ShellEnvironment.h"
#endif

using namespace ROOT::Minuit2;

int main()
{

   std::vector<double> positions;
   std::vector<double> measurements;
   std::vector<double> var;
   double nmeas = 0;

#ifdef USE_SEALBASE
   seal::Filename inputFile(
      seal::Filename("$SEAL/src/MathLibs/Minuit/tests/MnSim/paul2.txt").substitute(seal::ShellEnvironment()));
   std::ifstream in(inputFile.Name());
#else
   std::ifstream in("paul2.txt");
#endif
   if (!in) {
      std::cerr << "Error opening input data file" << std::endl;
      return 1;
   }

   // read input data
   {
      double x = 0., y = 0., width = 0., err = 0., un1 = 0., un2 = 0.;
      while (in >> x >> y >> width >> err >> un1 >> un2) {
         if (err < 1.e-8)
            continue;
         positions.push_back(x);
         measurements.push_back(y);
         var.push_back(err * err);
         nmeas += y;
      }
      std::cout << "size= " << var.size() << std::endl;
      assert(var.size() > 0);
      std::cout << "nmeas: " << nmeas << std::endl;
   }

   // create FCN function
   GaussFcn fFCN(measurements, positions, var);

   std::vector<double> meas = fFCN.Measurements();
   std::vector<double> pos = fFCN.Positions();

   // create initial starting values for parameters
   double x = 0.;
   double x2 = 0.;
   double norm = 0.;
   double area = 0.;
   double dx = pos[1] - pos[0];
   for (unsigned int i = 0; i < meas.size(); i++) {
      norm += meas[i];
      x += (meas[i] * pos[i]);
      x2 += (meas[i] * pos[i] * pos[i]);
      area += dx * meas[i];
   }
   double mean = x / norm;
   double rms2 = x2 / norm - mean * mean;

   std::cout << "initial mean: " << mean << std::endl;
   std::cout << "initial sigma: " << std::sqrt(rms2) << std::endl;
   std::cout << "initial area: " << area << std::endl;
   std::vector<double> init_val(3);
   init_val[0] = mean;
   init_val[1] = std::sqrt(rms2);
   init_val[2] = area;
   std::cout << "initial fval: " << fFCN(init_val) << std::endl;

   MnUserParameters upar;
   upar.Add("mean", mean, 1.);
   upar.Add("sigma", std::sqrt(rms2), 1.);
   upar.Add("area", area, 10.);

   MnMigrad migrad(fFCN, upar);
   std::cout << "start migrad " << std::endl;
   FunctionMinimum min = migrad();
   std::cout << "minimum: " << min << std::endl;

   std::cout << "start Minos" << std::endl;
   MnMinos Minos(fFCN, min);
   std::pair<double, double> e0 = Minos(0);
   std::pair<double, double> e1 = Minos(1);
   std::pair<double, double> e2 = Minos(2);

   std::cout << "par0: " << min.UserState().Value("mean") << " " << e0.first << " " << e0.second << std::endl;
   std::cout << "par1: " << min.UserState().Value(1) << " " << e1.first << " " << e1.second << std::endl;
   std::cout << "par2: " << min.UserState().Value(2) << " " << e2.first << " " << e2.second << std::endl;

   return 0;
}
