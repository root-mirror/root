// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#include "Quad12F.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnPrint.h"
#include <iostream>

// #include "TimingUtilities/PentiumTimer.h"

// StackAllocator gStackAllocator;

using namespace ROOT::Minuit2;

int main()
{

   Quad12F fcn;

   //   PentiumTimer stopwatch;
   //   stopwatch.start();

   //   long long int start = stopwatch.lap().ticks();
   //   long long int stop = stopwatch.lap().ticks();
   //   std::cout<<"stop-start: "<<stop - start<<std::endl;
   //   start = stopwatch.lap().ticks();
   {
      // test constructor
      MnUserParameters upar;
      upar.Add("x", 1., 0.1);
      upar.Add("y", 1., 0.1);
      upar.Add("z", 1., 0.1);
      upar.Add("w", 1., 0.1);
      upar.Add("x0", 1., 0.1);
      upar.Add("y0", 1., 0.1);
      upar.Add("z0", 1., 0.1);
      upar.Add("w0", 1., 0.1);
      upar.Add("x1", 1., 0.1);
      upar.Add("y1", 1., 0.1);
      upar.Add("z1", 1., 0.1);
      upar.Add("w1", 1., 0.1);

      MnMigrad migrad(fcn, upar);
      FunctionMinimum min = migrad();
      std::cout << "minimum: " << min << std::endl;
   }
   //   stop = stopwatch.lap().ticks();
   //   std::cout<<"stop-start: "<<stop - start<<std::endl;
   /*

     {
     //test constructor
     std::vector<double> par(8);
     std::vector<double> err(8);
     for(int i = 0; i < 8; i++) {
       par[i] = 1.;
       err[i] = 0.1;
     }
     MnMigrad migrad(fcn, par, err);
     FunctionMinimum min = migrad();
     std::cout<<"minimum: "<<min<<std::endl;
     }

     {
     //test edm Value
     std::vector<double> par(8);
     std::vector<double> err(8);
     for(int i = 0; i < 8; i++) {
       par[i] = 1.;
       err[i] = 0.1;
     }
     MnMigrad migrad(fcn, par, err);
     double edm = 1.e-1;
     FunctionMinimum min = migrad(20, edm);
     std::cout<<"minimum: "<<min<<std::endl;
     }

     {
     //test # of iterations
     std::vector<double> par(4);
     std::vector<double> err(4);
     for(int i = 0; i < 4; i++) {
       par[i] = 1.;
       err[i] = 0.1;
     }
     MnMigrad migrad(fcn, par, err);
     int niter = 2;
     FunctionMinimum min = migrad(niter, 1.e-5);
     std::cout<<"minimum: "<<min<<std::endl;
     }
   */

   return 0;
}
