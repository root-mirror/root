// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#include "Quad1F.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/VariableMetricMinimizer.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnMinos.h"
#include <iostream>

using namespace ROOT::Minuit2;

int main()
{

   {
      // test constructor
      {
         Quad1F fcn;
         MnUserParameters upar;
         upar.Add("x", 1., 0.1);
         MnMigrad migrad(fcn, upar);
         FunctionMinimum min = migrad();
         std::cout << "min= " << min << std::endl;
      }
      {
         // using VariableMetricMinimizer, analytical derivatives
         Quad1F fcn;
         std::vector<double> par(1, 1.);
         std::vector<double> err(1, 0.1);
         VariableMetricMinimizer mini;
         FunctionMinimum min = mini.Minimize(fcn, par, err);
         std::cout << "min= " << min << std::endl;
      }
      {
         // test Minos for one Parameter
         Quad1F fcn;
         std::vector<double> par(1, 1.);
         std::vector<double> err(1, 0.1);
         VariableMetricMinimizer mini;
         FunctionMinimum min = mini.Minimize(fcn, par, err);
         MnMinos Minos(fcn, min);
         std::pair<double, double> e0 = Minos(0);
         std::cout << "par0: " << min.UserState().Value(unsigned(0)) << " " << e0.first << " " << e0.second
                   << std::endl;
         fcn.SetErrorDef(4.);
         MnMinos minos2(fcn, min);
         std::pair<double, double> e02 = minos2(0);
         std::cout << "par0: " << min.UserState().Value(unsigned(0)) << " " << e02.first << " " << e02.second
                   << std::endl;
      }
   }

   return 0;
}
