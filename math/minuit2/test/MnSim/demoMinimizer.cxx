// @(#)root/minuit2:$Id$
// Author: L. Moneta    12/2005
/**
   test of a pure minimization passing a user function
   This is an example of running Minuit2 using the
   Minuit2Minimizer class (via the Minimizer interface)

*/
#include "Math/IFunction.h"
#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"
#include <string>
#include <iostream>

double Rosenbrock(const double *x)
{
   // Rosebrock function
   double tmp1 = (x[1] - x[0] * x[0]); // y -x^2
   double tmp2 = (1. - x[0]);          // 1-x
   return 100. * tmp1 * tmp1 + tmp2 * tmp2;
}

int demoMinimizer(const char *algoName, int printlevel)
{

   ROOT::Math::Minimizer *min = new ROOT::Minuit2::Minuit2Minimizer(algoName);

   // set tolerance , etc...
   min->SetMaxFunctionCalls(1000000);
   min->SetTolerance(0.001);
   min->SetPrintLevel(printlevel);

   // create funciton wrapper for minmizer
   // a IMultiGenFunction type
   ROOT::Math::Functor f(&Rosenbrock, 2);

   // starting point
   double variable[2] = {-1., 1.2};
   // initial spet sizes
   double step[2] = {0.01, 0.01};

   min->SetFunction(f);

   // Set the free variables to be minimized!
   min->SetVariable(0, "x", variable[0], step[0]);
   min->SetVariable(1, "y", variable[1], step[1]);

   // do the minimization
   min->Minimize();

   const double *xs = min->X();
   std::cout << "Minimum: f(" << xs[0] << "," << xs[1] << "): " << min->MinValue() << std::endl;

   // expected minimum is 0
   if (min->MinValue() < 1.E-4 && f(xs) < 1.E-4) {
      std::cout << "Minuit2 -  " << algoName << "   converged to the right minimum" << std::endl;
      return 0;
   } else
      std::cerr << "ERROR:  Minuit2 - " << algoName << "   failed to converge !!!" << std::endl;

   return -1;
}

int main(int argc, const char *argv[])
{

   int printLevel = 0;
   std::string algoName = ""; // use default (i.e. migrad)

   // Parse command line arguments
   for (int i = 1; i < argc; i++) {
      std::string arg = argv[i];

      if (arg == "-v") {
         printLevel = 1;
      } else if (arg == "-vv") {
         printLevel = 2;
      } else if (arg == "-vvv") {
         printLevel = 3;
      } else if (arg == "-n") {
         std::cout << "using method " << argv[i + 1] << " to minimize" << std::endl;
         algoName = argv[++i];
      } else if (arg == "-h") {
         std::cout << "usage: demoMinimize [ options ] " << std::endl;
         std::cout << "" << std::endl;
         std::cout << "       -n <algorithm> : use given algorithm (possible names: simplex, minimize, scan, fumili)"
                   << std::endl;
         std::cout << "       -v    : set minimul verbose mode to show final result" << std::endl;
         std::cout << "       -vv   : set medium verbose mode: show function value and edm at each minimization step"
                   << std::endl;
         std::cout << "       -vvv  : set very verbose mode: show full result at each minimization step" << std::endl;
         return 0;
      }
   }

   int iret = demoMinimizer(algoName.c_str(), printLevel);
   return iret;
}
