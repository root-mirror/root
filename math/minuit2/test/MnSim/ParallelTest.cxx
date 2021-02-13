// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#include "GaussRandomGen.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnPlot.h"
#include "Minuit2/MinosError.h"
#include "Minuit2/FCNBase.h"
#include <cmath>
#include <iostream>

// example of a multi dimensional fit where parallelization can be used
// to speed up the result
// define the environment variable OMP_NUM_THREADS to the number of desired threads
// By default it will have thenumber of core of the machine
// The default number of dimension is 20 (fit in 40 parameters) on 1000 data events.
// One can change the dimension and the number of events by doing:
// ./test_Minuit2_Parallel    ndim  nevents

using namespace ROOT::Minuit2;

const int default_ndim = 20;
const int default_ndata = 1000;

double GaussPdf(double x, double x0, double sigma)
{
   double tmp = (x - x0) / sigma;
   constexpr double two_pi = 2 * 3.14159265358979323846; // M_PI is not standard
   return (1.0 / (std::sqrt(two_pi) * std::fabs(sigma))) * std::exp(-tmp * tmp / 2);
}

double LogMultiGaussPdf(const std::vector<double> &x, const std::vector<double> &p)
{
   double f = 0;
   int ndim = x.size();
   for (int k = 0; k < ndim; ++k) {
      double y = GaussPdf(x[k], p[2 * k], p[2 * k + 1]);
      // std::cout << " k " << k << "  " << y << "  " << p[2*k] << "  " << p[2*k+1] << std::endl;
      y = std::max(y, 1.E-300);
      double lgaus = std::log(y);
      f += lgaus;
   }
   return f;
}

typedef std::vector<std::vector<double>> Data;

struct LogLikeFCN : public FCNBase {

   LogLikeFCN(const Data &data) : fData(data) {}

   double operator()(const std::vector<double> &p) const
   {
      double logl = 0;
      int ndata = fData.size();
      for (int i = 0; i < ndata; ++i) {
         logl -= LogMultiGaussPdf(fData[i], p);
      }
      return logl;
   }
   double Up() const { return 0.5; }
   const Data &fData;
};

int doFit(int ndim, int ndata)
{

   // generate the data (1000 data points) in 100 dimension

   Data data(ndata);
   std::vector<double> event(ndim);

   std::vector<double> mean(ndim);
   std::vector<double> sigma(ndim);
   for (int k = 0; k < ndim; ++k) {
      mean[k] = -double(ndim) / 2 + k;
      sigma[k] = 1. + 0.1 * k;
   }

   // generate the data
   for (int i = 0; i < ndata; ++i) {
      for (int k = 0; k < ndim; ++k) {
         GaussRandomGen rgaus(mean[k], sigma[k]);
         event[k] = rgaus();
      }
      data[i] = event;
   }

   // create FCN function
   LogLikeFCN fcn(data);

   // create initial starting values for parameters
   std::vector<double> init_par(2 * ndim);
   for (int k = 0; k < ndim; ++k) {
      init_par[2 * k] = 0;
      init_par[2 * k + 1] = 1;
   }

   std::vector<double> init_err(2 * ndim);
   for (int k = 0; k < 2 * ndim; ++k) {
      init_err[k] = 0.1;
   }
   // create minimizer (default constructor)
   VariableMetricMinimizer fMinimizer;

   // Minimize
   FunctionMinimum min = fMinimizer.Minimize(fcn, init_par, init_err);

   // output
   std::cout << "minimum: " << min << std::endl;

   //     // create MINOS Error factory
   //     MnMinos Minos(fFCN, min);

   //     {
   //       // 3-sigma MINOS errors (minimal interface)
   //       fFCN.SetErrorDef(9.);
   //       std::pair<double,double> e0 = Minos(0);
   //       std::pair<double,double> e1 = Minos(1);
   //       std::pair<double,double> e2 = Minos(2);

   //       // output
   //       std::cout<<"3-sigma Minos errors with limits: "<<std::endl;
   //       std::cout.precision(16);
   //       std::cout<<"par0: "<<min.UserState().Value("mean")<<" "<<e0.first<<" "<<e0.second<<std::endl;
   //       std::cout<<"par1: "<<min.UserState().Value(1)<<" "<<e1.first<<" "<<e1.second<<std::endl;
   //       std::cout<<"par2: "<<min.UserState().Value("area")<<" "<<e2.first<<" "<<e2.second<<std::endl;

   //     }

   //   }

   return 0;
}

int main(int argc, char **argv)
{
   int ndim = default_ndim;
   int ndata = default_ndata;
   if (argc > 1) {
      ndim = atoi(argv[1]);
   }
   if (argc > 2) {
      ndata = atoi(argv[2]);
   }
   std::cout << "do fit of " << ndim << " dimensional data on " << ndata << " events " << std::endl;
   doFit(ndim, ndata);
}
