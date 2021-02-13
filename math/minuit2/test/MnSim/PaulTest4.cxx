// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#include "Minuit2/FCNBase.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/SimplexMinimizer.h"

#include <iostream>
#include <fstream>
#include <cassert>

#ifdef USE_SEALBASE
#include "SealBase/Filename.h"
#include "SealBase/ShellEnvironment.h"
#endif

using namespace ROOT::Minuit2;

class PowerLawFunc {

public:
   PowerLawFunc(double _p0, double _p1) : fP0(_p0), fP1(_p1) {}

   ~PowerLawFunc() {}

   double operator()(double x) const { return p1() * std::exp(std::log(x) * p0()); }

   double p0() const { return fP0; }
   double p1() const { return fP1; }

private:
   double fP0;
   double fP1;
};

class PowerLawChi2FCN : public FCNBase {

public:
   PowerLawChi2FCN(const std::vector<double> &meas, const std::vector<double> &pos, const std::vector<double> &mvar)
      : fMeasurements(meas), fPositions(pos), fMVariances(mvar)
   {
   }

   ~PowerLawChi2FCN() {}

   double operator()(const std::vector<double> &par) const
   {
      assert(par.size() == 2);
      PowerLawFunc pl(par[0], par[1]);
      double chi2 = 0.;

      for (unsigned int n = 0; n < fMeasurements.size(); n++) {
         chi2 += ((pl(fPositions[n]) - fMeasurements[n]) * (pl(fPositions[n]) - fMeasurements[n]) / fMVariances[n]);
      }

      return chi2;
   }

   double Up() const { return 1.; }

private:
   std::vector<double> fMeasurements;
   std::vector<double> fPositions;
   std::vector<double> fMVariances;
};

class PowerLawLogLikeFCN : public FCNBase {

public:
   PowerLawLogLikeFCN(const std::vector<double> &meas, const std::vector<double> &pos)
      : fMeasurements(meas), fPositions(pos)
   {
   }

   ~PowerLawLogLikeFCN() {}

   double operator()(const std::vector<double> &par) const
   {
      assert(par.size() == 2);
      PowerLawFunc pl(par[0], par[1]);
      double logsum = 0.;

      for (unsigned int n = 0; n < fMeasurements.size(); n++) {
         double k = fMeasurements[n];
         double mu = pl(fPositions[n]);
         logsum += (k * std::log(mu) - mu);
      }

      return -logsum;
   }

   double Up() const { return 0.5; }

private:
   std::vector<double> fMeasurements;
   std::vector<double> fPositions;
};

int main()
{

   std::vector<double> positions;
   std::vector<double> measurements;
   std::vector<double> var;
   {

#ifdef USE_SEALBASE
      seal::Filename inputFile(
         seal::Filename("$SEAL/src/MathLibs/Minuit/tests/MnSim/paul4.txt").substitute(seal::ShellEnvironment()));
      std::ifstream in(inputFile.Name());
#else
      std::ifstream in("paul4.txt");
#endif
      if (!in) {
         std::cerr << "Error opening input data file" << std::endl;
         return 1;
      }

      double x = 0., y = 0., err = 0.;
      while (in >> x >> y >> err) {
         //       if(err < 1.e-8) continue;
         positions.push_back(x);
         measurements.push_back(y);
         var.push_back(err * err);
      }
      std::cout << "size= " << var.size() << std::endl;
   }
   {
      // create Chi2 FCN function
      std::cout << ">>> test Chi2" << std::endl;
      PowerLawChi2FCN fFCN(measurements, positions, var);

      MnUserParameters upar;
      upar.Add("p0", -2.3, 0.2);
      upar.Add("p1", 1100., 10.);

      MnMigrad migrad(fFCN, upar);
      std::cout << "start migrad " << std::endl;
      FunctionMinimum min = migrad();
      if (!min.IsValid()) {
         // try with higher strategy
         std::cout << "FM is invalid, try with strategy = 2." << std::endl;
         MnMigrad migrad2(fFCN, upar, 2);
         min = migrad2();
      }
      std::cout << "minimum: " << min << std::endl;
   }
   {
      std::cout << ">>> test log LikeliHood" << std::endl;
      // create LogLikelihood FCN function
      PowerLawLogLikeFCN fFCN(measurements, positions);

      MnUserParameters upar;
      upar.Add("p0", -2.1, 0.2);
      upar.Add("p1", 1000., 10.);

      MnMigrad migrad(fFCN, upar);
      std::cout << "start migrad " << std::endl;
      FunctionMinimum min = migrad();
      if (!min.IsValid()) {
         // try with higher strategy
         std::cout << "FM is invalid, try with strategy = 2." << std::endl;
         MnMigrad migrad2(fFCN, upar, 2);
         min = migrad2();
      }
      std::cout << "minimum: " << min << std::endl;
   }
   {
      std::cout << ">>> test Simplex" << std::endl;
      PowerLawChi2FCN chi2(measurements, positions, var);
      PowerLawLogLikeFCN mlh(measurements, positions);

      MnUserParameters upar;
      std::vector<double> par;
      par.push_back(-2.3);
      par.push_back(1100.);
      std::vector<double> err;
      err.push_back(1.);
      err.push_back(1.);

      SimplexMinimizer simplex;

      std::cout << "start simplex" << std::endl;
      FunctionMinimum min = simplex.Minimize(chi2, par, err);
      std::cout << "minimum: " << min << std::endl;
      FunctionMinimum min2 = simplex.Minimize(mlh, par, err);
      std::cout << "minimum: " << min2 << std::endl;
   }
   return 0;
}
