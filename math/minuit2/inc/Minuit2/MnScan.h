// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Minuit2_MnScan
#define ROOT_Minuit2_MnScan

#include "Minuit2/MnApplication.h"
#include "Minuit2/ScanMinimizer.h"

#include <vector>
#include <utility>

namespace ROOT {

namespace Minuit2 {

class FCNBase;

//_______________________________________________________________________
/**
    API class for minimization using a scan method to find the minimum;
    allows for user interaction: set/change parameters, do minimization,
    change parameters, re-do minimization etc.;

 */

class MnScan : public MnApplication {

public:
   /// construct from FCNBase + std::vector for parameters and errors
   MnScan(const FCNBase &fcn, const std::vector<double> &par, const std::vector<double> &err, unsigned int stra = 1)
      : MnApplication(fcn, MnUserParameterState(par, err), MnStrategy(stra)), fMinimizer(ScanMinimizer())
   {
   }

   /// construct from FCNBase + std::vector for parameters and covariance
   MnScan(const FCNBase &fcn, const std::vector<double> &par, unsigned int nrow, const std::vector<double> &cov,
          unsigned int stra = 1)
      : MnApplication(fcn, MnUserParameterState(par, cov, nrow), MnStrategy(stra)), fMinimizer(ScanMinimizer())
   {
   }

   /// construct from FCNBase + std::vector for parameters and MnUserCovariance
   MnScan(const FCNBase &fcn, const std::vector<double> &par, const MnUserCovariance &cov, unsigned int stra = 1)
      : MnApplication(fcn, MnUserParameterState(par, cov), MnStrategy(stra)), fMinimizer(ScanMinimizer())
   {
   }

   /// construct from FCNBase + MnUserParameters
   MnScan(const FCNBase &fcn, const MnUserParameters &par, unsigned int stra = 1)
      : MnApplication(fcn, MnUserParameterState(par), MnStrategy(stra)), fMinimizer(ScanMinimizer())
   {
   }

   /// construct from FCNBase + MnUserParameters + MnUserCovariance
   MnScan(const FCNBase &fcn, const MnUserParameters &par, const MnUserCovariance &cov, unsigned int stra = 1)
      : MnApplication(fcn, MnUserParameterState(par, cov), MnStrategy(stra)), fMinimizer(ScanMinimizer())
   {
   }

   /// construct from FCNBase + MnUserParameterState + MnStrategy
   MnScan(const FCNBase &fcn, const MnUserParameterState &par, const MnStrategy &str)
      : MnApplication(fcn, MnUserParameterState(par), str), fMinimizer(ScanMinimizer())
   {
   }

   MnScan(const MnScan &migr)
      : MnApplication(migr.Fcnbase(), migr.State(), migr.Strategy(), migr.NumOfCalls()), fMinimizer(migr.fMinimizer)
   {
   }

   ~MnScan() {}

   ModularFunctionMinimizer &Minimizer() { return fMinimizer; }
   const ModularFunctionMinimizer &Minimizer() const { return fMinimizer; }

   std::vector<std::pair<double, double>>
   Scan(unsigned int par, unsigned int maxsteps = 41, double low = 0., double high = 0.);

private:
   ScanMinimizer fMinimizer;

private:
   /// forbidden assignment (const FCNBase& = )
   MnScan &operator=(const MnScan &) { return *this; }
};

} // namespace Minuit2

} // namespace ROOT

#endif // ROOT_Minuit2_MnScan
