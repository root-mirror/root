// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Minuit2_MnMinimize
#define ROOT_Minuit2_MnMinimize

#include "Minuit2/MnApplication.h"
#include "Minuit2/CombinedMinimizer.h"

#include <vector>

namespace ROOT {

namespace Minuit2 {

class FCNBase;

/** API class for minimization using Variable Metric technology ("MIGRAD");
    allows for user interaction: set/change parameters, do minimization,
    change parameters, re-do minimization etc.;
    also used by MnMinos and MnContours;
 */

class MnMinimize : public MnApplication {

public:
   /// construct from FCNBase + std::vector for parameters and errors
   MnMinimize(const FCNBase &fcn, const std::vector<double> &par, const std::vector<double> &err, unsigned int stra = 1)
      : MnApplication(fcn, MnUserParameterState(par, err), MnStrategy(stra)), fMinimizer(CombinedMinimizer())
   {
   }

   /// construct from FCNBase + std::vector for parameters and covariance
   MnMinimize(const FCNBase &fcn, const std::vector<double> &par, unsigned int nrow, const std::vector<double> &cov,
              unsigned int stra = 1)
      : MnApplication(fcn, MnUserParameterState(par, cov, nrow), MnStrategy(stra)), fMinimizer(CombinedMinimizer())
   {
   }

   /// construct from FCNBase + std::vector for parameters and MnUserCovariance
   MnMinimize(const FCNBase &fcn, const std::vector<double> &par, const MnUserCovariance &cov, unsigned int stra = 1)
      : MnApplication(fcn, MnUserParameterState(par, cov), MnStrategy(stra)), fMinimizer(CombinedMinimizer())
   {
   }

   /// construct from FCNBase + MnUserParameters
   MnMinimize(const FCNBase &fcn, const MnUserParameters &par, unsigned int stra = 1)
      : MnApplication(fcn, MnUserParameterState(par), MnStrategy(stra)), fMinimizer(CombinedMinimizer())
   {
   }

   /// construct from FCNBase + MnUserParameters + MnUserCovariance
   MnMinimize(const FCNBase &fcn, const MnUserParameters &par, const MnUserCovariance &cov, unsigned int stra = 1)
      : MnApplication(fcn, MnUserParameterState(par, cov), MnStrategy(stra)), fMinimizer(CombinedMinimizer())
   {
   }

   /// construct from FCNBase + MnUserParameterState + MnStrategy
   MnMinimize(const FCNBase &fcn, const MnUserParameterState &par, const MnStrategy &str)
      : MnApplication(fcn, MnUserParameterState(par), str), fMinimizer(CombinedMinimizer())
   {
   }

   // interfaces using FCNGradientBase

   /// construct from FCNGradientBase + std::vector for parameters and errors
   MnMinimize(const FCNGradientBase &fcn, const std::vector<double> &par, const std::vector<double> &err,
              unsigned int stra = 1)
      : MnApplication(fcn, MnUserParameterState(par, err), MnStrategy(stra)), fMinimizer(CombinedMinimizer())
   {
   }

   /// construct from FCNGradientBase + std::vector for parameters and covariance
   MnMinimize(const FCNGradientBase &fcn, const std::vector<double> &par, unsigned int nrow,
              const std::vector<double> &cov, unsigned int stra = 1)
      : MnApplication(fcn, MnUserParameterState(par, cov, nrow), MnStrategy(stra)), fMinimizer(CombinedMinimizer())
   {
   }

   /// construct from FCNGradientBase + std::vector for parameters and MnUserCovariance
   MnMinimize(const FCNGradientBase &fcn, const std::vector<double> &par, const MnUserCovariance &cov,
              unsigned int stra = 1)
      : MnApplication(fcn, MnUserParameterState(par, cov), MnStrategy(stra)), fMinimizer(CombinedMinimizer())
   {
   }

   /// construct from FCNGradientBase + MnUserParameters
   MnMinimize(const FCNGradientBase &fcn, const MnUserParameters &par, unsigned int stra = 1)
      : MnApplication(fcn, MnUserParameterState(par), MnStrategy(stra)), fMinimizer(CombinedMinimizer())
   {
   }

   /// construct from FCNGradientBase + MnUserParameters + MnUserCovariance
   MnMinimize(const FCNGradientBase &fcn, const MnUserParameters &par, const MnUserCovariance &cov,
              unsigned int stra = 1)
      : MnApplication(fcn, MnUserParameterState(par, cov), MnStrategy(stra)), fMinimizer(CombinedMinimizer())
   {
   }

   /// construct from FCNGradientBase + MnUserParameterState + MnStrategy
   MnMinimize(const FCNGradientBase &fcn, const MnUserParameterState &par, const MnStrategy &str)
      : MnApplication(fcn, MnUserParameterState(par), str), fMinimizer(CombinedMinimizer())
   {
   }

   MnMinimize(const MnMinimize &migr)
      : MnApplication(migr.Fcnbase(), migr.State(), migr.Strategy(), migr.NumOfCalls()), fMinimizer(migr.fMinimizer)
   {
   }

   ~MnMinimize() {}

   ModularFunctionMinimizer &Minimizer() { return fMinimizer; }
   const ModularFunctionMinimizer &Minimizer() const { return fMinimizer; }

private:
   CombinedMinimizer fMinimizer;

private:
   // forbidden assignment operator
   MnMinimize &operator=(const MnMinimize &) { return *this; }
};

} // namespace Minuit2

} // namespace ROOT

#endif // ROOT_Minuit2_MnMinimize
