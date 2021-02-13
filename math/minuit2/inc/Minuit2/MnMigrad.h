// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Minuit2_MnMigrad
#define ROOT_Minuit2_MnMigrad

#include "Minuit2/MnApplication.h"
#include "Minuit2/VariableMetricMinimizer.h"

#include <vector>

namespace ROOT {

namespace Minuit2 {

class FCNBase;

//_____________________________________________________________________________
/**
   API class for minimization using Variable Metric technology ("MIGRAD");
    allows for user interaction: set/change parameters, do minimization,
    change parameters, re-do minimization etc.;
    also used by MnMinos and MnContours;
 */

class MnMigrad : public MnApplication {

public:
   /// construct from FCNBase + std::vector for parameters and errors
   MnMigrad(const FCNBase &fcn, const std::vector<double> &par, const std::vector<double> &err, unsigned int stra = 1)
      : MnApplication(fcn, MnUserParameterState(par, err), MnStrategy(stra)), fMinimizer(VariableMetricMinimizer())
   {
   }

   /// construct from FCNBase + std::vector for parameters and covariance
   MnMigrad(const FCNBase &fcn, const std::vector<double> &par, unsigned int nrow, const std::vector<double> &cov,
            unsigned int stra = 1)
      : MnApplication(fcn, MnUserParameterState(par, cov, nrow), MnStrategy(stra)),
        fMinimizer(VariableMetricMinimizer())
   {
   }

   /// construct from FCNBase + std::vector for parameters and MnUserCovariance
   MnMigrad(const FCNBase &fcn, const std::vector<double> &par, const MnUserCovariance &cov, unsigned int stra = 1)
      : MnApplication(fcn, MnUserParameterState(par, cov), MnStrategy(stra)), fMinimizer(VariableMetricMinimizer())
   {
   }

   /// construct from FCNBase + MnUserParameters
   MnMigrad(const FCNBase &fcn, const MnUserParameters &par, unsigned int stra = 1)
      : MnApplication(fcn, MnUserParameterState(par), MnStrategy(stra)), fMinimizer(VariableMetricMinimizer())
   {
   }

   /// construct from FCNBase + MnUserParameters + MnUserCovariance
   MnMigrad(const FCNBase &fcn, const MnUserParameters &par, const MnUserCovariance &cov, unsigned int stra = 1)
      : MnApplication(fcn, MnUserParameterState(par, cov), MnStrategy(stra)), fMinimizer(VariableMetricMinimizer())
   {
   }

   /// construct from FCNBase + MnUserParameterState + MnStrategy
   MnMigrad(const FCNBase &fcn, const MnUserParameterState &par, const MnStrategy &str)
      : MnApplication(fcn, MnUserParameterState(par), str), fMinimizer(VariableMetricMinimizer())
   {
   }

   // constructs from gradient FCN

   /// construct from FCNGradientBase + std::vector for parameters and errors
   MnMigrad(const FCNGradientBase &fcn, const std::vector<double> &par, const std::vector<double> &err,
            unsigned int stra = 1)
      : MnApplication(fcn, MnUserParameterState(par, err), MnStrategy(stra)), fMinimizer(VariableMetricMinimizer())
   {
   }

   /// construct from FCNGradientBase + std::vector for parameters and covariance
   MnMigrad(const FCNGradientBase &fcn, const std::vector<double> &par, unsigned int nrow,
            const std::vector<double> &cov, unsigned int stra = 1)
      : MnApplication(fcn, MnUserParameterState(par, cov, nrow), MnStrategy(stra)),
        fMinimizer(VariableMetricMinimizer())
   {
   }

   /// construct from FCNGradientBase + std::vector for parameters and MnUserCovariance
   MnMigrad(const FCNGradientBase &fcn, const std::vector<double> &par, const MnUserCovariance &cov,
            unsigned int stra = 1)
      : MnApplication(fcn, MnUserParameterState(par, cov), MnStrategy(stra)), fMinimizer(VariableMetricMinimizer())
   {
   }

   /// construct from FCNGradientBase + MnUserParameters
   MnMigrad(const FCNGradientBase &fcn, const MnUserParameters &par, unsigned int stra = 1)
      : MnApplication(fcn, MnUserParameterState(par), MnStrategy(stra)), fMinimizer(VariableMetricMinimizer())
   {
   }

   /// construct from FCNGradientBase + MnUserParameters + MnUserCovariance
   MnMigrad(const FCNGradientBase &fcn, const MnUserParameters &par, const MnUserCovariance &cov, unsigned int stra = 1)
      : MnApplication(fcn, MnUserParameterState(par, cov), MnStrategy(stra)), fMinimizer(VariableMetricMinimizer())
   {
   }

   /// construct from FCNGradientBase + MnUserParameterState + MnStrategy
   MnMigrad(const FCNGradientBase &fcn, const MnUserParameterState &par, const MnStrategy &str)
      : MnApplication(fcn, MnUserParameterState(par), str), fMinimizer(VariableMetricMinimizer())
   {
   }

   ~MnMigrad() {}

   /// Copy constructor, copy shares the reference to the same FCNBase in MnApplication
   MnMigrad(const MnMigrad &) = default;

   // Copy assignment deleted, since MnApplication has unassignable reference to FCNBase
   MnMigrad &operator=(const MnMigrad &) = delete;

   ModularFunctionMinimizer &Minimizer() { return fMinimizer; }
   const ModularFunctionMinimizer &Minimizer() const { return fMinimizer; }

private:
   VariableMetricMinimizer fMinimizer;
};

} // namespace Minuit2

} // namespace ROOT

#endif // ROOT_Minuit2_MnMigrad
