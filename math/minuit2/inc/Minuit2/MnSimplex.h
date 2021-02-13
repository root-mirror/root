// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Minuit2_MnSimplex
#define ROOT_Minuit2_MnSimplex

#include "Minuit2/MnApplication.h"
#include "Minuit2/SimplexMinimizer.h"

#include <vector>

namespace ROOT {

namespace Minuit2 {

class FCNBase;

//_________________________________________________________________________
/**
    API class for minimization using the Simplex method, which does not need and use
    the derivatives of the function, but only function values.
    More information on the minimization method is available
    <A HREF="http://seal.web.cern.ch/mathlibs/documents/minuit/mntutorial.pdf">here</A>.

    It allows for user interaction: set/change parameters, do minimization,
    change parameters, re-do minimization etc.;
 */

class MnSimplex : public MnApplication {

public:
   /// construct from FCNBase + std::vector for parameters and errors
   MnSimplex(const FCNBase &fcn, const std::vector<double> &par, const std::vector<double> &err, unsigned int stra = 1)
      : MnApplication(fcn, MnUserParameterState(par, err), MnStrategy(stra)), fMinimizer(SimplexMinimizer())
   {
   }

   /// construct from FCNBase + std::vector for parameters and covariance
   MnSimplex(const FCNBase &fcn, const std::vector<double> &par, unsigned int nrow, const std::vector<double> &cov,
             unsigned int stra = 1)
      : MnApplication(fcn, MnUserParameterState(par, cov, nrow), MnStrategy(stra)), fMinimizer(SimplexMinimizer())
   {
   }

   /// construct from FCNBase + std::vector for parameters and MnUserCovariance
   MnSimplex(const FCNBase &fcn, const std::vector<double> &par, const MnUserCovariance &cov, unsigned int stra = 1)
      : MnApplication(fcn, MnUserParameterState(par, cov), MnStrategy(stra)), fMinimizer(SimplexMinimizer())
   {
   }

   /// construct from FCNBase + MnUserParameters
   MnSimplex(const FCNBase &fcn, const MnUserParameters &par, unsigned int stra = 1)
      : MnApplication(fcn, MnUserParameterState(par), MnStrategy(stra)), fMinimizer(SimplexMinimizer())
   {
   }

   /// construct from FCNBase + MnUserParameters + MnUserCovariance
   MnSimplex(const FCNBase &fcn, const MnUserParameters &par, const MnUserCovariance &cov, unsigned int stra = 1)
      : MnApplication(fcn, MnUserParameterState(par, cov), MnStrategy(stra)), fMinimizer(SimplexMinimizer())
   {
   }

   /// construct from FCNBase + MnUserParameterState + MnStrategy
   MnSimplex(const FCNBase &fcn, const MnUserParameterState &par, const MnStrategy &str)
      : MnApplication(fcn, MnUserParameterState(par), str), fMinimizer(SimplexMinimizer())
   {
   }

   MnSimplex(const MnSimplex &migr)
      : MnApplication(migr.Fcnbase(), migr.State(), migr.Strategy(), migr.NumOfCalls()), fMinimizer(migr.fMinimizer)
   {
   }

   ~MnSimplex() {}

   ModularFunctionMinimizer &Minimizer() { return fMinimizer; }
   const ModularFunctionMinimizer &Minimizer() const { return fMinimizer; }

private:
   SimplexMinimizer fMinimizer;

private:
   // forbidden assignment of migrad (const FCNBase& = )
   MnSimplex &operator=(const MnSimplex &) { return *this; }
};

} // namespace Minuit2

} // namespace ROOT

#endif // ROOT_Minuit2_MnSimplex
