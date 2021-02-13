// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Minuit2_MnFumiliMinimize
#define ROOT_Minuit2_MnFumiliMinimize

#include "Minuit2/MnApplication.h"
#include "Minuit2/FumiliMinimizer.h"
#include "Minuit2/FumiliFCNBase.h"

#include <vector>

namespace ROOT {

namespace Minuit2 {

// class FumiliFCNBase;
// class FCNBase;

//___________________________________________________________________________
/**


API class for minimization using Fumili technology;
allows for user interaction: set/change parameters, do minimization,
change parameters, re-do minimization etc.;
also used by MnMinos and MnContours;


 */

class MnFumiliMinimize : public MnApplication {

public:
   /// construct from FumiliFCNBase + std::vector for parameters and errors
   MnFumiliMinimize(const FumiliFCNBase &fcn, const std::vector<double> &par, const std::vector<double> &err,
                    unsigned int stra = 1)
      : MnApplication(fcn, MnUserParameterState(par, err), MnStrategy(stra)), fMinimizer(FumiliMinimizer()), fFCN(fcn)
   {
   }

   /// construct from FumiliFCNBase + std::vector for parameters and covariance
   MnFumiliMinimize(const FumiliFCNBase &fcn, const std::vector<double> &par, unsigned int nrow,
                    const std::vector<double> &cov, unsigned int stra = 1)
      : MnApplication(fcn, MnUserParameterState(par, cov, nrow), MnStrategy(stra)), fMinimizer(FumiliMinimizer()),
        fFCN(fcn)
   {
   }

   /// construct from FumiliFCNBase + std::vector for parameters and MnUserCovariance
   MnFumiliMinimize(const FumiliFCNBase &fcn, const std::vector<double> &par, const MnUserCovariance &cov,
                    unsigned int stra = 1)
      : MnApplication(fcn, MnUserParameterState(par, cov), MnStrategy(stra)), fMinimizer(FumiliMinimizer()), fFCN(fcn)
   {
   }

   /// construct from FumiliFCNBase + MnUserParameters
   MnFumiliMinimize(const FumiliFCNBase &fcn, const MnUserParameters &par, unsigned int stra = 1)
      : MnApplication(fcn, MnUserParameterState(par), MnStrategy(stra)), fMinimizer(FumiliMinimizer()), fFCN(fcn)
   {
   }

   /// construct from FumiliFCNBase + MnUserParameters + MnUserCovariance
   MnFumiliMinimize(const FumiliFCNBase &fcn, const MnUserParameters &par, const MnUserCovariance &cov,
                    unsigned int stra = 1)
      : MnApplication(fcn, MnUserParameterState(par, cov), MnStrategy(stra)), fMinimizer(FumiliMinimizer()), fFCN(fcn)
   {
   }

   /// construct from FumiliFCNBase + MnUserParameterState + MnStrategy
   MnFumiliMinimize(const FumiliFCNBase &fcn, const MnUserParameterState &par, const MnStrategy &str)
      : MnApplication(fcn, MnUserParameterState(par), str), fMinimizer(FumiliMinimizer()), fFCN(fcn)
   {
   }

   MnFumiliMinimize(const MnFumiliMinimize &migr)
      : MnApplication(migr.Fcnbase(), migr.State(), migr.Strategy(), migr.NumOfCalls()), fMinimizer(migr.fMinimizer),
        fFCN(migr.Fcnbase())
   {
   }

   virtual ~MnFumiliMinimize() {}

   FumiliMinimizer &Minimizer() { return fMinimizer; }
   const FumiliMinimizer &Minimizer() const { return fMinimizer; }

   const FumiliFCNBase &Fcnbase() const { return fFCN; }

   /// overwrite Minimize to use FumiliFCNBase
   virtual FunctionMinimum operator()(unsigned int = 0, double = 0.1);

private:
   FumiliMinimizer fMinimizer;
   const FumiliFCNBase &fFCN;

private:
   // forbidden assignment of migrad (const FumiliFCNBase& = )
   MnFumiliMinimize &operator=(const MnFumiliMinimize &) { return *this; }
};

} // namespace Minuit2

} // namespace ROOT

#endif // ROOT_Minuit2_MnFumiliMinimize
