// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Minuit2_MinimumState
#define ROOT_Minuit2_MinimumState

#include "Minuit2/BasicMinimumState.h"

#include <memory>

namespace ROOT {

namespace Minuit2 {

class MinimumParameters;
class MinimumError;
class FunctionGradient;

/** MinimumState keeps the information (position, Gradient, 2nd deriv, etc)
    after one minimization step (usually in MinimumBuilder).
 */

class MinimumState {

public:
   /** invalid state */
   MinimumState(unsigned int n) : fData(std::make_shared<BasicMinimumState>(n, 0., 0., 0.)) {}
   /** state without parameters and errors (only function value an, edm and nfcn) */
   MinimumState(double fval, double edm, int nfcn) : fData(std::make_shared<BasicMinimumState>(0, fval, edm, nfcn)) {}
   /** state with parameters only (from stepping methods like Simplex, Scan) */
   MinimumState(const MinimumParameters &states, double edm, int nfcn)
      : fData(std::make_shared<BasicMinimumState>(states, edm, nfcn))
   {
   }

   /** state with parameters, Gradient and covariance (from Gradient methods
       such as Migrad) */
   MinimumState(const MinimumParameters &states, const MinimumError &err, const FunctionGradient &grad, double edm,
                int nfcn)
      : fData(std::make_shared<BasicMinimumState>(states, err, grad, edm, nfcn))
   {
   }

   const MinimumParameters &Parameters() const { return fData->Parameters(); }
   const MnAlgebraicVector &Vec() const { return fData->Vec(); }
   int size() const { return fData->size(); }

   const MinimumError &Error() const { return fData->Error(); }
   const FunctionGradient &Gradient() const { return fData->Gradient(); }
   double Fval() const { return fData->Fval(); }
   double Edm() const { return fData->Edm(); }
   int NFcn() const { return fData->NFcn(); }

   bool IsValid() const { return fData->IsValid(); }

   bool HasParameters() const { return fData->HasParameters(); }
   bool HasCovariance() const { return fData->HasCovariance(); }

private:
   std::shared_ptr<BasicMinimumState> fData;
};

} // namespace Minuit2

} // namespace ROOT

#endif // ROOT_Minuit2_MinimumState
