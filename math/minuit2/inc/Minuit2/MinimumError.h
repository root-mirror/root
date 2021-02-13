// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Minuit2_MinimumError
#define ROOT_Minuit2_MinimumError

#include "Minuit2/BasicMinimumError.h"

#include <memory>

namespace ROOT {

namespace Minuit2 {

/** MinimumError keeps the inv. 2nd derivative (inv. Hessian) used for
    calculating the Parameter step size (-V*g) and for the covariance Update
    (ErrorUpdator). The covariance matrix is equal to twice the inv. Hessian.
 */

class MinimumError {

public:
   using MnHesseFailed = BasicMinimumError::MnHesseFailed;
   using MnInvertFailed = BasicMinimumError::MnInvertFailed;
   using MnMadePosDef = BasicMinimumError::MnMadePosDef;
   using MnNotPosDef = BasicMinimumError::MnNotPosDef;

public:
   MinimumError(unsigned int n) : fData(std::make_shared<BasicMinimumError>(n)) {}

   MinimumError(const MnAlgebraicSymMatrix &mat, double dcov) : fData(std::make_shared<BasicMinimumError>(mat, dcov)) {}

   MinimumError(const MnAlgebraicSymMatrix &mat, MnHesseFailed)
      : fData(std::make_shared<BasicMinimumError>(mat, MnHesseFailed{}))
   {
   }

   MinimumError(const MnAlgebraicSymMatrix &mat, MnMadePosDef)
      : fData(std::make_shared<BasicMinimumError>(mat, MnMadePosDef{}))
   {
   }

   MinimumError(const MnAlgebraicSymMatrix &mat, MnInvertFailed)
      : fData(std::make_shared<BasicMinimumError>(mat, MnInvertFailed{}))
   {
   }

   MinimumError(const MnAlgebraicSymMatrix &mat, MnNotPosDef)
      : fData(std::make_shared<BasicMinimumError>(mat, MnNotPosDef{}))
   {
   }

   MnAlgebraicSymMatrix Matrix() const { return fData->Matrix(); }

   const MnAlgebraicSymMatrix &InvHessian() const { return fData->InvHessian(); }

   MnAlgebraicSymMatrix Hessian() const { return fData->Hessian(); }

   double Dcovar() const { return fData->Dcovar(); }
   bool IsAccurate() const { return fData->IsAccurate(); }
   bool IsValid() const { return fData->IsValid(); }
   bool IsPosDef() const { return fData->IsPosDef(); }
   bool IsMadePosDef() const { return fData->IsMadePosDef(); }
   bool HesseFailed() const { return fData->HesseFailed(); }
   bool InvertFailed() const { return fData->InvertFailed(); }
   bool IsAvailable() const { return fData->IsAvailable(); }

private:
   std::shared_ptr<BasicMinimumError> fData;
};

} // namespace Minuit2

} // namespace ROOT

#endif // ROOT_Minuit2_MinimumError
