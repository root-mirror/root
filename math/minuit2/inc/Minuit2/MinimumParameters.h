// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Minuit2_MinimumParameters
#define ROOT_Minuit2_MinimumParameters

#include "Minuit2/BasicMinimumParameters.h"

#include <memory>

namespace ROOT {

namespace Minuit2 {

class MinimumParameters {

public:
   MinimumParameters(unsigned int n, double fval = 0) : fData(std::make_shared<BasicMinimumParameters>(n, fval)) {}

   /** takes the Parameter vector */
   MinimumParameters(const MnAlgebraicVector &avec, double fval)
      : fData(std::make_shared<BasicMinimumParameters>(avec, fval))
   {
   }

   /** takes the Parameter vector plus step size x1 - x0 = dirin */
   MinimumParameters(const MnAlgebraicVector &avec, const MnAlgebraicVector &dirin, double fval)
      : fData(std::make_shared<BasicMinimumParameters>(avec, dirin, fval))
   {
   }

   const MnAlgebraicVector &Vec() const { return fData->Vec(); }
   const MnAlgebraicVector &Dirin() const { return fData->Dirin(); }
   double Fval() const { return fData->Fval(); }
   bool IsValid() const { return fData->IsValid(); }
   bool HasStepSize() const { return fData->HasStepSize(); }

private:
   std::shared_ptr<BasicMinimumParameters> fData;
};

} // namespace Minuit2

} // namespace ROOT

#endif // ROOT_Minuit2_MinimumParameters
