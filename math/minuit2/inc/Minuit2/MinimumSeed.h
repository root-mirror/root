// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Minuit2_MinimumSeed
#define ROOT_Minuit2_MinimumSeed

#include "Minuit2/BasicMinimumSeed.h"

#include <memory>

namespace ROOT {

namespace Minuit2 {

class MinimumState;
class MinimumParameters;
class MinimumError;
class FunctionGradient;
class MnUserTransformation;

/** MinimumSeed contains the starting values for the minimization produced
    by the SeedGenerator.
 */

class MinimumSeed {

public:
   MinimumSeed(const MinimumState &st, const MnUserTransformation &trafo)
      : fData(std::make_shared<BasicMinimumSeed>(st, trafo))
   {
   }

   const MinimumState &State() const { return fData->State(); }
   const MinimumParameters &Parameters() const { return fData->Parameters(); }
   const MinimumError &Error() const { return fData->Error(); }
   const FunctionGradient &Gradient() const { return fData->Gradient(); }
   const MnUserTransformation &Trafo() const { return fData->Trafo(); }
   const MnMachinePrecision &Precision() const { return fData->Precision(); }
   double Fval() const { return fData->Fval(); }
   double Edm() const { return fData->Edm(); }
   unsigned int NFcn() const { return fData->NFcn(); }
   bool IsValid() const { return fData->IsValid(); }

private:
   std::shared_ptr<BasicMinimumSeed> fData;
};

} // namespace Minuit2

} // namespace ROOT

#endif // ROOT_Minuit2_MinimumSeed
