// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Minuit2_FunctionGradient
#define ROOT_Minuit2_FunctionGradient

#include "Minuit2/BasicFunctionGradient.h"

#include <memory>

namespace ROOT {

namespace Minuit2 {

class FunctionGradient {

private:
public:
   explicit FunctionGradient(unsigned int n) : fData(std::make_shared<BasicFunctionGradient>(n)) {}

   // HD: this deep-copies, inconsistent to assignment?
   explicit FunctionGradient(const MnAlgebraicVector &grd) : fData(std::make_shared<BasicFunctionGradient>(grd)) {}

   FunctionGradient(const MnAlgebraicVector &grd, const MnAlgebraicVector &g2, const MnAlgebraicVector &gstep)
      : fData(std::make_shared<BasicFunctionGradient>(grd, g2, gstep))
   {
   }

   FunctionGradient(const FunctionGradient &grad) : fData(grad.fData) {}

   // HD: assignment shares the pointer
   FunctionGradient &operator=(const FunctionGradient & /*grad*/) = default;

   const MnAlgebraicVector &Grad() const { return fData->Grad(); }
   const MnAlgebraicVector &Vec() const { return fData->Vec(); }
   bool IsValid() const { return fData->IsValid(); }

   bool IsAnalytical() const { return fData->IsAnalytical(); }
   const MnAlgebraicVector &G2() const { return fData->G2(); }
   const MnAlgebraicVector &Gstep() const { return fData->Gstep(); }

private:
   std::shared_ptr<BasicFunctionGradient> fData;
};

} // namespace Minuit2

} // namespace ROOT

#endif // ROOT_Minuit2_FunctionGradient
