// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Minuit2_BasicFunctionGradient
#define ROOT_Minuit2_BasicFunctionGradient

#include "Minuit2/MnMatrix.h"

#include "Minuit2/StackAllocator.h"

namespace ROOT {

namespace Minuit2 {

// extern StackAllocator gStackAllocator;

class BasicFunctionGradient {

private:
public:
   explicit BasicFunctionGradient(unsigned int n)
      : fGradient(MnAlgebraicVector(n)), fG2ndDerivative(MnAlgebraicVector(n)), fGStepSize(MnAlgebraicVector(n)),
        fValid(false), fAnalytical(false)
   {
   }

   explicit BasicFunctionGradient(const MnAlgebraicVector &grd)
      : fGradient(grd), fG2ndDerivative(MnAlgebraicVector(grd.size())), fGStepSize(MnAlgebraicVector(grd.size())),
        fValid(true), fAnalytical(true)
   {
   }

   BasicFunctionGradient(const MnAlgebraicVector &grd, const MnAlgebraicVector &g2, const MnAlgebraicVector &gstep)
      : fGradient(grd), fG2ndDerivative(g2), fGStepSize(gstep), fValid(true), fAnalytical(false)
   {
   }

   ~BasicFunctionGradient() {}

   BasicFunctionGradient(const BasicFunctionGradient &grad)
      : fGradient(grad.fGradient), fG2ndDerivative(grad.fG2ndDerivative), fGStepSize(grad.fGStepSize),
        fValid(grad.fValid)
   {
   }

   BasicFunctionGradient &operator=(const BasicFunctionGradient &grad)
   {
      fGradient = grad.fGradient;
      fG2ndDerivative = grad.fG2ndDerivative;
      fGStepSize = grad.fGStepSize;
      fValid = grad.fValid;
      return *this;
   }

   void *operator new(size_t nbytes) { return StackAllocatorHolder::Get().Allocate(nbytes); }

   void operator delete(void *p, size_t /*nbytes */) { StackAllocatorHolder::Get().Deallocate(p); }

   const MnAlgebraicVector &Grad() const { return fGradient; }
   const MnAlgebraicVector &Vec() const { return fGradient; }
   bool IsValid() const { return fValid; }

   bool IsAnalytical() const { return fAnalytical; }
   const MnAlgebraicVector &G2() const { return fG2ndDerivative; }
   const MnAlgebraicVector &Gstep() const { return fGStepSize; }

private:
   MnAlgebraicVector fGradient;
   MnAlgebraicVector fG2ndDerivative;
   MnAlgebraicVector fGStepSize;
   bool fValid;
   bool fAnalytical;
};

} // namespace Minuit2

} // namespace ROOT

#endif // ROOT_Minuit2_BasicFunctionGradient
