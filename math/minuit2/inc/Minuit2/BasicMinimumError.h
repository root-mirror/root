// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Minuit2_BasicMinimumError
#define ROOT_Minuit2_BasicMinimumError

#include "Minuit2/MnConfig.h"
#include "Minuit2/MnMatrix.h"
#include "Minuit2/LaSum.h"
#include "Minuit2/StackAllocator.h"

namespace ROOT {

namespace Minuit2 {

// extern StackAllocator gStackAllocator;

/**
   Internal Class containing the error information on the
   estimated minimum :
   Error matrix + dcovar + additional flags for quality and validity checks
 */

class BasicMinimumError {

public:
   class MnNotPosDef {
   };
   class MnMadePosDef {
   };
   class MnHesseFailed {
   };
   class MnInvertFailed {
   };

public:
   BasicMinimumError(unsigned int n)
      : fMatrix(MnAlgebraicSymMatrix(n)), fDCovar(1.), fValid(false), fPosDef(false), fMadePosDef(false),
        fHesseFailed(false), fInvertFailed(false), fAvailable(false)
   {
   }

   BasicMinimumError(const MnAlgebraicSymMatrix &mat, double dcov)
      : fMatrix(mat), fDCovar(dcov), fValid(true), fPosDef(true), fMadePosDef(false), fHesseFailed(false),
        fInvertFailed(false), fAvailable(true)
   {
   }

   BasicMinimumError(const MnAlgebraicSymMatrix &mat, MnHesseFailed)
      : fMatrix(mat), fDCovar(1.), fValid(false), fPosDef(false), fMadePosDef(false), fHesseFailed(true),
        fInvertFailed(false), fAvailable(true)
   {
   }

   BasicMinimumError(const MnAlgebraicSymMatrix &mat, MnMadePosDef)
      : fMatrix(mat), fDCovar(1.), fValid(true), fPosDef(false), fMadePosDef(true), fHesseFailed(false),
        fInvertFailed(false), fAvailable(true)
   {
   }

   BasicMinimumError(const MnAlgebraicSymMatrix &mat, MnInvertFailed)
      : fMatrix(mat), fDCovar(1.), fValid(false), fPosDef(true), fMadePosDef(false), fHesseFailed(false),
        fInvertFailed(true), fAvailable(true)
   {
   }

   BasicMinimumError(const MnAlgebraicSymMatrix &mat, MnNotPosDef)
      : fMatrix(mat), fDCovar(1.), fValid(false), fPosDef(false), fMadePosDef(false), fHesseFailed(false),
        fInvertFailed(false), fAvailable(true)
   {
   }

   ~BasicMinimumError() {}

   BasicMinimumError(const BasicMinimumError &e)
      : fMatrix(e.fMatrix), fDCovar(e.fDCovar), fValid(e.fValid), fPosDef(e.fPosDef), fMadePosDef(e.fMadePosDef),
        fHesseFailed(e.fHesseFailed), fInvertFailed(e.fInvertFailed), fAvailable(e.fAvailable)
   {
   }

   BasicMinimumError &operator=(const BasicMinimumError &err)
   {
      fMatrix = err.fMatrix;
      fDCovar = err.fDCovar;
      fValid = err.fValid;
      fPosDef = err.fPosDef;
      fMadePosDef = err.fMadePosDef;
      fHesseFailed = err.fHesseFailed;
      fInvertFailed = err.fInvertFailed;
      fAvailable = err.fAvailable;
      return *this;
   }

   void *operator new(size_t nbytes) { return StackAllocatorHolder::Get().Allocate(nbytes); }

   void operator delete(void *p, size_t /*nbytes */) { StackAllocatorHolder::Get().Deallocate(p); }

   MnAlgebraicSymMatrix Matrix() const { return 2. * fMatrix; }

   const MnAlgebraicSymMatrix &InvHessian() const { return fMatrix; }

   MnAlgebraicSymMatrix Hessian() const;

   double Dcovar() const { return fDCovar; }
   bool IsAccurate() const { return fDCovar < 0.1; }
   bool IsValid() const { return fValid; }
   bool IsPosDef() const { return fPosDef; }
   bool IsMadePosDef() const { return fMadePosDef; }
   bool HesseFailed() const { return fHesseFailed; }
   bool InvertFailed() const { return fInvertFailed; }
   bool IsAvailable() const { return fAvailable; }

private:
   MnAlgebraicSymMatrix fMatrix;
   double fDCovar;
   bool fValid;
   bool fPosDef;
   bool fMadePosDef;
   bool fHesseFailed;
   bool fInvertFailed;
   bool fAvailable;
};

} // namespace Minuit2

} // namespace ROOT

#endif // ROOT_Minuit2_BasicMinimumError
