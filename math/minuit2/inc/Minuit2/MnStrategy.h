// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Minuit2_MnStrategy
#define ROOT_Minuit2_MnStrategy

namespace ROOT {

namespace Minuit2 {

//_________________________________________________________________________
/**
    API class for defining three levels of strategies: low (0), medium (1),
    high (>=2);
    acts on: Migrad (behavioural),
             Minos (lowers strategy by 1 for Minos-own minimization),
             Hesse (iterations),
             Numerical2PDerivative (iterations)
 */

class MnStrategy {

public:
   // default strategy
   MnStrategy();

   // user defined strategy (0, 1, >=2)
   explicit MnStrategy(unsigned int);

   ~MnStrategy() {}

   unsigned int Strategy() const { return fStrategy; }

   unsigned int GradientNCycles() const { return fGradNCyc; }
   double GradientStepTolerance() const { return fGradTlrStp; }
   double GradientTolerance() const { return fGradTlr; }

   unsigned int HessianNCycles() const { return fHessNCyc; }
   double HessianStepTolerance() const { return fHessTlrStp; }
   double HessianG2Tolerance() const { return fHessTlrG2; }
   unsigned int HessianGradientNCycles() const { return fHessGradNCyc; }

   int StorageLevel() const { return fStoreLevel; }

   bool IsLow() const { return fStrategy == 0; }
   bool IsMedium() const { return fStrategy == 1; }
   bool IsHigh() const { return fStrategy >= 2; }

   void SetLowStrategy();
   void SetMediumStrategy();
   void SetHighStrategy();

   void SetGradientNCycles(unsigned int n) { fGradNCyc = n; }
   void SetGradientStepTolerance(double stp) { fGradTlrStp = stp; }
   void SetGradientTolerance(double toler) { fGradTlr = toler; }

   void SetHessianNCycles(unsigned int n) { fHessNCyc = n; }
   void SetHessianStepTolerance(double stp) { fHessTlrStp = stp; }
   void SetHessianG2Tolerance(double toler) { fHessTlrG2 = toler; }
   void SetHessianGradientNCycles(unsigned int n) { fHessGradNCyc = n; }

   // set storage level of iteration quantities
   // 0 = store only last iterations 1 = full storage (default)
   void SetStorageLevel(unsigned int level) { fStoreLevel = level; }

private:
   unsigned int fStrategy;

   unsigned int fGradNCyc;
   double fGradTlrStp;
   double fGradTlr;
   unsigned int fHessNCyc;
   double fHessTlrStp;
   double fHessTlrG2;
   unsigned int fHessGradNCyc;
   int fStoreLevel;
};

} // namespace Minuit2

} // namespace ROOT

#endif // ROOT_Minuit2_MnStrategy
