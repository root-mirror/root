// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#include "Minuit2/MnStrategy.h"

namespace ROOT {

namespace Minuit2 {

MnStrategy::MnStrategy() : fStoreLevel(1)
{
   // default strategy
   SetMediumStrategy();
}

MnStrategy::MnStrategy(unsigned int stra) : fStoreLevel(1)
{
   // user defined strategy (0, 1, >=2)
   if (stra == 0)
      SetLowStrategy();
   else if (stra == 1)
      SetMediumStrategy();
   else
      SetHighStrategy();
}

void MnStrategy::SetLowStrategy()
{
   // set low strategy (0) values
   fStrategy = 0;
   SetGradientNCycles(2);
   SetGradientStepTolerance(0.5);
   SetGradientTolerance(0.1);
   SetHessianNCycles(3);
   SetHessianStepTolerance(0.5);
   SetHessianG2Tolerance(0.1);
   SetHessianGradientNCycles(1);
}

void MnStrategy::SetMediumStrategy()
{
   // set minimum strategy (1) the default
   fStrategy = 1;
   SetGradientNCycles(3);
   SetGradientStepTolerance(0.3);
   SetGradientTolerance(0.05);
   SetHessianNCycles(5);
   SetHessianStepTolerance(0.3);
   SetHessianG2Tolerance(0.05);
   SetHessianGradientNCycles(2);
}

void MnStrategy::SetHighStrategy()
{
   // set high strategy (2)
   fStrategy = 2;
   SetGradientNCycles(5);
   SetGradientStepTolerance(0.1);
   SetGradientTolerance(0.02);
   SetHessianNCycles(7);
   SetHessianStepTolerance(0.1);
   SetHessianG2Tolerance(0.02);
   SetHessianGradientNCycles(6);
}

} // namespace Minuit2

} // namespace ROOT
