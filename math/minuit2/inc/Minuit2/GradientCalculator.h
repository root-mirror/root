// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Minuit2_GradientCalculator
#define ROOT_Minuit2_GradientCalculator

namespace ROOT {

namespace Minuit2 {

class MinimumParameters;
class FunctionGradient;

/**
   interface class for gradient calculators
 */
class GradientCalculator {

public:
   virtual ~GradientCalculator() {}

   virtual FunctionGradient operator()(const MinimumParameters &) const = 0;

   virtual FunctionGradient operator()(const MinimumParameters &, const FunctionGradient &) const = 0;
};

} // namespace Minuit2

} // namespace ROOT

#endif // ROOT_Minuit2_GradientCalculator
