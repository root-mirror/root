// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Minuit2_NegativeG2LineSearch
#define ROOT_Minuit2_NegativeG2LineSearch

namespace ROOT {

namespace Minuit2 {

class MnFcn;
class MinimumState;
class GradientCalculator;
class MnMachinePrecision;
class FunctionGradient;

/** In case that one of the components of the second derivative g2 calculated
    by the numerical Gradient calculator is negative, a 1dim line search in
    the direction of that component is done in order to find a better position
    where g2 is again positive.
 */

class NegativeG2LineSearch {

public:
   NegativeG2LineSearch() {}

   ~NegativeG2LineSearch() {}

   MinimumState
   operator()(const MnFcn &, const MinimumState &, const GradientCalculator &, const MnMachinePrecision &) const;

   bool HasNegativeG2(const FunctionGradient &, const MnMachinePrecision &) const;

private:
};

} // namespace Minuit2

} // namespace ROOT

#endif // ROOT_Minuit2_NegativeG2LineSearch
