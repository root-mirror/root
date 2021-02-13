// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Minuit2_MinimumErrorUpdator
#define ROOT_Minuit2_MinimumErrorUpdator

namespace ROOT {

namespace Minuit2 {

class MinimumState;
class MinimumError;
class MinimumParameters;
class FunctionGradient;

class MinimumErrorUpdator {

public:
   virtual ~MinimumErrorUpdator() {}

   virtual MinimumError Update(const MinimumState &, const MinimumParameters &, const FunctionGradient &) const = 0;
};

} // namespace Minuit2

} // namespace ROOT

#endif // ROOT_Minuit2_MinimumErrorUpdator
