// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Minuit2_SimplexSeedGenerator
#define ROOT_Minuit2_SimplexSeedGenerator

#include "Minuit2/MinimumSeedGenerator.h"

namespace ROOT {

namespace Minuit2 {

class MinimumSeed;
class MnFcn;
class MnUserParameterState;
class MnStrategy;

/**
   generate Simplex starting point (state)
 */
class SimplexSeedGenerator : public MinimumSeedGenerator {

public:
   SimplexSeedGenerator() {}

   ~SimplexSeedGenerator() {}

   virtual MinimumSeed
   operator()(const MnFcn &, const GradientCalculator &, const MnUserParameterState &, const MnStrategy &) const;

   virtual MinimumSeed operator()(const MnFcn &, const AnalyticalGradientCalculator &, const MnUserParameterState &,
                                  const MnStrategy &) const;

private:
};

} // namespace Minuit2

} // namespace ROOT

#endif // ROOT_Minuit2_SimplexSeedGenerator
