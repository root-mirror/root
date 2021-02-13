// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Minuit2_CombinedMinimizer
#define ROOT_Minuit2_CombinedMinimizer

#include "Minuit2/ModularFunctionMinimizer.h"
#include "Minuit2/MnSeedGenerator.h"
#include "Minuit2/CombinedMinimumBuilder.h"

namespace ROOT {

namespace Minuit2 {

//__________________________________________________________________________
/**
   Combined minimizer: combination of Migrad and Simplex. I
   If the Migrad method fails at first attempt, a simplex
   minimization is performed and then migrad is tried again.


*/

class CombinedMinimizer : public ModularFunctionMinimizer {

public:
   CombinedMinimizer() : fMinSeedGen(MnSeedGenerator()), fMinBuilder(CombinedMinimumBuilder()) {}

   ~CombinedMinimizer() {}

   const MinimumSeedGenerator &SeedGenerator() const { return fMinSeedGen; }
   const MinimumBuilder &Builder() const { return fMinBuilder; }
   MinimumBuilder &Builder() { return fMinBuilder; }

private:
   MnSeedGenerator fMinSeedGen;
   CombinedMinimumBuilder fMinBuilder;
};

} // namespace Minuit2

} // namespace ROOT

#endif // ROOT_Minuit2_CombinedMinimizer
