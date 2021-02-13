// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Minuit2_FumiliBuilder
#define ROOT_Minuit2_FumiliBuilder

#include "Minuit2/MinimumBuilder.h"
#include "Minuit2/VariableMetricEDMEstimator.h"
#include "Minuit2/FumiliErrorUpdator.h"
#include "Minuit2/MnFcn.h"
#include "Minuit2/FunctionMinimum.h"

#include <vector>

namespace ROOT {

namespace Minuit2 {

/**

Builds the FunctionMinimum using the Fumili method.

@author Andras Zsenei, Creation date: 29 Sep 2004

@see <A HREF="http://www.cern.ch/winkler/minuit/tutorial/mntutorial.pdf">MINUIT Tutorial</A> on function minimization,
section 5

@ingroup Minuit

\todo the role of the strategy in Fumili

*/

class FumiliBuilder : public MinimumBuilder {

public:
   FumiliBuilder() : fEstimator(VariableMetricEDMEstimator()), fErrorUpdator(FumiliErrorUpdator()) {}

   ~FumiliBuilder() {}

   /**

   Class the member function calculating the Minimum and verifies the result
   depending on the strategy.

   @param fMnFcn the function to be minimized.

   @param fGradienCalculator not used in Fumili.

   @param fMinimumSeed the seed generator.

   @param fMnStrategy the strategy describing the number of function calls
   allowed for Gradient calculations.

   @param maxfcn maximum number of function calls after which the calculation
   will be stopped even if it has not yet converged.

   @param edmval expected vertical distance to the Minimum.

   @return Returns the function Minimum found.


   \todo Complete the documentation by understanding what is the reason to
   have two Minimum methods.

   */

   virtual FunctionMinimum Minimum(const MnFcn &fMnFcn, const GradientCalculator &fGradienCalculator,
                                   const MinimumSeed &fMinimumSeed, const MnStrategy &fMnStrategy, unsigned int maxfcn,
                                   double edmval) const;

   /**

   Calculates the Minimum based on the Fumili method

   @param fMnFcn the function to be minimized.

   @param fGradienCalculator not used in Fumili

   @param fMinimumSeed the seed generator.

   @param states vector containing the state result of each iteration

   @param maxfcn maximum number of function calls after which the calculation
   will be stopped even if it has not yet converged.

   @param edmval expected vertical distance to the Minimum

   @return Returns the function Minimum found.

   @see <A HREF="http://www.cern.ch/winkler/minuit/tutorial/mntutorial.pdf">MINUIT Tutorial</A> on function
   minimization, section 5


   \todo some nice Latex based formula here...

   */

   FunctionMinimum Minimum(const MnFcn &fMnFcn, const GradientCalculator &fGradienCalculator,
                           const MinimumSeed &fMinimumSeed, std::vector<MinimumState> &states, unsigned int maxfcn,
                           double edmval) const;

   /**

   Accessor to the EDM (expected vertical distance to the Minimum) estimator.

   @return The EDM estimator used in the builder.

   \todo Maybe a little explanation concerning EDM in all relevant classes.

   */

   const VariableMetricEDMEstimator &Estimator() const { return fEstimator; }

   /**

   Accessor to the Error updator of the builder.

   @return The FumiliErrorUpdator used by the FumiliBuilder.

   */

   const FumiliErrorUpdator &ErrorUpdator() const { return fErrorUpdator; }

private:
   VariableMetricEDMEstimator fEstimator;
   FumiliErrorUpdator fErrorUpdator;
};

} // namespace Minuit2

} // namespace ROOT

#endif // ROOT_Minuit2_FumiliBuilder
