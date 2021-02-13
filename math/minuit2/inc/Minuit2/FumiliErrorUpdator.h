// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Minuit2_FumiliErrorUpdator
#define ROOT_Minuit2_FumiliErrorUpdator

#include "Minuit2/MinimumErrorUpdator.h"

namespace ROOT {

namespace Minuit2 {

class MinimumState;
class MinimumParameters;
class GradientCalculator;
class FumiliFCNBase;
class FunctionGradient;

/**

In the case of the Fumili algorithm the Error matrix (or the Hessian
matrix containing the (approximate) second derivatives) is calculated
using a linearization of the model function negleting second
derivatives. (In some sense the Name Updator is a little bit misleading
as the Error matrix is not calculated by iteratively updating, like
in Davidon's or other similar variable metric methods, but by
recalculating each time).


@author  Andras Zsenei and Lorenzo Moneta, Creation date: 28 Sep 2004

@see <A HREF="http://www.cern.ch/winkler/minuit/tutorial/mntutorial.pdf">MINUIT Tutorial</A> on function minimization,
section 5

@see DavidonErrorUpdator

@ingroup Minuit

*/

class FumiliErrorUpdator : public MinimumErrorUpdator {

public:
   FumiliErrorUpdator() {}

   ~FumiliErrorUpdator() {}

   /**

   Member function that calculates the Error matrix (or the Hessian
   matrix containing the (approximate) second derivatives) using a
   linearization of the model function negleting second derivatives.

   @param fMinimumState used to calculate the change in the covariance
   matrix between the two iterations

   @param fMinimumParameters the parameters at the present iteration

   @param fGradientCalculator the Gradient calculator used to retrieved the Parameter transformation

   @param lambda the Marquard lambda factor


   \todo Some nice latex mathematical formuli...

   */

   virtual MinimumError Update(const MinimumState &fMinimumState, const MinimumParameters &fMinimumParameters,
                               const GradientCalculator &fGradientCalculator, double lambda) const;

   /**

   Member function which is only present due to the design already in place
   of the software. As all classes calculating the Error matrix are supposed
   inherit from the MinimumErrorUpdator they must inherit this method. In some
   methods calculating the aforementioned matrix some of these parameters are
   not needed and other parameters are necessary... Hopefully, a more elegant
   solution will be found in the future.

   \todo How to get rid of this dummy method which is only due to the inheritance

   */

   virtual MinimumError Update(const MinimumState &, const MinimumParameters &, const FunctionGradient &) const;

private:
};

} // namespace Minuit2

} // namespace ROOT

#endif // ROOT_Minuit2_FumiliErrorUpdator
