// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Minuit2_FCNBase
#define ROOT_Minuit2_FCNBase

#include "Minuit2/MnConfig.h"

#include <vector>

#include "Minuit2/GenericFunction.h"

namespace ROOT {

namespace Minuit2 {

/**

  \defgroup Minuit Minuit2 Minimization Library

  New Object-oriented implementation of the MINUIT minimization package.
  More information is available at the home page of the \ref Minuit2Page "Minuit2" minimization package".

  \ingroup Math
*/

//______________________________________________________________________________
/**


Interface (abstract class) defining the function to be minimized, which has to be implemented by the user.

@author Fred James and Matthias Winkler; modified by Andras Zsenei and Lorenzo Moneta

\ingroup Minuit

 */

class FCNBase : public GenericFunction {

public:
   virtual ~FCNBase() {}

   /**

      The meaning of the vector of parameters is of course defined by the user,
      who uses the values of those parameters to calculate their function Value.
      The order and the position of these parameters is strictly the one specified
      by the user when supplying the starting values for minimization. The starting
      values must be specified by the user, either via an std::vector<double> or the
      MnUserParameters supplied as input to the MINUIT minimizers such as
      VariableMetricMinimizer or MnMigrad. Later values are determined by MINUIT
      as it searches for the Minimum or performs whatever analysis is requested by
      the user.

      @param v function parameters as defined by the user.

      @return the Value of the function.

      @see MnUserParameters
      @see VariableMetricMinimizer
      @see MnMigrad

   */

   virtual double operator()(const std::vector<double> &v) const = 0;

   /**

      Error definition of the function. MINUIT defines Parameter errors as the
      change in Parameter Value required to change the function Value by up. Normally,
      for chisquared fits it is 1, and for negative log likelihood, its Value is 0.5.
      If the user wants instead the 2-sigma errors for chisquared fits, it becomes 4,
      as Chi2(x+n*sigma) = Chi2(x) + n*n.

      Comment a little bit better with links!!!!!!!!!!!!!!!!!

   */

   virtual double ErrorDef() const { return Up(); }

   /**

      Error definition of the function. MINUIT defines Parameter errors as the
      change in Parameter Value required to change the function Value by up. Normally,
      for chisquared fits it is 1, and for negative log likelihood, its Value is 0.5.
      If the user wants instead the 2-sigma errors for chisquared fits, it becomes 4,
      as Chi2(x+n*sigma) = Chi2(x) + n*n.

      \todo Comment a little bit better with links!!!!!!!!!!!!!!!!! Idem for ErrorDef()

   */

   virtual double Up() const = 0;

   /**
       add interface to set dynamically a new error definition
       Re-implement this function if needed.
   */
   virtual void SetErrorDef(double){};
};

} // namespace Minuit2

} // namespace ROOT

#endif // ROOT_Minuit2_FCNBase
