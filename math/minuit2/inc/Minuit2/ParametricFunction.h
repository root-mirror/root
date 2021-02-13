// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Minuit2_ParametricFunction
#define ROOT_Minuit2_ParametricFunction

#include "Minuit2/MnConfig.h"
#include <vector>
#include <cassert>

#include "Minuit2/FCNBase.h"

namespace ROOT {

namespace Minuit2 {

/**

Function which has parameters. For example, one could define
a one-dimensional Gaussian, by considering x as an input coordinate
for the evaluation of the function, and the mean and the square root
of the variance as parameters.
<p>
AS OF NOW PARAMETRICFUNCTION INHERITS FROM FCNBASE INSTEAD OF
GENERICFUNCTION. THIS IS ONLY BECAUSE NUMERICAL2PGRADIENTCALCULATOR
NEEDS AN FCNBASE OBJECT AND WILL BE CHANGED!!!!!!!!!!!!!!!!

@ingroup Minuit

\todo ParametricFunction and all the classes that inherit from it
are inheriting also FCNBase so that the Gradient calculation has
the Up() member function. That is not really good...


 */

class ParametricFunction : public FCNBase {

public:
   /**

   Constructor which initializes the ParametricFunction with the
   parameters given as input.

   @param params vector containing the initial Parameter values

   */

   ParametricFunction(const std::vector<double> &params) : par(params) {}

   /**

   Constructor which initializes the ParametricFunction by setting
   the number of parameters.

   @param nparams number of parameters of the parametric function

   */

   ParametricFunction(int nparams) : par(nparams) {}

   virtual ~ParametricFunction() {}

   /**

   Sets the parameters of the ParametricFunction.

   @param params vector containing the Parameter values

   */

   virtual void SetParameters(const std::vector<double> &params) const
   {

      assert(params.size() == par.size());
      par = params;
   }

   /**

   Accessor for the state of the parameters.

   @return vector containing the present Parameter settings

   */

   virtual const std::vector<double> &GetParameters() const { return par; }

   /**

   Accessor for the number of  parameters.

   @return the number of function parameters

   */
   virtual unsigned int NumberOfParameters() const { return par.size(); }

   // Why do I need to declare it here, it should be inherited without
   // any problems, no?

   /**

   Evaluates the function with the given coordinates.

   @param x vector containing the input coordinates

   @return the result of the function evaluation with the given
   coordinates.

   */

   virtual double operator()(const std::vector<double> &x) const = 0;

   /**

   Evaluates the function with the given coordinates and Parameter
   values. This member function is useful to implement when speed
   is an issue as it is faster to call only one function instead
   of two (SetParameters and operator()). The default implementation,
   provided for convenience, does the latter.

   @param x vector containing the input coordinates

   @param params vector containing the Parameter values

   @return the result of the function evaluation with the given
   coordinates and parameters

   */

   virtual double operator()(const std::vector<double> &x, const std::vector<double> &params) const
   {
      SetParameters(params);
      return operator()(x);
   }

   /**

   Member function returning the Gradient of the function with respect
   to its variables (but without including gradients with respect to
   its internal parameters).

   @param x vector containing the coordinates of the point where the
   Gradient is to be calculated.

   @return the Gradient vector of the function at the given point.

   */

   virtual std::vector<double> GetGradient(const std::vector<double> &x) const;

protected:
   /**

   The vector containing the parameters of the function
   It is mutable for "historical reasons" as in the hierarchy
   methods and classes are const and all the implications of changing
   them back to non-const are not clear.

   */

   mutable std::vector<double> par;
};

} // namespace Minuit2

} // namespace ROOT

#endif // ROOT_Minuit2_ParametricFunction
