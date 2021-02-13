// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef MN_GaussianModelFunction_H_
#define MN_GaussianModelFunction_H_

#define _USE_MATH_DEFINES

#include "Minuit2/ParametricFunction.h"

#include "Minuit2/MnFcn.h"
#include "Minuit2/MnStrategy.h"
#include "Minuit2/MnUserParameterState.h"

#include <cmath>
#include <vector>
#include <cassert>

namespace ROOT {

namespace Minuit2 {

/**

Sample implementation of a parametric function. It can be used for
example for the Fumili method when minimizing with Minuit.
In the present case the function is a one-dimensional Gaussian,
which is described by its mean, standard deviation and the constant
term describing the amplitude. As it is used for
function minimization, the role of the variables (or coordinates) and
parameters is inversed! I.e. in the case of a one-dimensional
Gaussian it is x that will be the Parameter and the mean, standard
deviation etc will be the variables.

@author Andras Zsenei and Lorenzo Moneta, Creation date: 26 Oct 2004

@see <A HREF="http://mathworld.wolfram.com/NormalDistribution.html"> Definition of
the Normal/Gaussian distribution </A> (note: this Gaussian is normalized).

@see ParametricFunction

@see FumiliFCNBase

@see FumiliMaximumLikelihoodFCN

@ingroup Minuit

*/

class GaussianModelFunction : public ParametricFunction {

public:
   /**

   Constructor which initializes the normalized Gaussian with x = 0.0.

   */

   GaussianModelFunction() : ParametricFunction(1)
   {

      // setting some default values for the parameters
      std::vector<double> param;
      param.push_back(0.0);
      SetParameters(param);
   }

   /**

   Constructor which initializes the ParametricFunction with the
   parameters given as input.

   @param params vector containing the initial Parameter Value.

   */

   GaussianModelFunction(const std::vector<double> &params) : ParametricFunction(params) { assert(params.size() == 1); }

   ~GaussianModelFunction() {}

   /**

   Calculates the Gaussian as a function of the given input.

   @param x vector containing the mean, standard deviation and amplitude.

   @return the Value of the Gaussian for the given input.

   @see <A HREF="http://mathworld.wolfram.com/NormalDistribution.html"> Definition of
   the Normal/Gaussian distribution </A> (note: this Gaussian is normalized).

   */

   double operator()(const std::vector<double> &x) const
   {

      assert(x.size() == 3);
      // commented out for speed-up (even though that is the object-oriented
      // way to do things)
      // std::vector<double> par = GetParameters();

      constexpr double two_pi = 2 * 3.14159265358979323846; // M_PI is not standard

      return x[2] * std::exp(-0.5 * (par[0] - x[0]) * (par[0] - x[0]) / (x[1] * x[1])) /
             (std::sqrt(two_pi) * std::fabs(x[1]));
   }

   /**

   Calculates the Gaussian as a function of the given input.

   @param x vector containing the mean, the standard deviation and the constant
   describing the Gaussian.

   @param par vector containing the x coordinate (which is the Parameter in
   the case of a minimization).

   @return the Value of the Gaussian for the given input.

   @see <A HREF="http://mathworld.wolfram.com/NormalDistribution.html"> Definition of
   the Normal/Gaussian distribution </A> (note: this Gaussian is normalized).

   */

   double operator()(const std::vector<double> &x, const std::vector<double> &param) const
   {

      constexpr double two_pi = 2 * 3.14159265358979323846; // M_PI is not standard

      assert(param.size() == 1);
      assert(x.size() == 3);
      return x[2] * std::exp(-0.5 * (param[0] - x[0]) * (param[0] - x[0]) / (x[1] * x[1])) /
             (std::sqrt(two_pi) * std::fabs(x[1]));
   }

   /**

   THAT SHOULD BE REMOVED, IT IS ONLY HERE, BECAUSE AT PRESENT FOR GRADIENT
   CALCULATION ONE NEEDS TO INHERIT FROM FCNBASE WHICH NEEDS THIS METHOD

   */

   virtual double Up() const { return 1.0; }

   std::vector<double> GetGradient(const std::vector<double> &x) const
   {

      const std::vector<double> &param = GetParameters();
      assert(param.size() == 1);
      std::vector<double> grad(x.size());

      constexpr double two_pi = 2 * 3.14159265358979323846; // M_PI is not standard

      double y = (param[0] - x[0]) / x[1];
      double gaus = std::exp(-0.5 * y * y) / (std::sqrt(two_pi) * std::fabs(x[1]));

      grad[0] = y / (x[1]) * gaus * x[2];
      grad[1] = x[2] * gaus * (y * y - 1.0) / x[1];
      grad[2] = gaus;
      // std::cout << "GRADIENT" << y << "  " << gaus << "  " << x[0] << "  " << x[1] << "  " << grad[0] << "   " <<
      // grad[1] << std::endl;

      return grad;
   }
};

} // namespace Minuit2

} // namespace ROOT

#endif // MN_GaussianModelFunction_H_
