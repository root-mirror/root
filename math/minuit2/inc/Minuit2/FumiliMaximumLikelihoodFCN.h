// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Minuit2_FumiliMaximumLikelihoodFCN
#define ROOT_Minuit2_FumiliMaximumLikelihoodFCN

#include "FumiliFCNBase.h"
#include "Minuit2/ParametricFunction.h"
#include "Math/Util.h"
#include <vector>
#include <cassert>

namespace ROOT {

namespace Minuit2 {

/**

Extension of the FCNBase for the Fumili method. Fumili applies only to
minimization problems used for fitting. The method is based on a
linearization of the model function negleting second derivatives.
User needs to provide the model function. In this cased the function
to be minimized is the sum of the logarithms of the model function
for the different measurements times -1.


@author Andras Zsenei and Lorenzo Moneta, Creation date: 3 Sep 2004

@see <A HREF="http://www.cern.ch/winkler/minuit/tutorial/mntutorial.pdf">MINUIT Tutorial</A> on function minimization,
section 5

@see FumiliStandardMaximumLikelihoodFCN

@ingroup Minuit

\todo Insert a nice latex formula...

*/

class FumiliMaximumLikelihoodFCN : public FumiliFCNBase {

public:
   FumiliMaximumLikelihoodFCN() {}

   virtual ~FumiliMaximumLikelihoodFCN() {}

   /**

   Sets the model function for the data (for example gaussian+linear for a peak)

   @param modelFCN a reference to the model function.

   */

   void SetModelFunction(const ParametricFunction &modelFCN) { fModelFunction = &modelFCN; }

   /**

   Returns the model function used for the data.

   @return Returns a pointer to the model function.

   */

   const ParametricFunction *ModelFunction() const { return fModelFunction; }

   /**

   Evaluates the model function for the different measurement points and
   the Parameter values supplied, calculates a figure-of-merit for each
   measurement and returns a vector containing the result of this
   evaluation.

   @param par vector of Parameter values to feed to the model function.

   @return A vector containing the figures-of-merit for the model function evaluated
   for each set of measurements.

   */

   virtual std::vector<double> Elements(const std::vector<double> &par) const = 0;

   /**

   Accessor to the parameters of a given measurement. For example in the
   case of a chi-square fit with a one-dimensional Gaussian, the Parameter
   characterizing the measurement will be the position. It is the Parameter
   that is feeded to the model function.

   @param Index Index of the measueremnt the parameters of which to return
   @return A vector containing the values characterizing a measurement

   */

   virtual const std::vector<double> &GetMeasurement(int Index) const = 0;

   /**

   Accessor to the number of measurements used for calculating the
   present figure of merit.

   @return the number of measurements

   */

   virtual int GetNumberOfMeasurements() const = 0;

   /**

   Calculates the function for the maximum likelihood method. The user must
   implement in a class which inherits from FumiliChi2FCN the member function
   Elements() which will supply the Elements for the sum.


   @param par vector containing the Parameter values for the model function

   @return The sum of the natural logarithm of the Elements multiplied by -1

   @see FumiliFCNBase#elements

   */

   double operator()(const std::vector<double> &par) const
   {

      double sumoflogs = 0.0;
      std::vector<double> vecElements = Elements(par);
      unsigned int vecElementsSize = vecElements.size();

      for (unsigned int i = 0; i < vecElementsSize; ++i) {
         double tmp = vecElements[i];
         // for max likelihood probability have to be positive
         assert(tmp >= 0);
         sumoflogs -= ROOT::Math::Util::EvalLog(tmp);
         // std::cout << " i " << tmp << " lik " << sumoflogs << std::endl;
      }

      return sumoflogs;
   }

   /**

   !!!!!!!!!!!! to be commented

   */

   virtual double Up() const { return 0.5; }

private:
   // A pointer to the model function which describes the data
   const ParametricFunction *fModelFunction;
};

} // namespace Minuit2

} // namespace ROOT

#endif // ROOT_Minuit2_FumiliMaximumLikelihoodFCN
