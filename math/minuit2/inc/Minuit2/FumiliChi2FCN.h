// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Minuit2_FumiliChi2FCN
#define ROOT_Minuit2_FumiliChi2FCN

#include "FumiliFCNBase.h"
#include <vector>
#include "Minuit2/ParametricFunction.h"

namespace ROOT {

namespace Minuit2 {

/**

Extension of the FCNBase for the Fumili method. Fumili applies only to
minimization problems used for fitting. The method is based on a
linearization of the model function negleting second derivatives.
User needs to provide the model function. The figure-of-merit describing
the difference between the model function and the actual measurements in
the case of chi-square is the sum of the squares of the figures-of-merit
calculated for each measurement point, which is implemented by the
operator() member function. The user still has to implement the calculation
of the individual figures-of-merit (which in the majority of the cases
will be the (measured Value - the Value predicted by the model)/standard deviation
implemeted by the FumiliStandardChi2FCN;
however this form can become more complicated (see for an example Numerical Recipes'
section on "Straight-Line Data with Errors in Both Coordinates")).


@author Andras Zsenei and Lorenzo Moneta, Creation date: 24 Aug 2004

@see <A HREF="http://www.cern.ch/winkler/minuit/tutorial/mntutorial.pdf">MINUIT Tutorial</A> on function minimization,
section 5

@see FumiliStandardChi2FCN

@ingroup Minuit

*/

class FumiliChi2FCN : public FumiliFCNBase {

public:
   FumiliChi2FCN() {}

   virtual ~FumiliChi2FCN() {}

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
   @return A reference to a vector containing the values characterizing a measurement

   */

   virtual const std::vector<double> &GetMeasurement(int Index) const = 0;

   /**

   Accessor to the number of measurements used for calculating the
   present figure of merit.

   @return the number of measurements

   */

   virtual int GetNumberOfMeasurements() const = 0;

   /**

   Calculates the sum of Elements squared, ie the chi-square. The user must
   implement in a class which inherits from FumiliChi2FCN the member function
   Elements() which will supply the Elements for the sum.


   @param par vector containing the Parameter values for the model function

   @return The sum of Elements squared

   @see FumiliFCNBase#elements

   */

   double operator()(const std::vector<double> &par) const
   {

      double chiSquare = 0.0;
      std::vector<double> vecElements = Elements(par);
      unsigned int vecElementsSize = vecElements.size();

      for (unsigned int i = 0; i < vecElementsSize; ++i)
         chiSquare += vecElements[i] * vecElements[i];

      return chiSquare;
   }

   /**

   !!!!!!!!!!!! to be commented

   */

   virtual double Up() const { return 1.0; }

private:
   // A pointer to the model function which describes the data
   const ParametricFunction *fModelFunction = nullptr;
};

} // namespace Minuit2

} // namespace ROOT

#endif // ROOT_Minuit2_FumiliChi2FCN
