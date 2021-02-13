// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Minuit2_FumiliStandardMaximumLikelihoodFCN
#define ROOT_Minuit2_FumiliStandardMaximumLikelihoodFCN

#include "Minuit2/FumiliMaximumLikelihoodFCN.h"
#include "Minuit2/ParametricFunction.h"
#include <vector>

namespace ROOT {

namespace Minuit2 {

/**

Class implementing the Elements member function for the standard
maximum likelihood method.

@author Andras Zsenei and Lorenzo Moneta, Creation date: 4 Sep 2004

@see FumiliMaximumLikelihoodFCN

@ingroup Minuit

*/

class FumiliStandardMaximumLikelihoodFCN : public FumiliMaximumLikelihoodFCN {

public:
   /**

   Constructor which initializes the measurement points for the one dimensional model function.

   @param modelFCN the model function used for describing the data.

   @param pos vector containing the x values corresponding to the
   measurements

   */

   FumiliStandardMaximumLikelihoodFCN(const ParametricFunction &modelFCN, const std::vector<double> &pos)
   {
      this->SetModelFunction(modelFCN);
      unsigned int n = pos.size();
      fPositions.reserve(n);
      std::vector<double> x(1);
      for (unsigned int i = 0; i < n; ++i) {
         x[0] = pos[i];
         fPositions.push_back(x);
      }
   }

   /**

   Constructor which initializes the measurement points for the multi dimensional model function.

   @param modelFCN the model function used for describing the data.

   @param pos vector containing the x values corresponding to the
   measurements

   */

   FumiliStandardMaximumLikelihoodFCN(const ParametricFunction &modelFCN, const std::vector<std::vector<double>> &pos)
   {
      this->SetModelFunction(modelFCN);
      fPositions = pos;
   }

   ~FumiliStandardMaximumLikelihoodFCN() {}

   /**

   Evaluates the model function for the different measurement points and
   the Parameter values supplied.

   @param par vector of Parameter values to feed to the model function.

   @return A vector containing the model function evaluated
   for each measurement point.

   */

   std::vector<double> Elements(const std::vector<double> &par) const;

   /**

   Accessor to the position of the measurement (x coordinate).

   @param Index Index of the measuerement the position of which to return.

   @return the position of the measurement.

   */

   virtual const std::vector<double> &GetMeasurement(int Index) const;

   /**

   Accessor to the number of measurements used for calculating
   the maximum likelihood.

   @return the number of measurements.

   */

   virtual int GetNumberOfMeasurements() const;

   /**

   Evaluate function Value, Gradient and Hessian using Fumili approximation, for values of parameters p
   The resul is cached inside and is return from the FumiliFCNBase::Value ,  FumiliFCNBase::Gradient and
   FumiliFCNBase::Hessian methods

   @param par vector of parameters

   **/

   virtual void EvaluateAll(const std::vector<double> &par);

private:
   std::vector<std::vector<double>> fPositions;
};

} // namespace Minuit2

} // namespace ROOT

#endif // ROOT_Minuit2_FumiliStandardMaximumLikelihoodFCN
