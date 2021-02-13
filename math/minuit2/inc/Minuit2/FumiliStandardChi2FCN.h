// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Minuit2_FumiliStandardChi2FCN
#define ROOT_Minuit2_FumiliStandardChi2FCN

#include "Minuit2/FumiliChi2FCN.h"
#include "Minuit2/ParametricFunction.h"
#include <cassert>
#include <vector>
#include <cmath>

namespace ROOT {

namespace Minuit2 {

/**

Class implementing the standard chi square function, which
is the sum of the squares of the figures-of-merit calculated for each measurement
point, the individual figures-of-merit being: (the Value predicted by the
model-measured Value)/standard deviation.

@author Andras Zsenei and Lorenzo Moneta, Creation date: 31 Aug 2004

@see FumiliChi2FCN

@ingroup Minuit

\todo nice formula for the documentation...

*/

class FumiliStandardChi2FCN : public FumiliChi2FCN {

public:
   /**

   Constructor which initializes chi square function for one-dimensional model function

   @param modelFCN the model function used for describing the data.

   @param meas vector containing the measured values.

   @param pos vector containing the x values corresponding to the
   measurements

   @param mvar vector containing the variances corresponding to each
   measurement (where the variance equals the standard deviation squared).
   If the variances are zero, a Value of 1 is used (as it is done in ROOT/PAW)

   */

   FumiliStandardChi2FCN(const ParametricFunction &modelFCN, const std::vector<double> &meas,
                         const std::vector<double> &pos, const std::vector<double> &mvar)
   { // this->fModelFCN = &modelFunction;
      this->SetModelFunction(modelFCN);

      assert(meas.size() == pos.size());
      assert(meas.size() == mvar.size());
      fMeasurements = meas;
      std::vector<double> x(1);
      unsigned int n = mvar.size();
      fPositions.reserve(n);
      // correct for variance == 0
      fInvErrors.resize(n);
      for (unsigned int i = 0; i < n; ++i) {
         x[0] = pos[i];
         fPositions.push_back(x);
         // PAW/ROOT hack : use 1 for 0 entries bins
         if (mvar[i] == 0)
            fInvErrors[i] = 1;
         else
            fInvErrors[i] = 1.0 / std::sqrt(mvar[i]);
      }
   }

   /**

   Constructor which initializes the multi-dimensional model function.

   @param modelFCN the model function used for describing the data.

   @param meas vector containing the measured values.

   @param pos vector containing the x values corresponding to the
   measurements

   @param mvar vector containing the variances corresponding to each
   measurement (where the variance equals the standard deviation squared).
   If the variances are zero, a Value of 1 is used (as it is done in ROOT/PAW)

   */

   FumiliStandardChi2FCN(const ParametricFunction &modelFCN, const std::vector<double> &meas,
                         const std::vector<std::vector<double>> &pos, const std::vector<double> &mvar)
   { // this->fModelFCN = &modelFunction;
      this->SetModelFunction(modelFCN);

      assert(meas.size() == pos.size());
      assert(meas.size() == mvar.size());
      fMeasurements = meas;
      fPositions = pos;
      // correct for variance == 0
      unsigned int n = mvar.size();
      fInvErrors.resize(n);
      for (unsigned int i = 0; i < n; ++i) {
         // PAW/ROOT hack : use 1 for 0 entries bins
         if (mvar[i] == 0)
            fInvErrors[i] = 1;
         else
            fInvErrors[i] = 1.0 / std::sqrt(mvar[i]);
      }
   }

   ~FumiliStandardChi2FCN() {}

   /**

   Evaluates the model function for the different measurement points and
   the Parameter values supplied, calculates a figure-of-merit for each
   measurement and returns a vector containing the result of this
   evaluation. The figure-of-merit is (Value predicted by the model
   function-measured Value)/standard deviation.

   @param par vector of Parameter values to feed to the model function.

   @return A vector containing the figures-of-merit for the model function evaluated
   for each set of measurements.

   \todo What to do when the variances are 0???!! (right now just pushes back 0...)

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
   the chi-square.

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
   std::vector<double> fMeasurements;
   // support multi dim coordinates
   std::vector<std::vector<double>> fPositions;
   std::vector<double> fInvErrors;
};

} // namespace Minuit2

} // namespace ROOT

#endif // ROOT_Minuit2_FumiliStandardChi2FCN
