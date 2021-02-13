// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#include "Minuit2/MnUserFcn.h"
#include "Minuit2/FCNBase.h"
#include "Minuit2/MnUserTransformation.h"

namespace ROOT {

namespace Minuit2 {

double MnUserFcn::operator()(const MnAlgebraicVector &v) const
{
   // call Fcn function transforming from a MnAlgebraicVector of internal values to a std::vector of external ones
   fNumCall++;

   // calling fTransform() like here was not thread safe because it was using a cached vector
   // return Fcn()( fTransform(v) );
   // make a new thread-safe implementation creating a vector each time
   // a bit slower few% in stressFit and 10% in Rosenbrock function but it is negligible in big fits

   // get first initial values of parameter (in case some one is fixed)
   std::vector<double> vpar(fTransform.InitialParValues().begin(), fTransform.InitialParValues().end());

   const std::vector<MinuitParameter> &parameters = fTransform.Parameters();
   unsigned int n = v.size();
   for (unsigned int i = 0; i < n; i++) {
      int ext = fTransform.ExtOfInt(i);
      if (parameters[ext].HasLimits()) {
         vpar[ext] = fTransform.Int2ext(i, v(i));
      } else {
         vpar[ext] = v(i);
      }
   }
   return Fcn()(vpar);
}

} // namespace Minuit2

} // namespace ROOT
