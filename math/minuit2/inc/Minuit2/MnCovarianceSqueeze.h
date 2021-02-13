// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Minuit2_MnCovarianceSqueeze
#define ROOT_Minuit2_MnCovarianceSqueeze

#include "Minuit2/MnMatrix.h"

namespace ROOT {

namespace Minuit2 {

class MnUserCovariance;
class MinimumError;

/**
   class to reduce the covariance matrix when a parameter is fixed by
   removing the corresponding row and index
 */
class MnCovarianceSqueeze {

public:
   MnCovarianceSqueeze() {}

   ~MnCovarianceSqueeze() {}

   MnUserCovariance operator()(const MnUserCovariance &, unsigned int) const;

   MinimumError operator()(const MinimumError &, unsigned int) const;

   MnAlgebraicSymMatrix operator()(const MnAlgebraicSymMatrix &, unsigned int) const;

private:
};

} // namespace Minuit2

} // namespace ROOT

#endif // ROOT_Minuit2_MnCovarianceSqueeze
