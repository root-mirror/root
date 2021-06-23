// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei, E.G.P. Bos   2003-2017

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 * Copyright (c) 2017 Patrick Bos, Netherlands eScience Center        *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Minuit2_FunctionMinimizer
#define ROOT_Minuit2_FunctionMinimizer

#include "Minuit2/MnConfig.h"
#include "Minuit2/FCNGradientBase.h"
#include <vector>

namespace ROOT {

namespace Minuit2 {

class FCNBase;
class FCNGradientBase;
class FunctionMinimum;

//_____________________________________________________________________________________
/** base class for function minimizers; user may give FCN or FCN with Gradient,
    Parameter starting values and initial Error guess (sigma) (or "step size"),
    or Parameter starting values and initial covariance matrix;
    covariance matrix is stored in Upper triangular packed storage format,
    e.g. the Elements in the array are arranged like
    {a(0,0), a(0,1), a(1,1), a(0,2), a(1,2), a(2,2), ...},
    the size is nrow*(nrow+1)/2 (see also MnUserCovariance.h);
 */

class FunctionMinimizer {

public:
   virtual ~FunctionMinimizer() {}

   // starting values for parameters and errors
   virtual FunctionMinimum Minimize(const FCNBase &, const std::vector<double> &par, const std::vector<double> &err,
                                    unsigned int strategy, unsigned int maxfcn, double toler) const = 0;

   // starting values for parameters and errors and FCN with Gradient
   virtual FunctionMinimum Minimize(const FCNGradientBase &, const std::vector<double> &par,
                                    const std::vector<double> &err, unsigned int strategy, unsigned int maxfcn,
                                    double toler) const = 0;

   // starting values for parameters and covariance matrix
   virtual FunctionMinimum Minimize(const FCNBase &, const std::vector<double> &par, unsigned int nrow,
                                    const std::vector<double> &cov, unsigned int strategy, unsigned int maxfcn,
                                    double toler) const = 0;

   // starting values for parameters and covariance matrix and FCN with Gradient
   virtual FunctionMinimum Minimize(const FCNGradientBase &, const std::vector<double> &par, unsigned int nrow,
                                    const std::vector<double> &cov, unsigned int strategy, unsigned int maxfcn,
                                    double toler) const = 0;
};

} // namespace Minuit2

} // namespace ROOT

#endif // ROOT_Minuit2_FunctionMinimizer
