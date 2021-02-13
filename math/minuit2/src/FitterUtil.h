// @(#)root/minuit2:$Id$
// Author: L. Moneta    10/2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 ROOT Foundation,  CERN/PH-SFT                   *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_FitterUtil_H_
#define ROOT_FitterUtil_H_

#include <vector>

class TF1;

/// utility functions to be used in the fitter classes

namespace FitterUtil {

/**
   Evaluate integral of model function around the bin
   To use when fitting with integral option
*/
double
EvalIntegral(TF1 *func, const std::vector<double> &x1, const std::vector<double> &x2, const std::vector<double> &par);

} // namespace FitterUtil

#endif
