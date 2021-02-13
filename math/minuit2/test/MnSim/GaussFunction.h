// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef MN_GaussFunction_H_
#define MN_GaussFunction_H_

#include <cmath>

namespace ROOT {

namespace Minuit2 {

class GaussFunction {

public:
   GaussFunction(double mean, double sig, double constant) : fMean(mean), fSigma(sig), fConstant(constant) {}

   ~GaussFunction() {}

   double m() const { return fMean; }
   double s() const { return fSigma; }
   double c() const { return fConstant; }

   double operator()(double x) const
   {
      constexpr double two_pi = 2 * 3.14159265358979323846; // M_PI is not standard
      return c() * std::exp(-0.5 * (x - m()) * (x - m()) / (s() * s())) / (std::sqrt(two_pi) * s());
   }

private:
   double fMean;
   double fSigma;
   double fConstant;
};

} // namespace Minuit2

} // namespace ROOT

#endif // MN_GaussFunction_H_
