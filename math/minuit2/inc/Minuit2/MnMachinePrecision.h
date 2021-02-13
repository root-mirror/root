// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Minuit2_MnMachinePrecision
#define ROOT_Minuit2_MnMachinePrecision

#include <cmath>

namespace ROOT {

namespace Minuit2 {

/**
    Sets the relative floating point (double) arithmetic precision.
    By default the precision values are obtained from the standard functions
    std::numeric_limits<double>::epsilon.
    The values can optionally be computed directly using the ComputePrecision()
    member function. For a IEEE standard floating point arithmetic the computed values and
    the one from std::numeric_limits<double>::epsilon  are the same.

    SetPrecision() method can instead be used to override Minuit's own determination,
    when the user knows that the {FCN} function Value is not calculated to
    the nominal machine accuracy.
 */

class MnMachinePrecision {

public:
   MnMachinePrecision();

   /// eps returns the smallest possible number so that 1.+eps > 1.
   double Eps() const { return fEpsMac; }

   /// eps2 returns 2*sqrt(eps)
   double Eps2() const { return fEpsMa2; }

   /// override Minuit's own determination
   void SetPrecision(double prec)
   {
      fEpsMac = prec;
      fEpsMa2 = 2. * std::sqrt(fEpsMac);
   }

   /// compute Machine precision directly instead
   /// of using default values from std::numeric_limits
   void ComputePrecision();

private:
   double fEpsMac;
   double fEpsMa2;
};

} // namespace Minuit2

} // namespace ROOT

#endif // ROOT_Minuit2_MnMachinePrecision
