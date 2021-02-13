// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Minuit2_MnParameterScan
#define ROOT_Minuit2_MnParameterScan

#include "Minuit2/MnConfig.h"
#include "Minuit2/MnUserParameters.h"

#include <vector>
#include <utility>

namespace ROOT {

namespace Minuit2 {

class FCNBase;

/** Scans the values of FCN as a function of one Parameter and retains the
    best function and Parameter values found.
 */

class MnParameterScan {

public:
   MnParameterScan(const FCNBase &, const MnUserParameters &);

   MnParameterScan(const FCNBase &, const MnUserParameters &, double);

   ~MnParameterScan() {}

   // returns pairs of (x,y) points, x=parameter Value, y=function Value of FCN
   std::vector<std::pair<double, double>>
   operator()(unsigned int par, unsigned int maxsteps = 41, double low = 0., double high = 0.);

   const MnUserParameters &Parameters() const { return fParameters; }
   double Fval() const { return fAmin; }

private:
   const FCNBase &fFCN;
   MnUserParameters fParameters;
   double fAmin;
};

} // namespace Minuit2

} // namespace ROOT

#endif // ROOT_Minuit2_MnParameterScan
