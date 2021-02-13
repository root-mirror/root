// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef MN_GaussDataGen_H_
#define MN_GaussDataGen_H_

#include <vector>

namespace ROOT {

namespace Minuit2 {

class GaussDataGen {

public:
   GaussDataGen(unsigned int npar = 100);

   ~GaussDataGen() {}

   std::vector<double> Positions() const { return fPositions; }
   std::vector<double> Measurements() const { return fMeasurements; }
   std::vector<double> Variances() const { return fVariances; }

   double Sim_Mean() const { return fSimMean; }
   double Sim_var() const { return fSimVar; }
   double Sim_const() const { return 1.; }

private:
   double fSimMean;
   double fSimVar;
   std::vector<double> fPositions;
   std::vector<double> fMeasurements;
   std::vector<double> fVariances;
};

} // namespace Minuit2

} // namespace ROOT

#endif // MN_GaussDataGen_H_
