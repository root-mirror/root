// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef MN_GaussFcn2_H_
#define MN_GaussFcn2_H_

#include "Minuit2/FCNBase.h"

#include <vector>

namespace ROOT {

namespace Minuit2 {

class GaussFcn2 : public FCNBase {

public:
   GaussFcn2(const std::vector<double> &meas, const std::vector<double> &pos, const std::vector<double> &mvar)
      : fMeasurements(meas), fPositions(pos), fMVariances(mvar), fMin(0.)
   {
      Init();
   }
   ~GaussFcn2() {}

   virtual void Init();

   virtual double Up() const { return 1.; }
   virtual double operator()(const std::vector<double> &) const;
   virtual double ErrorDef() const { return Up(); }

   std::vector<double> Measurements() const { return fMeasurements; }
   std::vector<double> Positions() const { return fPositions; }
   std::vector<double> Variances() const { return fMVariances; }

private:
   std::vector<double> fMeasurements;
   std::vector<double> fPositions;
   std::vector<double> fMVariances;
   double fMin;
};

} // namespace Minuit2

} // namespace ROOT

#endif // MN_GaussFcn2_H_
