// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef MN_GaussFcn_H_
#define MN_GaussFcn_H_

#include "Minuit2/FCNBase.h"

#include <vector>

namespace ROOT {

namespace Minuit2 {

class GaussFcn : public FCNBase {

public:
   GaussFcn(const std::vector<double> &meas, const std::vector<double> &pos, const std::vector<double> &mvar)
      : fMeasurements(meas), fPositions(pos), fMVariances(mvar), fErrorDef(1.)
   {
   }

   ~GaussFcn() {}

   virtual double Up() const { return fErrorDef; }
   virtual double operator()(const std::vector<double> &) const;

   std::vector<double> Measurements() const { return fMeasurements; }
   std::vector<double> Positions() const { return fPositions; }
   std::vector<double> Variances() const { return fMVariances; }

   void SetErrorDef(double def) { fErrorDef = def; }

private:
   std::vector<double> fMeasurements;
   std::vector<double> fPositions;
   std::vector<double> fMVariances;
   double fErrorDef;
};

} // namespace Minuit2

} // namespace ROOT

#endif // MN_GaussFcn_H_
