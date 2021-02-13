// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#include "Minuit2/FCNBase.h"

#include <vector>

namespace ROOT {

namespace Minuit2 {

class Quad8F : public FCNBase {

public:
   double operator()(const std::vector<double> &par) const
   {

      double x = par[0];
      double y = par[1];
      double z = par[2];
      double w = par[3];
      double x0 = par[4];
      double y0 = par[5];
      double z0 = par[6];
      double w0 = par[7];

      return ((1. / 70.) * (21 * x * x + 20 * y * y + 19 * z * z - 14 * x * z - 20 * y * z) + w * w +
              (1. / 70.) * (21 * x0 * x0 + 20 * y0 * y0 + 19 * z0 * z0 - 14 * x0 * z0 - 20 * y0 * z0) + w0 * w0);
   }

   double Up() const { return 1.; }

private:
};

} // namespace Minuit2

} // namespace ROOT
