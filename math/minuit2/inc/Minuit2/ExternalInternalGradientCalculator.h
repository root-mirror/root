// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei, E.G.P. Bos   2003-2017

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 * Copyright (c) 2017 Patrick Bos, Netherlands eScience Center        *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Minuit2_ExternalInternalGradientCalculator
#define ROOT_Minuit2_ExternalInternalGradientCalculator

#include "Minuit2/AnalyticalGradientCalculator.h"

namespace ROOT {

  namespace Minuit2 {

    class FCNGradientBase;
    class MnUserTransformation;

    class ExternalInternalGradientCalculator : public AnalyticalGradientCalculator {

    public:

      ExternalInternalGradientCalculator(const FCNGradientBase &fcn, const MnUserTransformation &state) :
          AnalyticalGradientCalculator(fcn, state)
      {}

      ~ExternalInternalGradientCalculator() {}

      virtual FunctionGradient operator()(const MinimumParameters &) const;

      virtual FunctionGradient operator()(const MinimumParameters &,
                                          const FunctionGradient &) const;

    };

  }  // namespace Minuit2

}  // namespace ROOT

#endif  // ROOT_Minuit2_ExternalInternalGradientCalculator
