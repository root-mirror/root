// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Minuit2_VariableMetricBuilder
#define ROOT_Minuit2_VariableMetricBuilder

#include "Minuit2/MnConfig.h"
#include "Minuit2/MinimumBuilder.h"
#include "Minuit2/VariableMetricEDMEstimator.h"
#include "Minuit2/DavidonErrorUpdator.h"
#include "Minuit2/BFGSErrorUpdator.h"

#include <vector>
#include <memory>

namespace ROOT {

namespace Minuit2 {

/**
   Build (find) function minimum using the Variable Metric method (MIGRAD)
   Two possible error updators can be choosen
    - Davidon : this is the standard formula used in Migrad
    - BFGS this is the new formula based on BFGS algorithm
      (see Broyden–Fletcher–Goldfarb–Shanno algorithm
      https://en.wikipedia.org/wiki/Broyden–Fletcher–Goldfarb–Shanno_algorithm )
 */
class VariableMetricBuilder : public MinimumBuilder {

public:
   enum ErrorUpdatorType { kDavidon, kBFGS };

   VariableMetricBuilder(ErrorUpdatorType type = kDavidon) : fEstimator(VariableMetricEDMEstimator())
   {
      if (type == kBFGS)
         fErrorUpdator = std::unique_ptr<MinimumErrorUpdator>(new BFGSErrorUpdator());
      else
         fErrorUpdator = std::unique_ptr<MinimumErrorUpdator>(new DavidonErrorUpdator());
   }

   ~VariableMetricBuilder() {}

   virtual FunctionMinimum Minimum(const MnFcn &, const GradientCalculator &, const MinimumSeed &, const MnStrategy &,
                                   unsigned int, double) const;

   FunctionMinimum Minimum(const MnFcn &, const GradientCalculator &, const MinimumSeed &, std::vector<MinimumState> &,
                           unsigned int, double) const;

   const VariableMetricEDMEstimator &Estimator() const { return fEstimator; }
   const MinimumErrorUpdator &ErrorUpdator() const { return *fErrorUpdator; }

   void AddResult(std::vector<MinimumState> &result, const MinimumState &state) const;

private:
   VariableMetricEDMEstimator fEstimator;
   std::shared_ptr<MinimumErrorUpdator> fErrorUpdator;
};

} // namespace Minuit2

} // namespace ROOT

#endif // ROOT_Minuit2_VariableMetricBuilder
