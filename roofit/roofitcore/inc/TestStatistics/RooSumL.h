/*
 * Project: RooFit
 * Authors:
 *   PB, Patrick Bos, Netherlands eScience Center, p.bos@esciencecenter.nl
 *
 * Copyright (c) 2016-2020, Netherlands eScience Center
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 */
#ifndef ROOT_ROOFIT_TESTSTATISTICS_RooSumL
#define ROOT_ROOFIT_TESTSTATISTICS_RooSumL

#include <TestStatistics/RooAbsL.h>
#include <TestStatistics/optional_parameter_types.h>

namespace RooFit {
namespace TestStatistics {

class RooSumL : public RooAbsL {
public:
   // main constructor
   RooSumL(RooAbsPdf* pdf, RooAbsData* data, std::vector<std::unique_ptr<RooAbsL>> components,
           RooAbsL::Extended extended = RooAbsL::Extended::Auto);
   // Note: when above ctor is called without std::moving components, you get a really obscure error. Pass as std::move(components)!

   double evaluate_partition(Section events, std::size_t components_begin,
                             std::size_t components_end) override;

   // necessary only for legacy offsetting mode in LikelihoodWrapper; TODO: remove this if legacy mode is ever removed
   std::tuple<double, double> get_subsidiary_value();

   void constOptimizeTestStatistic(RooAbsArg::ConstOpCode opcode, bool doAlsoTrackingOpt) override;

private:
   bool processEmptyDataSets() const;

   std::vector<std::unique_ptr<RooAbsL>> components_;
};

}
}

#endif // ROOT_ROOFIT_TESTSTATISTICS_RooSumL
