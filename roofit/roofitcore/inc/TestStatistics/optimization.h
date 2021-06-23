/*
 * Project: RooFit
 * Authors:
 *   PB, Patrick Bos, Netherlands eScience Center, p.bos@esciencecenter.nl
 *
 * Copyright (c) 2016-2021, Netherlands eScience Center
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 */
#ifndef ROOT_ROOFIT_TESTSTATISTICS_optimization
#define ROOT_ROOFIT_TESTSTATISTICS_optimization

// forward declarations
class RooAbsReal;
class RooArgSet;
class RooAbsData;

namespace RooFit {
namespace TestStatistics {

// this is a class only for convenience: it saves multiple friend definitions in RooAbsData for otherwise free functions
struct ConstantTermsOptimizer {
   static void enable_constant_terms_optimization(RooAbsReal *function, RooArgSet *norm_set, RooAbsData *dataset,
                                           bool applyTrackingOpt);
   static void optimize_caching(RooAbsReal *function, RooArgSet *norm_set, RooArgSet* observables, RooAbsData *dataset);
   static void disable_constant_terms_optimization(RooAbsReal *function, RooArgSet *norm_set, RooArgSet* observables, RooAbsData *dataset);
   static RooArgSet requiredExtraObservables();
};

}
}

#endif // ROOT_ROOFIT_TESTSTATISTICS_optimization
